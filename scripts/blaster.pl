#!/usr/bin/perl

use strict;
use warnings;
$| = 1;

print "

This script performs BLAST in a batch mode with all the proteins given
in the protein list against the NCBI nr database.

Usage:
	Place input files in input/ folder in working directory, then execute:
	perl blaster.pl <working directory> <input file>

Types of input file:
	1. NCBI protein summary file containing multiple protein records.
	2. NCBI RefSeq genome record *.gb
	3. A list of accession numbers (one number per line)

Output:
	blast/<name>/<accn>.bla
	blast/<name>.log

" and exit unless @ARGV and $#ARGV;

print "Executing blaster...\n";


## public variables ##

my $i; my $j; my $s; my $t; my @a; my @b; my @c; my %h;

## modules ##

use LWP::Simple;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);
# use URI::Escape;

## global variables ##

my $wkDir = $ARGV[0];								# working directory
my $title = $ARGV[1];								# list of proteins (proteome)
$title =~ s/\.[^.]+$//;
my $query;											# current query protein
my @queries;										# list of query proteins
my @failed;											# list of failed queries
my $iRetry = 0;										# current number of retries
my $isError = 0;									# if error
my $ua = LWP::UserAgent->new;

## subroutines ##

sub blast ();
sub retry ();

## program parameters ##

my $blastMode = 0;									# BLAST mode (0: via http connection, 1: standalone BLAST, remote database, 2: standalone BLAST, local database)
my $taxaMode = 0;									# taxonomy retriever mode (0: http, 1: local, via blastdbcmd, 2: local, from source file)

my $nRetry = 10;									# maximum number of retries
my $nHit = 250;										# number of hits to return
my $maxHits = 0;									# maximum number of valid hits to preserve. if 0 then = nHit
my $evalue = 0.01;									# maximum E-value cutoff
my $percIdent = 0;									# minimum percent identity cutoff

my $blastServer = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi";
my $eSearchServer = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
my $eFetchServer = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
my $eSummaryServer = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi";

my $dbBlast = "nr";
my $exUncultured = 1;								# Exclude uncultured and environmental samples
my @searchTaxids = ();								# search under the following taxon groups (taxids)
my @ignoreTaxids = ();								# ignore organisms under the following taxids
my $eqText = "";									# entrez query parameter

my $seqBlast = 0;									# retrieve hit sequences
my $taxBlast = 0;									# retrieve taxonomy report
my $alnBlast = 0;									# retrieve multiple sequence alignment (conflicts seqBlast)

my $mergeDuplicates = 1;							# ignore hits with same taxon names and bit scores
my $taxonUCASE = 1;									# ignore taxon names that do not start with a capital letter
my @ignoreTaxa = ();								# ignore taxon names containing the following words
my @trimTaxa = ();									# trim the following words from the beginning of taxon names
													## not used in this script ##

my $blastdbcmd = "blastdbcmd";
my $blastp = "blastp";
my $nThreads = 1;									# Multiple threads for local BLAST program

## read configurations ##

if (-e "$wkDir/config.txt"){
	open IN, "<$wkDir/config.txt";
	while (<IN>){
		s/^\s+//g; s/\s+$//g; next if /^#/; next unless $_;
		$blastMode = $1 if /^blastMode=(\d)$/;
		$nHit = $1 if /^nHit=(\d+)$/;
		$evalue = $1 if /^evalue=(.+)$/;
		$percIdent = $1 if /^percIdent=(.+)$/;

		$nRetry = $1 if /^nRetry=(\d+)$/;
		$maxHits = $1 if /^maxHits=(\d+)$/;
		$dbBlast = $1 if /^dbBlast=(.+)$/;
		$eqText = $1 if /^eqText=(.+)$/;
		
		$blastServer = $1 if /^blastServer=(.+)$/;
		$blastdbcmd = $1 if /^blastdbcmd=(.+)$/;
		$blastp = $1 if /^blastp=(.+)$/;
		$nThreads = $1 if /^nThreads=(\d+)$/;
		
		$exUncultured = $1 if /^exUncultured=([01])$/;
		@searchTaxids = split(/,/, $1) if /^searchTaxids=(.+)$/;
		@ignoreTaxids = split(/,/, $1) if /^ignoreTaxids=(.+)$/;
		push @ignoreTaxa, split(/,/, $1) if /^ignoreTaxa=(.+)$/;
		$taxonUCASE = $1 if /^taxonUCASE=([01])$/;
		$mergeDuplicates = $1 if /^mergeDuplicates=([01])$/;
		@trimTaxa = split(/,/, $1) if /^trimTaxa=(.+)$/;

		$seqBlast = $1 if /^seqBlast=([01])$/;		
		$taxBlast = $1 if /^taxBlast=([01])$/;
		$alnBlast = $1 if /^alnBlast=([01])$/;
	}
	close IN;
}


## generate Entrez query text ##

unless ($eqText){
	if ($exUncultured){
		$eqText = "all [filter] NOT(environmental samples[filter] OR metagenomes[orgn])";
	}
	if (@searchTaxids){
		for ($i=0; $i<=$#searchTaxids; $i++){
			next unless ($searchTaxids[$i] =~ /^\d+$/);
			$eqText .= " OR " if $eqText;
			$eqText .= "txid$searchTaxids[$i]\[orgn\]";
		}
	}
	if (@ignoreTaxids){
		$eqText .= " NOT (";
		for ($i=0; $i<=$#ignoreTaxids; $i++){
			next unless $ignoreTaxids[$i] =~ /^\d+$/;
			$eqText .= " OR " if $i;
			$eqText .= "txid$ignoreTaxids[$i]\[orgn\]";
		}
		$eqText .= ")";
	}
}

## read query list ##

print "Reading protein list $title... ";
open QUERY, "<$wkDir/input/$ARGV[1]" or die "Error: input file $ARGV[1] not accessible.\n";
if ($ARGV[1] =~ /\.gbk?$/i){							# GenBank file
	while(<QUERY>){
		s/\s+$//;
		push @queries, $1 if /^\s+\/protein_id="(.+)"$/;
	}
}else{													# protein summary
	while(<QUERY>){
		s/\s+$//; next unless $_; next if /^\d/;
		@a = split /\s/;
		push @queries, $a[0];
	}
}
close QUERY;
die "Error: no protein found.\n" unless @queries;
print "done. ". scalar (@queries) ." protein(s) read.\n";

mkdir "$wkDir/blast" unless -d "$wkDir/blast";
mkdir "$wkDir/blast/$title" unless -d "$wkDir/blast/$title";

@a = localtime(time);
$s = (++$a[4])."/$a[3]/".($a[5]+=1900)." $a[2]:$a[1]:$a[0]";
open LOG, ">>$wkDir/blast/$title.log";
print LOG "Program started at $s.\n";
print LOG "Number of queries: ". scalar (@queries) .".\n";
close LOG;


## perform batch BLAST ##

foreach (@queries){
	$query = $_;
	$query =~ s/\.[\d]+$//;
	next if -e "$wkDir/blast/$title/$query.bla";
	print "Blasting $query ...";
	blast;
	print " done.\n" unless $isError;
	print " error.\n" if $isError;
}

@a = localtime(time);
$s = (++$a[4])."/$a[3]/".($a[5]+=1900)." $a[2]:$a[1]:$a[0]";
open LOG, ">>$wkDir/blast/$title.log";
print LOG "\nProgram ended at $s.\n";
close LOG;

exit 0;


## perform single BLAST ##

sub blast (){
	if ($blastMode){

	## BLAST using standalone ncbi-blast+ program ##

		my @hits;
		# the report contains: accession, gi, length, taxid, sequence, title
		my @out = `$blastdbcmd -db $dbBlast -entry $query -outfmt \"%a %g %l %T %s %t\"`;
		foreach (@out){
			s/\s+$//; @b = split (/\s+/);
			last if ($b[0] eq $query or $b[0] =~/^$query\.\d+$/);
		}
		open (OUT, ">$wkDir/blast/$title/$query.bla");
		print OUT "#NEXUS\nBEGIN QUERY;\n";
		print OUT "\tGI=$b[1];\n\tAccession=$b[0];\n\tLength=$b[2];\n";
		my $length = $b[2];
		$s = join (" ", @b[5..$#b]);
		$s =~ /^\s*(.+\S)\s*\[(.+)\]$/;
		print OUT "\tProduct=$1;\n\tOrganism=$2;\nEND;\n\n";
		open TMP, ">blast.seq"; print TMP $b[4]; close TMP;
		# the report contains: subject accessions (all), E-value, bit score, aligned part of subject sequence
		$s = "$blastp -query blast.seq -db $dbBlast -evalue $evalue -max_target_seqs $nHit -outfmt \"6 sallacc evalue bitscore pident sseq\"";
		if ($blastMode == 1){
			$s .= " -remote";
			$s .= " -entrez_query $eqText" if $eqText;
		}else{
			$s .= " -num_threads $nThreads" if ($nThreads > 1);
		}
		@out = `$s`;
		foreach (@out){
			s/\s+$//; @a = split (/\t/);
			next if ($percIdent and $percIdent > $a[3]); # % identity cutoff
			my %hit = ();
			$hit{'expect'} = $a[1];
			$hit{'score'} = $a[2];
			$hit{'identity'} = $a[3];
			$hit{'sequence'} = $a[4];
			my @accns = split (/;/, $a[0]);
			my %hittaxa = (); # accn -> organism
			foreach my $accn (@accns){
				my @out2 = `$blastdbcmd -db $dbBlast -entry $accn -outfmt \"%a %T %t\"`;
				foreach (@out2){
					s/\s+$//; @b = split (/\s+/);
					last if ($b[0] eq $accn or $b[0] =~/^$accn\.\d+$/);
				}
				%h = ('taxid', $b[1], 'organism', "");
				$s = join (" ", @b[1..$#b]);
				$h{'organism'} = $1 if $s =~ /\[(.+)\]$/;
				$hittaxa{$accn} = {%h};
			}
			my %taxahit = (); # organism -> hits
			foreach my $accn (keys %hittaxa){
				if (exists $taxahit{$hittaxa{$accn}{'taxid'}}){
					$taxahit{$hittaxa{$accn}{'taxid'}}{'accns'} .= "/$accn";
				}else{
					%h = ('organism', $hittaxa{$accn}{'organism'}, 'accns', $accn);
					$taxahit{$hittaxa{$accn}{'taxid'}} = {%h};
				}
			}
			foreach my $taxid (keys %taxahit){
				$hit{'organism'} = $taxahit{$taxid}{'organism'};
				$hit{'accn'} = $taxahit{$taxid}{'accns'};
				$hit{'taxid'} = $taxid;
				if (@ignoreTaxa){
					$i = 0;
					foreach (@ignoreTaxa){
						if ($hit{'organism'} =~ /$_/){ $i = 1; last; }
					}
					next if $i;
				}
				push @hits, {%hit};
			}
			# if (@accns > 10){ for ($i=0; $i<=$#hits; $i++){ print "$hits[$i]{'accn'}\t$hits[$i]{'organism'}\n"; } die; }
			last if scalar @hits >= $nHit;
		}
		print OUT "BEGIN ORGANISM;\n";
		for ($i=0; $i<=$#hits; $i++){
			print OUT $hits[$i]{'accn'}."\t".$hits[$i]{'organism'}."\t".$hits[$i]{'identity'}."\t".$hits[$i]{'taxid'}."\t".$hits[$i]{'score'}."\t".$hits[$i]{'expect'}."\n";
		}
		print OUT ";\nEND;\n\n";
		print OUT "BEGIN DATA;\n";
		print OUT "\tDIMENSIONS NTAX=". scalar (@hits) ." NCHAR=$length;\n\tFORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\n\tMATRIX\n";
		for ($i=0; $i<=$#hits; $i++){
			@a = split (/\//, $hits[$i]{'accn'});
			print OUT $a[0]."\t".$hits[$i]{'sequence'}."\n";
			last if $i > $nHit;
		}		
		print OUT ";\nEND;\n\n";
		close OUT;
		
	}else{

	## BLAST via http connection to NCBI server ##

		# examples of commands:
		# search:
		## http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&PROGRAM=blastp&DATABASE=nr&QUERY=CBI80307&EXPECT=1e-5&FILTER=m S&MAX_NUM_SEQ=200
		# search info:
		## http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID=K2NW0S8U016&FORMAT_OBJECT=SearchInfo
		# tabular report:
		## http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID=K2NW0S8U016&ALIGNMENT_VIEW=Tabular&FORMAT_TYPE=Text&ALIGNMENTS=200
		# taxonomy report:
		## http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID=K2NW0S8U016&FORMAT_OBJECT=TaxBlast&ALIGNMENTS=100000
		# multiple sequence alignment:
		## http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID=K2NW0S8U016&ALIGNMENT_VIEW=FlatQueryAnchoredNoIdentities&FORMAT_TYPE=Text&ALIGNMENTS=100000
		# search protein for GI:
		## http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=CBI80307,CBI80501
		# get summary
		## http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=319406667,319406866
		# get sequence
		## http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=FASTA&id=319406667,319406866

		$isError = 0;
		my $args = "CMD=Put&PROGRAM=blastp&DATABASE=$dbBlast&QUERY=$query&EXPECT=$evalue&FILTER=m S";
		$args .= "&MAX_NUM_SEQ=$nHit" if ($nHit - 100);
		$args .= "&EQ_TEXT=$eqText" if $eqText;
		my $req = new HTTP::Request POST => $blastServer;
		$req->content_type('application/x-www-form-urlencoded');
		$req->content($args);
		my $response = $ua->request($req);
		my $rid;
		if ($response->content =~ /^    RID = (.*$)/m){ $rid = $1; }else{ retry; return; };
		if ($response->content =~ /^    RTOE = (.*$)/m){ sleep $1; }else{ retry; return; };
		while (1){
			sleep 5;
			$args = "$blastServer?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
			$req = new HTTP::Request GET => $args;
			$response = $ua->request($req);
			$s = $response->content;
			if ($s =~ /\s+Status=WAITING/m){ next; }
			if ($s =~ /\s+Status=FAILED/m){ $isError = 1; last; }
			if ($s =~ /\s+Status=UNKNOWN/m){ $isError = 1; last; }
			if ($s =~ /\s+Status=READY/m){
				if ($s =~ /\s+ThereAreHits=yes/m){ last; }
				else{ last; } # no hits;
			}
			$isError = 1;
			last;
		}
		if ($isError){ print "error"; retry; return; }
		
		# retrieve tabular report #
		
		$args = "$blastServer?CMD=Get&ALIGNMENT_VIEW=Tabular&FORMAT_TYPE=Text&RID=$rid";
		$args .= "&ALIGNMENTS=$nHit" if ($nHit - 100);
		$req = new HTTP::Request GET => $args;
		$response = $ua->request($req);
		$s = $response->content;
		if (($s !~ /# blastp/) or ($s !~ /# Query:\s/)){ print "tab"; retry; return; }
		my %hits = ();
		my $hitid = 0;
		my ($gi, $accession, $product, $organism, $length);
		my ($expect, $score, $identity);
		foreach (split(/\n/, $s)){
			if (/^# Query/){
				s/>.*//;
				@a = split (/\|/);
				$gi = $a[1];
				$accession = $a[3];
				$a[4] =~ /^\s*(.+\S)\s*\[(.+)\]/;
				$product = $1;
				$organism = $2;
			}
			next unless (/^gi/);
			@a = split (/\t/);
			next if ($percIdent and $percIdent > $a[2]); # % identity cutoff
			# next unless ($#a);
			@b = split (/;/, $a[1]);
			$hitid ++;
			foreach (@b){
				next unless (/\|$/);
				@c = split (/\|/);
				$c[3] =~ s/\.\d+$//;
				my %hit = ('id', $hitid, 'accn', $c[3], 'expect', $a[11], 'score', $a[12], 'identity', $a[2]);
				$hits{$c[1]} = {%hit};
			}
		}

		# retrieve taxonomy information #

		$i = 0; # count
		@a = (); # all results
		@b = (); # subset of GIs
		foreach (keys %hits){
			$i ++;
			push (@b, $_);
			if ($i == 190){ # a URI should not exceed ~2000 characters, which is approximately 190 GIs.
				while (1){
					$s = get $eSummaryServer."?db=protein&id=".join (",", @b);
					last if (defined $s);
					sleep 10;
				}
				push (@a, $1) while ($s =~ s/<DocSum>(.+?)<\/DocSum>//s);
				$i = 0; @b = ();
			}
		}
		if (@b){
			while (1){
				$s = get $eSummaryServer."?db=protein&id=".join (",", @b);
				last if (defined $s);
				sleep 10;
			}
			push (@a, $1) while ($s =~ s/<DocSum>(.+?)<\/DocSum>//s);	
		}

		foreach (@a){
			/<Id>(\d+)<\/Id>/;
			$s = $1;
			/<Item Name=\"TaxId\" Type=\"Integer\">(\d+?)<\/Item>/;
			$hits{$s}{'taxid'} = $1;
			/<Item Name=\"Title\" Type=\"String\">.*\[([^\[\]]+?)\]<\/Item>/;
			$hits{$s}{'organism'} = $1;
			/<Item Name=\"Length\" Type=\"Integer\">(\d+)<\/Item>/;
			$hits{$s}{'length'} = $1;
		}
		
		# read query information # (??????????)
		
		foreach (keys %hits){
			if ($query eq $hits{$_}{'accn'}){
				$length = $hits{$_}{'length'};
				last;
			}
		}
		unless ($length){
			while (1){
				$s = get "$eSearchServer?db=protein&term=$query";
				last if (defined $s);
				sleep 10;
			}
			$s =~ /<Id>(\d+)<\/Id>/;
			while (1){
				$s = get "$eSummaryServer?db=protein&id=$1";
				last if (defined $s);
				sleep 10;
			}
			$s =~ /<Item Name=\"Length\" Type=\"Integer\">(\d+)<\/Item>/;
			$length = $1;
		}

		# discard hits whose taxonomy information is unidentified #

		foreach (keys %hits){
			delete $hits{$_} unless exists $hits{$_}{'taxid'};
			delete $hits{$_} unless $hits{$_}{'taxid'};
			delete $hits{$_} unless exists $hits{$_}{'organism'};
			delete $hits{$_} unless $hits{$_}{'organism'};
			delete $hits{$_} if ($hits{$_}{'organism'} =~ /^Unresolved/);
			delete $hits{$_} if ($taxonUCASE and ($hits{$_}{'organism'} !~ /^[A-Z]/));
		}

		# discard hits whose organism names contain designated strings #
		
		if (@ignoreTaxa){
			foreach (keys %hits){
				foreach $s (@ignoreTaxa){
					if ($hits{$_}{'organism'} =~ /$s/){
						delete $hits{$_};
						last;
					}
				}
			}
		}

		# merge identical proteins (defined by NCBI) from one organism #

		$i = 0; # current hit ID
		$j = 0; # current taxid
		$s = ""; # current key
		foreach ( sort {$hits{$a}{'id'} <=> $hits{$b}{'id'}} keys %hits){
			if ($hits{$_}{'id'} != $i){
				$i = $hits{$_}{'id'};
				$j = $hits{$_}{'taxid'};
				$s = $_;
			}else{
				if ($hits{$_}{'taxid'} == $j){
					$hits{$s}{'accn'} .= "/".$hits{$_}{'accn'};
					delete $hits{$_};
				}
			}
		}

		# merge duplicated hits (proteins with same bit score from one organism) #
		
		if ($mergeDuplicates){
			$i = 0; # current bit score
			$j = 0; # current taxid
			$s = ""; # current key
			foreach ( sort {($hits{$b}{'score'} <=> $hits{$a}{'score'}) or ($hits{$a}{'taxid'} <=> $hits{$b}{'taxid'})} keys %hits){
				if (($hits{$_}{'score'} != $i) or ($hits{$_}{'taxid'} != $j)){
					$i = $hits{$_}{'score'};
					$j = $hits{$_}{'taxid'};
					$s = $_;
				}else{
					$hits{$s}{'accn'} .= "/".$hits{$_}{'accn'};
					delete $hits{$_};
				}
			}
		}

		# remove excessive hits #
		
		if ($maxHits){
			$i = 0;
			foreach (sort {$hits{$b}{'score'} <=> $hits{$a}{'score'}} keys %hits){
				$i ++;
				delete $hits{$_} if ($i > $maxHits);
			}
		}

		# reorder accession numbers #
		# order: NP_ > XP_ = YP_ = ZP_ = AP_ > all else
		
		foreach (keys %hits){
			next unless $hits{$_}{'accn'} =~ /\//;
			@a = split (/\//, $hits{$_}{'accn'});
			@b = ();
			foreach (@a){
				if (/^.P_/){ unshift (@b, $_); }
				else{ push (@b, $_); }
			}
			$hits{$_}{'accn'} = join ("/", @b);
		}

		# create output file #
		
		open (OUT, ">$wkDir/blast/$title/$query.bla");
		print OUT "#NEXUS\nBEGIN QUERY;\n";
		print OUT "\tGI=$gi;\n\tAccession=$accession;\n\tLength=$length;\n\tProduct=$product;\n\tOrganism=$organism;\nEND;\n\n";
		
		# retrieve taxonomy report (using TaxBLAST)
		
		if ($taxBlast){
			$args = "$blastServer?CMD=Get&FORMAT_TYPE=HTML&FORMAT_OBJECT=TaxBlast&RID=$rid&ALIGNMENTS=100000";
			$req = new HTTP::Request GET => $args;
			$response = $ua->request($req);
			$s = $response->content;
			if (($s !~ /Tax BLAST Report/) or ($s !~ /Lineage Report/)){
				print OUT "BEGIN ERROR;\nRetrieval of taxonomy report failed.\n;\nEND;\n\n";
			}else{
				$t = 0; # reading status
				foreach (split(/\n/, $s)){
					if ($_ eq "<B><A NAME=lineage>Lineage Report</A></B><BR>"){
						print OUT "BEGIN LINEAGE;\n";
						$t = 1; next;
					}
					if (($_ eq "</FONT></PRE><HR>") and $t){
						print OUT ";\nEND;\n\n";
						$t = 0; next;
					}
					if ($t){
						s/\<.*?>//g;
						s/( hit[ s] \[.*?\])(.*)$/$1/;
						if (@ignoreTaxa){
							$i = 0;
							foreach $s (@ignoreTaxa){
								if (/$s/){ $i = 1; last; }
							}
							next if $i;
						}
						print OUT "$_\n";
						next;
					}
				}
			}
		}

		# output hit table #
		
		print OUT "BEGIN ORGANISM;\n";
		foreach (sort {($hits{$b}{'score'} <=> $hits{$a}{'score'}) or ($hits{$a}{'id'} <=> $hits{$b}{'id'})} keys %hits){
			print OUT $hits{$_}{'accn'}."\t".$hits{$_}{'organism'}."\t".$hits{$_}{'identity'}."\t".$hits{$_}{'taxid'}."\t".$hits{$_}{'score'}."\t".$hits{$_}{'expect'}."\n";
		}
		print OUT ";\nEND;\n\n";

		# retrieve hit sequences #
		
		if ($seqBlast){
			$i = 0; # count
			$s = ""; # all results
			@b = (); # subset of GIs
			foreach (keys %hits){
				$i ++;
				push (@b, $_);
				if ($i == 190){ # a URI should not exceed ~2000 characters, which is approximately 190 GIs.
					while (1){
						$t = get $eFetchServer."?db=protein&rettype=FASTA&id=".join (",", @b);
						last if (defined $t);
						sleep 10;
					}
					$s .= $t;
					$i = 0; @b = ();
				}
			}
			if (@b){
				while (1){
					$t = get $eFetchServer."?db=protein&rettype=FASTA&id=".join (",", @b);
					last if (defined $t);
					sleep 10;
				}
				$s .= $t;
			}
			$i = 0; # current GI
			foreach (split (/\n/, $s)){
				next unless $_;
				if (/^>gi\|(\d+)\|/){
					$i = $1;
					$hits{$i}{'sequence'} = "";
				}else{
					$hits{$i}{'sequence'} .= $_ if $i;
				}
			}
		}

		# retrieve multiple sequence alignment (conflicts seqBlast)

		if ($alnBlast and not $seqBlast){
			$args = "$blastServer?CMD=Get&ALIGNMENT_VIEW=FlatQueryAnchoredNoIdentities&FORMAT_TYPE=Text&RID=$rid&ALIGNMENTS=100000";
			$req = new HTTP::Request GET => $args;
			$response = $ua->request($req);
			$s = $response->content;
			if (($s !~ /blastp/i) or ($s !~ /\nQuery=\s/)){
				print OUT "BEGIN ERROR;\nRetrieval of multiple sequence alignment failed.\n;\nEND;\n\n";
			}else{
				@c = split(/\n/, $s);
				$t = 0; # reading status
				my $iBlock = -1; # block ID
				my %seqs; # sequence alignment
				foreach (@c){
					# Read alignment
					if ($_ eq "ALIGNMENTS"){ $t = 1; next; }
					if ($t and /^\s/){ $t = 0; last; }
					next unless $t;
					unless ($_){ $iBlock ++; next; } # Start a new block if empty line
					@a = split(/\s+/,substr($_,0,19)); # id
					#next if ($a[0] eq "Query");
					$_ =~ /(\s\s\d*$)/;
					$s = substr($_,19,length($_)-19-length($1)); # sequence
					$s =~ s/\s/-/g;
					$seqs{$a[0]} = "-"x($iBlock*60) unless (exists $seqs{$a[0]}); # add new sequence
					$seqs{$a[0]} .= $s; # add new sequence
				}
				$i = 0;
				foreach $s(keys %seqs){
					$i = length($seqs{$s}) if (length($seqs{$s}) > $i);
				}
				foreach $s(keys %seqs){
					$seqs{$s} .= "-"x($i - length($seqs{$s})) if (length($seqs{$s}) < $i);
				}
				foreach $t (keys %hits){
					@a = split /\//, $hits{$t}{'accn'};
					foreach (@a){
						if (exists $seqs{$_}){
							$hits{$t}{'sequence'} = $seqs{$_};
							last;
						}
					}
				}
			}
		}

		# output sequences #

		if ($seqBlast or $alnBlast){
			$i = 0; # number of sequences
			$j = 0; # maximum length of sequence
			foreach (keys %hits){
				next unless exists $hits{$_}{'sequence'};
				$i ++;
				$j = length ($hits{$_}{'sequence'}) if (length ($hits{$_}{'sequence'}) > $j);
			}
			print OUT "BEGIN DATA;\n";
			print OUT "\tDIMENSIONS NTAX=$i NCHAR=$j;\n\tFORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\n\tMATRIX\n";
			foreach (sort {($hits{$a}{'id'} <=> $hits{$b}{'id'})} keys %hits){
				next unless exists $hits{$_}{'sequence'};
				@a = split (/\//, $hits{$_}{'accn'});
				print OUT $a[0]."\t".$hits{$_}{'sequence'}."\n";
			}		
			print OUT ";\nEND;\n\n";
		}

		close OUT;

	}
	open LOG, ">>$wkDir/blast/$title.log";
	print LOG "$query\t.\n";
	close LOG;
}

## retry BLAST ##

sub retry (){
	if ($iRetry < $nRetry){ # retry
		print ".";
		$iRetry ++;
		sleep 10;
		blast;
	}else{ # fail
		$iRetry = 0;
		push @failed, $query;
		open LOG, ">>$wkDir/blast/$title.log";
		print LOG "$query\tfailed.\n";
		close LOG;
	}
}


