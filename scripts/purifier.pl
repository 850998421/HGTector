#!/usr/bin/perl

use warnings;
use strict;
$| = 1;

print "

This script performs post-BLAST treatments, including: removal of invalid or redundant
hits, multiple sequence alignment, alignment clean-up, tree reconstruction, and distance
matrix generation.

Usage:
	purifier.pl <working directory>

Output:
	append to each blast file (*.bla)

" and exit unless @ARGV;

print "Executing purifier...\n";


## global variables ##

my $i; my $j; my $n; my $s; my $t; my @a; my @b; my @c; my %h;

my @sets;											# protein sets (genomes)
my @hits;											# blast hits, each element is a hash, including:
														# accn, accns, organism, group, taxid, score, distance, ignore

my %self;											# query information
my %seqs;											# accn -> sequence
my %names;											# accn -> organism name
my $nSeq;											# number of hits with sequences
my $nChar;											# number of aa per sequence in alignment
my $tree;											# phylogenetic tree with annotations


## program parameters ##

my $wkDir = $ARGV[0];							# working directory

my $minHits = 3;									# minimal number of valid hits to consider for HGT
my $minSize = 0;									# minimal size (aa) of a valid protein (0 means infinite)

my @ignoreTaxa = ();								# ignore taxon names containing the following words
my $taxonUCASE = 0;								# ignore taxon names that do not start with a capital letter
my $ignoreParalogs = 1;							# ignore potential paralogs (hits with same taxon names but different bit scores)
my $ignoreSeqRepeats = 1;						# ignore repeated sequences (hits targetting different regions of the same protein)
my $ignoreSubspecies = 0;						# ignore more than one subspecies from the same species

my $trimSeq = 0;									# trim unreliable regions from sequence alignment (options: 0, gblocks)
my $realign = 0;									# realign sequences (options: 0, clustalw, mafft, muscle, clustalo)
my $buildTree = 0;								# Build phylogenetic tree using specified program (options: 0, clustalw, mafft, phyml, raxml)
my $bsTree = 0;									# perform bootstrap for designated times.
my $distance = 0;									# use distance matrix instead of BLAST bit scores to rank hits (options: 0, clustalw, mafft, raxml)
my $aaModel = "WAG";								# protein substitution model for RAxML (options: LG, WAG)

my $gblocks = "Gblocks";
my $clustalw = "clustalw";
my $mafft = "mafft";
my $muscle = "muscle";
my $raxml = "raxmlHPC";
my $phyml = "phyml";
my $fasttree = "FastTree";


## public subroutines ##

sub read_fasta ($);
sub write_fasta ($);
sub write_phylip ($);

## read configurations ##

if (-e "$wkDir/config.txt"){
	open IN, "<$wkDir/config.txt";
	my $readMonitor = 0; my $readMiddle = 0;
	while (<IN>){
		s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_;
		$minHits = $1 if /^minHits=(\d+)$/;
		$minSize = $1 if /^minSize=(\d+)$/;
		@ignoreTaxa = split(/,/, $1) if /^ignoreTaxa=(.+)$/;
		$taxonUCASE = $1 if /^taxonUCASE=([01])$/;
		$ignoreParalogs = $1 if /^ignoreParalogs=([01])$/;
		$ignoreSeqRepeats = $1 if /^ignoreSeqRepeats=([01])$/;
		$ignoreSubspecies = $1 if /^ignoreSubspecies=([01])$/;
		$trimSeq = $1 if /^trimSeq=(.+)$/;
		$realign = $1 if /^realign=(.+)$/;
		$buildTree = $1 if /^buildTree=(.+)$/;
		$bsTree = $1 if /^bsTree=(.+)$/;
		$distance = $1 if /^distance=(.+)$/;
		$aaModel = $1 if /^aaModel=(.+)$/;
		$gblocks = $1 if /^gblocks=(.+)$/;
		$clustalw = $1 if /^clustalw=(.+)$/;
		$mafft = $1 if /^mafft=(.+)$/;
		$muscle = $1 if /^muscle=(.+)$/;
		$raxml = $1 if /^raxml=(.+)$/;
		$phyml = $1 if /^phyml=(.+)$/;
		$fasttree = $1 if /^fasttree=(.+)$/;
	}
	close IN;
}


## read protein sets ##

print "\nReading protein sets...";
opendir (DIR, "$wkDir/blast") or die "Blast folder not exist.\n";
@a = readdir(DIR);
close DIR;

foreach (@a){
	next if (/^\./);
	push @sets, $_ if -d "$wkDir/blast/$_";
}
print "done. ";
die "No protein sets detected.\n" unless @sets;
print @sets." protein sets detected.\n";

print "\n0-------------25-------------50------------75------------100%\n";

foreach my $set (@sets){
	my $dir = "$wkDir/blast/$set";
	opendir (DIR, "$dir") or next;
	my @files = grep (/\.bla$/, readdir(DIR));
	close DIR;
	next unless @files;
	
	my $iProtein = 0;
	my $iProgress = 0;
	my $nProtein = $#files+1;
	print "$set has $nProtein proteins. Purifying...\n";
	
	foreach my $file (@files){
		$iProtein ++;
		@hits = (); %seqs = (); %self = (); $nSeq = 0; $nChar = 0; $tree = ""; %names = ();
		open IN, "<$dir/$file" or next;
		my $reading = 0;
		while (<IN>) {
			s/\s+$//g;
			if (/^BEGIN QUERY/){ $reading = "query"; next; }
			if (/^BEGIN ORGANISM/){ $reading = "organism"; next; }
			if (/^BEGIN DATA/){ $reading = "data"; next; }
			if (/^END;/){ $reading = 0; next; }
			if ($reading eq "query"){
				if (/^\tAccession=(.+);$/){ $self{'accn'} = $1; $self{'accn'} =~ s/\.[\d]+$//; }
				$self{'length'} = $1 if /^\tLength=(.+);$/;
				$self{'product'} = $1 if /^\tProduct=(.+);$/;
			}
			if ($reading eq "organism"){ # read organisms
				last unless exists $self{'accn'};
				next if /^;/;
				@a = split (/\t/);
				my %hit = ();
				$hit{'organism'} = $a[1];
				$hit{'group'} = $a[2];
				$hit{'taxid'} = $a[3];
				$hit{'score'} = $a[4];
				$hit{'expect'} = $a[5];
				$hit{'accns'} = $a[0];
				@a = split(/\//, $a[0]);
				$hit{'accn'} = $a[0];
				push @hits, {%hit};
				unless (exists $self{'taxid'}){ # identify self
					foreach (@a){
						if ($self{'accn'} eq $_){
							$self{'id'} = $#hits;
							$self{'taxid'} = $hit{'taxid'};
							$self{'score'} = $hit{'score'};
							$self{'expect'} = $hit{'expect'};
							$self{'organism'} = $hit{'organism'};
							last;
						}
					}
				}
			}
			if ($reading eq "data"){ # read sequences
				next if /^;/;
				next if /^\t/;
				@a = split (/\t/);
				$seqs{$a[0]} = $a[1];
				$nSeq ++;
				$nChar = length($a[1]) unless $nChar;
			}
		}
		close IN;
		print "Blast report of $file is incomplete.\n" and next unless exists $self{'accn'};
		next if @hits < $minHits;


		## ignore lower lettered taxon names ##
		
		if ($taxonUCASE){
			for ($i=0; $i<=$#hits; $i++){
				if ($hits[$i]{'organism'} !~ /^[A-Z]/){
					$hits[$i]{'ignore'} = 1;
				}
			}
		}

		## ignore user-defined taxon names ##
		
		if (@ignoreTaxa){
			for ($i=0; $i<=$#hits; $i++){
				$j = 0;
				foreach (@ignoreTaxa){
					if ($hits[$i]{'organism'} =~ /$_/){ $j = 1; last; }
				}
				$hits[$i]{'ignore'} = 1 if $j;
			}
		}

		## ignore paralogs ##
		
		if ($ignoreParalogs){
			my %usedTaxids = ();
			for ($i=0; $i<=$#hits; $i++){
				next if exists $hits[$i]{'ignore'};
				if (exists $usedTaxids{$hits[$i]{'taxid'}}){
					$hits[$i]{'ignore'} = 1;
				}else{
					$usedTaxids{$hits[$i]{'taxid'}} = 1;
				}
			}
		}

		## ignore repeated sequences ##
		
		if ($ignoreSeqRepeats){
			for ($i=1; $i<=$#hits; $i++){
				next if exists $hits[$i]{'ignore'};
				if ($hits[$i]{'accn'} eq $hits[$i-1]{'accn'}){
					$hits[$i]{'ignore'} = 1;
				}
			}
		}

		## ignore subspecies ##
		
		if ($ignoreSubspecies and -e "$wkDir/taxonomy/taxa.db"){
			my %speciesDB = ();
			open IN, "<$wkDir/taxonomy/taxa.db";
			while (<IN>){
				@a = split (/\t/);
				$speciesDB{$a[0]} = $a[3];
			}
			close IN;
			my %usedSpecies = ();
			for ($i=0; $i<=$#hits; $i++){
				next if exists $hits[$i]{'ignore'};
				if (exists $speciesDB{$hits[$i]{'taxid'}} and $speciesDB{$hits[$i]{'taxid'}}){
					$s = $speciesDB{$hits[$i]{'taxid'}};
					if (exists $usedSpecies{$s}){
						$hits[$i]{'ignore'} = 1;
					}else{
						$usedSpecies{$s} = 1;
					}
				}else{
					$hits[$i]{'ignore'} = 1;
				}
			}
		}
		
		## realign with ClustalW ##
		
		if ($realign eq "clustalw"){
			next if write_fasta ("$dir/temp") < 2;
			system "$clustalw -infile=$dir/temp -quicktree -output=fasta -quiet > $dir/buffer";
			die "Execution of ClustalW failed. Please check.\n" unless -s "$dir/temp.fasta";
			read_fasta "$dir/temp.fasta";
			unlink "$dir/temp", "$dir/temp.fasta", "$dir/temp.dnd", "$dir/buffer";
		}

		## trim alignment with GBlocks ##
		
		if ($trimSeq eq "gblocks"){
			write_fasta "$dir/temp";
			system "$gblocks $dir/temp -b2=".(int($nSeq/2)+1)." -b3=3 -b4=6 -b5=a > $dir/buffer";
			die "Execution of Gblocks failed. Please check.\n" unless -s "$dir/temp-gb";
			read_fasta "$dir/temp-gb";
			foreach (keys %seqs){
				$seqs{$_} =~ s/\s//g;
				delete $seqs{$_} unless ($seqs{$_} =~ /[a-zA-Z]/);
			}
			unlink "$dir/temp", "$dir/temp-gb", "$dir/temp-gb.htm";
		}

		## create translation table for tree annotation ##
		
		if ($buildTree){
			for ($i=0; $i<=$#hits; $i++){
				$s = $hits[$i]{'organism'};
				$s =~ s/[^a-zA-Z0-9,\.\-]/ /g;
				$names{$hits[$i]{'accn'}} = $s;
			}
		}

		## build phylogenetic tree and compute distance matrix with ClustalW (Neighbor-Joining) ##
	
		if ($buildTree eq "clustalw"){
			if (write_fasta("$dir/temp")>= 3){
				system "$clustalw -infile=$dir/temp -tree -outputtree=dist -quiet > $dir/buffer";
				die "Execution of ClustalW failed. Please check.\n" unless -s "$dir/temp.ph";
				open IN, "<$dir/temp.ph";
				while (<IN>){
					s/\s+$//;
					$tree .= $_;
				}
				close IN;
				$tree =~ s/([,\(])$_:/$1'$names{$_}':/ for (keys %names);
				if ($distance){
					$s = "";
					open IN, "<$dir/temp.dst";
					$_ = <IN>; $s = <IN>; $s =~ s/\s+$//;
					while (<IN>){
						last unless (/^\s/);
						s/\s+$//; s/^\s+//;
						$s .= " $_";
					}
					close IN;
					@a = split (/\s+/, $s);
					shift @a;
					for ($i=0; $i<=$#hits; $i++){
						next if exists $hits[$i]{'ignore'};
						next unless exists $seqs{$hits[$i]{'accn'}};
						$hits[$i]{'distance'} = shift (@a);
					}
				}
				unlink "$dir/temp", "$dir/temp.ph", "$dir/temp.dst", "$dir/buffer";
			}
		}

		## realign sequences, build phylogenetic tree and compute distance matrix with MAFFT (Neighbor-Joining) ##

		if ($realign eq "mafft" or $buildTree eq "mafft"){
			if (write_fasta("$dir/temp")>=3){
				system "$mafft --retree 1 --treeout --distout --quiet $dir/temp > $dir/buffer";
				die "Execution of MAFFT failed. Please check.\n" unless -s "$dir/buffer";
				read_fasta "$dir/buffer" if ($realign eq "mafft");
				if ($buildTree eq "mafft"){
					open IN, "<$dir/temp.tree";
					$tree = <IN>;
					close IN;
					$tree =~ s/\s+$/;/;
					$tree =~ s/_+:/:/g;
					$tree =~ s/([,\(])\d{1,4}_/$1/g;
					$tree =~ s/([,\(])$_:/$1'$names{$_}':/ for (keys %names);
					if ($distance){
						my %mafftIDs = ();
						open IN, "<$dir/temp.hat2";
						$_ = <IN>; $_ = <IN>; $_ = <IN>;
						@a = ();
						while (<IN>){
							s/\s+$//; s/^\s+//; next if (/\. =/);
							@a = (@a, split (/\s/));
							last if @a >= $nSeq;
						}
						close IN;
						$j = 0;
						$hits[0]{'distance'} = "0.000";
						for ($i=1; $i<=$#hits; $i++){
							next if exists $hits[$i]{'ignore'};
							next unless exists $seqs{$hits[$i]{'accn'}};
							$hits[$i]{'distance'} = $a[$j++];
						}
					}
				}
				unlink "$dir/temp", "$dir/buffer", "$dir/temp.tree", "$dir/temp.hat2";
			}
		}

		## build phylogenetic tree with PhyML (Maximum Likelihood) ##

		if ($buildTree eq "phyml"){
			if (write_phylip("$dir/temp")>=3){
				$s = "$phyml -i $dir/temp -d aa --quiet";
				$s .= "-b $bsTree" if $bsTree;
				system $s;
				open IN, "<$dir/temp_phyml_tree.txt"; $tree = <IN>; close IN;
				$tree =~ s/\s+$//;
				$tree =~ s/([,\(])t(\d+):/$1$hits[$2]{'accn'}:/g;
				$tree =~ s/([,\(])$_:/$1'$names{$_}':/ for (keys %names);
				unlink "$dir/temp", "$dir/temp_phyml_stats.txt", "$dir/temp_phyml_tree.txt";
				unlink "$dir/temp_phyml_boot_stats.txt", "$dir/temp_phyml_boot_trees.txt" if $bsTree;
			}
		}

		## build phylogenetic tree with FastTree (Maximum Likelihood) and perform SH test ##

		if ($buildTree eq "fasttree"){
			if (write_fasta("$dir/temp")>=3){
				system "$fasttree < $dir/temp > fasttree.tmp -quiet";
				open IN, "<fasttree.tmp"; $tree = <IN>; close IN;
				$tree =~ s/\s+$//;
				$tree =~ s/([,\(])$_:/$1'$names{$_}':/ for (keys %names);
				unlink "$dir/temp", "fasttree.tmp";
			}
		}

		## build phylogenetic tree and compute distance matrix with RAxML (Maximum Likelihood) ##

		if ($buildTree eq "raxml"){
			open OUT, ">$dir/temp";
			foreach (keys %seqs){
				$i = length ($seqs{$_});
				last;
			}
			print OUT " $nSeq $i\n";
			for ($i=0; $i<=$#hits; $i++){
				next if exists $hits[$i]{'ignore'};
				next unless exists $seqs{$hits[$i]{'accn'}};
				print OUT "t$i".(" " x (13-length("t$i"))).$seqs{$hits[$i]{'accn'}}."\n";
			}
			close OUT;
			$s = "$raxml -m PROTGAMMA$aaModel -p 12345 -s $dir/temp -n tmp";
			$s .= " -f a -x 12345 -# $bsTree" if $bsTree;
			$s .= " > RAxML_screen.tmp";
			system $s;
			die "Execution of RAxML failed. Please check.\n" unless ((-s "RAxML_result.tmp") or (-s "RAxML_bipartitions.tmp"));
			if ($bsTree){ open IN, "RAxML_bipartitions.tmp"; }
			else { open IN, "RAxML_result.tmp"; }
			$tree = <IN>; close IN;
			$tree =~ s/\s+$//;
			$tree =~ s/([,\(])t(\d+):/$1$hits[$2]{'accn'}:/g;
			$tree =~ s/([,\(])$_:/$1'$names{$_}':/ for (keys %names);
			if ($distance){
				$s = "$raxml -f x -m PROTGAMMA$aaModel -s $dir/temp -n tmp2";
				if ($bsTree){ $s .= " -t RAxML_bipartitions.tmp"; }
				else { open IN, " -t RAxML_result.tmp"; }
				$s .= " > RAxML_screen.tmp";
				system $s;
				open IN, "<RAxML_distances.tmp2";
				@a = ();
				while (<IN>){
					s/\s+$//; next unless $_;
					next unless (/^t0/);
					@b = split (/\s+/);
					push @a, $b[2];
				}
				close IN;
				$hits[0]{'distance'} = "0.000000";
				for ($i=1; $i<=$#hits; $i++){
					next if exists $hits[$i]{'ignore'};
					next unless exists $seqs{$hits[$i]{'accn'}};
					$hits[$i]{'distance'} = shift (@a);
				}
				unlink "RAxML_info.tmp2", "RAxML_distances.tmp2";
			}
			unlink "$dir/temp", "$dir/temp.reduced";
			unlink "RAxML_bestTree.tmp", "RAxML_info.tmp", "RAxML_log.tmp", "RAxML_parsimonyTree.tmp", "RAxML_result.tmp", "RAxML_screen.tmp";
			unlink "RAxML_bootstrap.tmp", "RAxML_bipartitions.tmp", "RAxML_bipartitionsBranchLabels.tmp" if $bsTree;
		}

		## print realigned sequences ##
		
		if ($realign or $trimSeq or $ignoreParalogs or $ignoreSubspecies or $distance){
			foreach (keys %seqs){
				$nChar = length ($seqs{$_});
				last;
			}
			$s = "";
			$reading = 1;
			open IN, "<$dir/$file" or next;
			while (<IN>){
				if (/^BEGIN TREE/ and $buildTree){ $reading = 0; next; }
				if ($reading) { $s .= $_; }
				if (/^END;/ and not $reading) { $reading = 1; next; }
				if (/^BEGIN ORGANISM/ and ($ignoreParalogs or $ignoreSubspecies or $distance)){
					$reading = 0;
					for ($i=0; $i<=$#hits; $i++){
						$s .= $hits[$i]{'accns'}."\t".$hits[$i]{'organism'}."\t".$hits[$i]{'group'}."\t".$hits[$i]{'taxid'}."\t".$hits[$i]{'score'}."\t".$hits[$i]{'expect'};
						if (exists $hits[$i]{'ignore'}){
							$s .= "\tx\n";
						}else{
							if ($distance and exists $hits[$i]{'distance'}){
								$s .= "\t".$hits[$i]{'distance'}."\n";
							}else{
								$s .= "\t\n";
							}
						}
					}
					$s .= ";\nEND;\n";
				}
				if (/^BEGIN DATA/){
					$reading = 0;
					$s .= "\tDIMENSIONS NTAX=$nSeq NCHAR=$nChar;\n";
					$s .= "\tFORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\n";
					$s .= "\tMATRIX\n";
					for ($i=0; $i<=$#hits; $i++){
						next if exists $hits[$i]{'ignore'};
						next unless exists $seqs{$hits[$i]{'accn'}};
						$s .= $hits[$i]{'accn'}."\t".$seqs{$hits[$i]{'accn'}}."\n";
					}
					$s .= ";\nEND;\n";
				}
			}
			close IN;
			open OUT, ">$dir/$file" or next;
			print OUT $s;
			close OUT;
		}

		## print tree ##
		
		if ($buildTree){
			open OUT, ">>$dir/$file" or next;
			print OUT "BEGIN TREES;\n\tTREE 1 = $tree\nEND;\n";
			close OUT;
		}
		
		## show progress ##
		
		if ($iProtein/$nProtein >= $iProgress/60){ 
			print ".";
			$iProgress++;
		}
	}
	
	unlink "RAxML_screen.tmp" if -e "RAxML_screen.tmp";
	
	print " done.\n";
}
print "Execution of purifier completed.\n";
exit 0;


## public subroutines ##

sub write_fasta ($){
	my $count = 0;
	open OUT, ">$_[0]";
	for ($i=0; $i<=$#hits; $i++){
		next if exists $hits[$i]{'ignore'};
		next unless exists $seqs{$hits[$i]{'accn'}};
		print OUT ">".$hits[$i]{'accn'}."\n".$seqs{$hits[$i]{'accn'}}."\n";
		$count ++;
	}
	close OUT;
	return $count;
}

sub write_phylip ($){
	my $count = 0;
	open OUT, ">$_[0]";
	foreach (keys %seqs){
		$i = length ($seqs{$_});
		last;
	}
	print OUT " $nSeq $i\n";
	for ($i=0; $i<=$#hits; $i++){
		next if exists $hits[$i]{'ignore'};
		next unless exists $seqs{$hits[$i]{'accn'}};
		print OUT "t$i".(" " x (13-length("t$i"))).$seqs{$hits[$i]{'accn'}}."\n";
		$count ++;
	}
	close OUT;
	return $count;
}


sub read_fasta ($){
	open IN, "<$_[0]";
	my $name = "";
	while (<IN>){
		s/\s+$//g;
		next unless $_;
		if (s/^>//){
			$name = $_;
			$seqs{$name} = "";
		}else{
			$seqs{$name} .= $_;
		}
	}
	close IN;
}

