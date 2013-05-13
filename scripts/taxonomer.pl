#!/usr/bin/perl

use strict;
use warnings;
$| = 1;

print "

This script looks up taxonomy information for organisms mentioned in BLAST reports

Usage:
	taxonomer.pl <working directory>

Output:
	taxonomy/self.info
	taxonomy/taxa.db
	taxonomy/ranks.db

" and exit unless @ARGV;

print "Executing taxonomer...\n";

my $i; my $s; my @a; my @b; my %h;

use LWP::Simple;

# use File::Basename;
# use LWP::UserAgent;
# use HTTP::Request;
# my $ua = LWP::UserAgent->new;

sub remote_query(@);
sub remote_query2();

my $wkDir = $ARGV[0];					# working directory
my $srcTaxa = 0;						# source taxonomy information
my $nRetry = 10;						# number of retries of remote look-up
my $delay = 1;							# time delay (seconds) between remote queries
my @ranks = ('species', 'genus', 'family', 'order', 'class', 'phylum');

if (-e "$wkDir/config.txt"){
	open IN, "<$wkDir/config.txt";
	while (<IN>){
		s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_; 
		$srcTaxa = $1 if /^nGenus=([012])$/;
		$nRetry = $1 if /^nRetry=(\d+)$/;
		$delay = $1 if /^delay=(\d+)$/;
		@ranks = split (/,/, $1) if /^ranks=(.+)$/;
	}
}

my %taxa = ();							# master storage of taxonomy information for each organism mentioned
										# in the BLAST reports. Keys are tax IDs. Each element is a hash
										# containing these elements:
										# organism (name),
										# lineage, in a form of "2/1224/28211/356/..." (high to low)
										# species, genus, family, order, class, phylum
										# each rank's value is its tax ID.
										# other intermediate ranks are not represented.

my %ranks = ();							# stores names of taxonomy ranks
my %sets = ();							# stores 'self' information


## Collect 'self' information ##

print "Identifying input organisms...";
opendir (DIR, "$wkDir/blast");
@a = readdir(DIR);
close DIR;
foreach my $set (@a){
	next if $set =~ /^\./;
	next unless -d "$wkDir/blast/$set";
	opendir (DIR, "$wkDir/blast/$set");
	@b = grep(/\.bla$/, readdir(DIR));
	close DIR;
	next unless @b;
	$b[0] =~ /^(.+)\.bla$/;
	%h = ('accn',$1,'accn0',$1);
	$sets{$set} = {%h};
}
@a = (); push (@a, $sets{$_}{'accn'}) for (keys %sets);
$s = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=". join (",", @a);
$s = get $s;
die "Retrieving taxonomic information failed.\n" unless $s =~ /<Count>\d+<\/Count>/s;
@a = (); push (@a, $1) while $s =~ s/<Id>(\d+)<\/Id>//s;
$s = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=";
$s .= $_."," for (@a); $s =~ s/,$//;
$s = get $s;
@a = (); push (@a, $1) while $s =~ s/<DocSum>(.+?)<\/DocSum>//s;
foreach $s (@a){
	$s =~ /<Item Name=\"Caption\" Type=\"String\">(\S+?)<\/Item>/;
	foreach (keys %sets){
		if ($sets{$_}{'accn0'} eq $1){ $i = $_; last; }
	}
	$s =~ /<Item Name=\"TaxId\" Type=\"Integer\">(\d+?)<\/Item>/;
	$sets{$i}{'taxid'} = $1;
	$s =~ /<Item Name=\"Title\" Type=\"String\">.*\[([^\[\]]+?)\]<\/Item>/;
	$sets{$i}{'organism'} = $1;
	%h = ('organism', $1);
	$taxa{$sets{$i}{'taxid'}} = {%h};
}
%h = (); $h{$sets{$_}{'taxid'}} = 1 for (keys %sets);
print " done. ". scalar (keys %h) ." organisms identified.\n";


## Summarize organisms mentioned in BLAST reports ##

print "Reading organisms mentioned in BLAST reports...\n  ";
foreach my $set (keys %sets){
	print "$set ";
	opendir (DIR, "$wkDir/blast/$set");
	my @files = grep(/\.bla$/, readdir(DIR));
	close DIR;
	next unless @b;
	foreach my $file (@files){
		open IN, "<$wkDir/blast/$set/$file" or next;
		my $reading = 0;
		while (<IN>) {
			s/\s+$//;
			if (/^BEGIN ORGANISM/){ $reading = "organism"; next; }
			if (/^END;/ or /^;/){ $reading = 0; next; }
			if ($reading eq "organism"){
				@a = split (/\t/);
				next if exists $taxa{$a[3]};
				%h = ('organism', $a[1]);
				$taxa{$a[3]} = {%h};
			}
		}
	}
}
my $n = scalar keys %taxa;
print "\ndone. $n organisms read.\n";


## Look up taxonomy information ##

my $iTaxon = 0; my $iProgress = 0;
print "Retrieving organism information...\n";
print "0-------------25-------------50------------75------------100%\n";

if ($srcTaxa == 0){
	my $count = 0; @a = ();
	foreach my $id (keys %taxa){
		next if exists $taxa{$id}{'lineage'};
		$taxa{$id}{'lineage'} = "";
		$taxa{$id}{$_} = "" for (@ranks);
		push @a, $id;
		$count ++; $iTaxon ++;
		if ($count == 20){
			remote_query @a;
			if ($iTaxon/$n >= $iProgress/60){ 
				print "."; $iProgress ++;
			}
			@a = (); $count = 0; sleep $delay;
		}
	}
	remote_query @a;
	print "done.\n";
	# for (0..$nRetry){
	#	remote_query2;
	# }
}

## Write taxonomy information ##
print "Writing taxonomy information...";
mkdir "$wkDir/taxonomy" unless -d "$wkDir/taxonomy";
open OUT, ">$wkDir/taxonomy/taxa.db";
foreach my $id (keys %taxa){
	print OUT $id."\t".$taxa{$id}{'organism'}."\t".$taxa{$id}{'lineage'};
	print OUT "\t".$taxa{$id}{$_} for (@ranks);
	print OUT "\n";
}
close OUT;
open OUT, ">$wkDir/taxonomy/ranks.db";
foreach my $id (keys %ranks){
	print OUT $id."\t".$ranks{$id}."\n";
}
close OUT;
print " done.\n";
open OUT, ">$wkDir/taxonomy/self.info";
foreach (keys %sets){
	print OUT $_."\t".$sets{$_}{'accn'}."\t".$sets{$_}{'taxid'}."\t".$sets{$_}{'organism'}."\n";
}
foreach my $rank (@ranks){
	$i = 1; $s = 0;
	foreach (keys %sets){
		unless ($s){
			$s = $taxa{$sets{$_}{'taxid'}}{$rank};
		}else{
			if ($s != $taxa{$sets{$_}{'taxid'}}{$rank}){
				$i = 0; last;
			}
		}
	}
	if ($i){
		print "The protein sets all belong to $rank $ranks{$s} ($s).\n";
		last;
	}
}
close OUT;

exit 0;


## Remote query (from NCBI server)

sub remote_query(@){
	$s = "http://www.ncbi.nlm.nih.gov/taxonomy/?report=info&term=".join (",", @_);
	$s = get $s;
	unless (defined $s) { sleep 10; remote_query @_; return; }
	die "Error: some taxa are not identified.\n" if ($s =~ /terms were not found/);
	my @results;
	push (@results, $1) while ($s =~ s/<div class=\"rslt\">(.+?)<\/div>//s);
	unless (@results){ sleep 10; remote_query @_; return; }
	foreach my $result (@results){
		$result =~ /<dt>Taxonomy ID:<\/dt> <dd><[^<>]*>(\d+)<\/span><\/dd>/;
		my $id = $1;
#		$_ =~ s/<.+?>//g;
		if ($result =~ /<dt>Rank:<\/dt> <dd>([^<>]+?)<\/dd>/){
			foreach (@ranks){
				if ($1 eq $_){ $taxa{$id}{$_} = $id; last; }
			}
		}
		while ($result =~ s/title="(.+?)".*?link_uid=(\d+).*?>(.+?)<\/a>//){
			$taxa{$id}{'lineage'} .= "/$2";
			foreach (@ranks){
				if ($1 eq $_){ $taxa{$id}{$_} = $2; last; }
			}
			$ranks{$2} = $3 unless exists $ranks{$2};
		}
	}
}

sub remote_query2(){ # old approach. slower.
	foreach my $id (keys %taxa){
		next if exists $taxa{$id}{'lineage'};
		select (undef, undef, undef, $delay); # delay certain amount of time
		print $id."\t".$taxa{$id}{'organism'}."\n";
		# my $request = new HTTP::Request GET => "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$id";
		# my $response = $ua->request($request);
		# my @taxReport = split (/\n/, $response->content);
		my @taxReport = split (/\n/, "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$id");
		my $type = 0;
		$taxa{$id}{'lineage'} = "";
		$taxa{$id}{$_} = "" for (@ranks);
		for ($i = 0; $i <= $#taxReport; $i ++){
			if ($taxReport[$i] =~ /<em>Rank: <\/em>(.+?)<br>/){
				next if $1 eq "no rank";
				foreach (@ranks){
					if ($1 eq $_){ $taxa{$id}{$_} = $id; last; }
				}
			}
			if ($taxReport[$i] =~ /^<dd>/){ $type = 1; last; }
			if ($taxReport[$i] =~ /^<A HREF=.+Lineage.+\(full\)/){ $type = 2; last; }
		}
		$s = $taxReport[$i];
		$s =~ s/<[\/]?dd>// if ($type == 1);
		$s =~ s/^.+full// if ($type == 2);
		print "  Identification of taxid $id failed. Queued to retry.\n" and next unless $type;
		@a = split (/; /, $s);
		foreach (@a){
			/id=(\d+).+TITLE="(.+)">(.+)<\/a>/ if ($type == 1);
			/id=(\d+).+TITLE="(.+)" .+>(.+)<\/A>/ if ($type == 2);
			$taxa{$id}{'lineage'} .= "/$1";
			foreach (@ranks){
				if ($2 eq $_){ $taxa{$id}{$_} = $1; last; }
			}
			$ranks{$1} = $3 unless exists $ranks{$1};
		}
	}
}
