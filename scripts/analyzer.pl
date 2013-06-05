#!/usr/bin/perl

use warnings;
use strict;
$| = 1;

print "

This program analyzes BLAST results of proteins from each genome
and identify ones that may have undergone horizontal gene transfer
events and their donor and/or recipients, as well as any potential
gene origination and loss events.

The program considers every folder in the working directory as a
genome containing BLAST results, and analyze its contents. The
result is saved as a tab-delimited file, which can be visualized
by programs such as Excel.

Usage:
	analyzer.pl <working directory>

Output:
	result/<name>.txt

" and exit unless @ARGV;

print "Executing analyzer...\n";


## all-purpose variables ##

my $i; my $j; my $n; my $s; my $t; my @a; my @b; my @c; my %h;


## subroutines ##

sub median(@);
sub mad(@);
sub quantiles(@);
sub z_scores(@);
sub z_test(@);
sub modified_z(@);
sub boxplot(@);
sub adjusted_boxplot(@);
sub recurse_deOutlier(@);
sub recurse_Z(@);


## global variables ##

my @sets;											# proteins sets to analyze
my @lvs;												# taxonomic ranks to analyze
my %selves = ();									# taxonomic information of self group
my %lvList = ();									# list of taxids of each level

my %taxadb = ();									# taxa.db
my %ranksdb = ();									# ranks.db
my %selfinfo = ();								# self.info

my %proteins = ();								# accn -> name (identified by COG)

my @files;											# blast report files for each protein set

## the master table storing everything. It's an array of hashes. Each row is a record.
my %results = ();

## the phyletic pattern of the whole genome. I name it as "fingerprint"
my %fpN = ();										# number of hits per blast
# my %fpS = ();									# score of each hit of each blast
														# not available in this version

## program parameters ##

my $wkDir = $ARGV[0];							# working directory
my $interactive = 1;								# interactive or automatic mode

my $minHits = 0;									# minimal number of hits a valid blast report should contain
my $maxHits = 0;									# maximal number of hits to retain from one blast, 0 means infinite
my $minSize = 0;									# minimal size (aa) of a valid protein (0 means infinite)
my $evalue = 1e-5;								# E-value cutoff
my $deHypo = 0;									# ignore hypothetical proteins
my $smOrthology;									# file containing scheme of orthology

# algorithm
my $selfRank = 0;									# taxonomic rank(s) on which the program analyzes
my $normalize = 1;								# use relative bit score (bit score of subject / bit score of query)
my $unite = 1;										# Blast pattern (0: each genome has own pattern, 1: one pattern for all genomes)

my $useDistance = 0;								# use phylogenetic distance instead of BLAST bit scores
my $useWeight = 1;								# use weight (sum of scores) instead of number of hits

# fingerprints
my $outRaw = 1;									# output raw number/weight data
my $outFp = 1;										# output fingerprint
my $graphFp = 0;									# graph fingerprint (requires R)

my $boxPlot = 1;									# box plot
my $histogram = 1;								# histogram
my $densityPlot = 1;								# density plot
my $histDensity = 0;								# overlap histogram and density plot
my $merge3Groups = 0;							# merge 3 groups in density plot
my $scatterPlot = 1;								# scatter plot
my $plot3D = 0;									# 3-way scatter plot

# cutoffs
my $howCO = 4;										# how to determine cutoffs (0: user-defined global cutoff (%), 1: user-defined individual
														  # cutoffs, 2: wait for user input, 3: histogram, 4: kernel density estimation,
														  # 5: hierarchical clustering)
    
my $globalCO = 0.25;								# arbitrary global cutoff (%)
my ($selfCO, $closeCO, $distalCO) = (0, 0, 0);		# user-defined cutoffs for individual groups

my $exOutlier = 0;								# exclude outliers, hits distant from other hits of the same group

my $nBin = 20;										# number of bins in histogram
my $plotHist = 0;									# plot histogram on screen

my $toolKDE = 0;									# computational tool for kernel density estimation
my $bwF = 1;										# bandwidth selection factor
my $plotKDE = 0;									# plot density function on screen
my $toolExtrema = 0;								# computational tool for identifying local extrema of density function (0: Perl code, 1: R package "pastecs")
my $modKCO = 1;									# location of cutoff (0: 1st pit, 1: midpoint of x-coordinates between 1st peak and 1st pit, 2/3: quantile
my $qKCO = 0.5;									# horizontal/vertical quantile from 1st pit toward 1st peak

my $dipTest = 0;									# perform non-unimodality test (Hartigan's dip test) and report p-value
my $dipSig = 0;									# use global cutoff if dip test's result is not significant

my $selfLow = 0;									# HGT-derived genes must have low self weight (an optional criterion)

my $BBH = 0;										# use conventional best match method instead
my $loss = 0;										# also report gene loss events
my $POE = 0;										# also report POE

my @ranks = ('species', 'genus', 'family', 'order', 'class', 'phylum');

my @selfGroup = ();
my @closeGroup = ();
my @inSets = ();
my @exSets = ();

my $R;												# Statistics::R instance

## read configurations ##

if (-e "$wkDir/config.txt"){
	open IN, "<$wkDir/config.txt";
	my $readMonitor = 0; my $readMiddle = 0;
	while (<IN>){
		s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_;
		$interactive = $1 if /^interactive=([01])$/;
		$minHits = $1 if /^minHits=(\d+)$/;
		$maxHits = $1 if /^maxHits=(\d+)$/;
		$minSize = $1 if /^minSize=(\d+)$/;
		$evalue = $1 if /^evalue=(.+)$/;
		$deHypo = $1 if /^deHypo=([012])$/;
		$smOrthology = $1 if /^smOrthology=(.+)$/;

		$selfRank = $1 if /^selfRank=(.+)$/;
		$normalize = $1 if /^normalize=(.+)$/;
		$unite = $1 if /^unite=([01])$/;
		$useDistance = $1 if /^useDistance=([01])$/;
		$useWeight = $1 if /^useWeight=([01])$/;

		$howCO = $1 if /^howCO=(\d)$/;
		$globalCO = $1 if /^globalCO=(.+)$/;
		$selfCO = $1 if /^selfCO=(.+)$/;
		$closeCO = $1 if /^closeCO=(.+)$/;
		$distalCO = $1 if /^distalCO=(.+)$/;

		$exOutlier = $1 if /^exOutlier=([0123])$/;

		$nBin = $1 if /^nBin=(\d+)$/;
		$plotHist = $1 if /^plotHist=([01])$/;

		$toolKDE = $1 if /^toolKDE=([01])$/;
		$bwF = $1 if /^bwF=(.+)$/;
		$plotKDE = $1 if /^plotKDE=([01])$/;
		$toolExtrema = $1 if /^toolExtrema=([01])$/;
		$modKCO = $1 if /^modKCO=([0123])$/;
		$qKCO = $1 if /^qKCO=(.+)$/;

		$dipTest = $1 if /^dipTest=([01])$/;
		$dipSig = $1 if /^dipSig=(.+)$/;

		$selfLow = $1 if /^selfLow=([01])$/;
		
		$outRaw = $1 if /^outRaw=([01])$/;
		$outFp = $1 if /^outFp=([01])$/;
		$graphFp = $1 if /^graphFp=([01])$/;
		$boxPlot = $1 if /^boxPlot=([01])$/;
		$histogram = $1 if /^histogram=([01])$/;
		$densityPlot = $1 if /^densityPlot=([01])$/;
		$histDensity = $1 if /^histDensity=([01])$/;
		$merge3Groups = $1 if /^merge3Groups=([01])$/;
		$scatterPlot = $1 if /^scatterPlot=([01])$/;
		$plot3D = $1 if /^plot3D=([01])$/;
	
		$BBH = $1 if /^BBH=([01])$/;
		$loss = $1 if /^loss=([01])$/;
		$POE = $1 if /^POE=([01])$/;
		
		@ranks = split (/,/, $1) if /^ranks=(.+)$/;

		@selfGroup = split (/,/, $1) if /^selfGroup=(.+)$/;
		@closeGroup = split (/,/, $1) if /^closeGroup=(.+)$/;
		@inSets = split (/,/, $1) if /^inSets=(.+)$/;
		@exSets = split (/,/, $1) if /^exSets=(.+)$/;
	}
	close IN;
}

## verify configurations ##

if ($globalCO){
	$globalCO = $globalCO / 100 if ($globalCO =~ s/%$//);
	die "Error: Global cutoff must be between 0 and 1.\n" if ($globalCO <= 0 or $globalCO >=1);
}

# conditionally use Statistics::R for Perl-R communication;
if ($graphFp or ($howCO == 4 and ($toolKDE or $toolExtrema)) or ($howCO == 5) or $dipTest){
	eval{ require Statistics::R; Statistics::R->import() };
	die "Error: Perl module Statistics::R is not available.\n" if ($@);
	$R = Statistics::R->new();
}

## initiate global variables ##

if ($unite){
	@b = (); @a = ('0','1','2'); #,'01','12','012'); 
	$fpN{'0'}{$_}{'data'} = [@b] for (@a);
	# $fpS{'0'}{$_}{'data'} = [@b] for (@a);
}


## Identify taxonomic levels ##

print "Reading taxonomic information...";
open IN, "<$wkDir/taxonomy/taxa.db";
while (<IN>){
	s/\s+$//; next if /^#/; next unless $_;
	@a = split /\t/;
	%h = ('name',$a[1],'rank',$a[2]);
	$i = 3; $h{$_} = $a[$i++] for (@ranks);
	$taxadb{$a[0]} = {%h};
}
close IN;
open IN, "<$wkDir/taxonomy/ranks.db";
while (<IN>){
	s/\s+$//; next if /^#/; next unless $_;
	@a = split /\t/;
	$ranksdb{$a[0]} = $a[1];
}
close IN;
open IN, "<$wkDir/taxonomy/self.info";
while (<IN>){
	s/\s+$//; next if /^#/; next unless $_;
	@a = split /\t/;
	next if ($#a < 3);
	%h = ('accn',$a[1],'taxid',$a[2],'name',$a[3]);
	$selfinfo{$a[0]} = {%h};
}
close IN;
print " done.\n";

print "Analyzing taxonomic information...";
for ($i=0; $i<=$#ranks; $i++){
	$j = 1; # if all names are the same
	$s = 0; # last name of this rank
	foreach my $key (keys %selfinfo){
		$t = $taxadb{$selfinfo{$key}{'taxid'}}{$ranks[$i]};
		$s = $t unless $s;
		if ($t != $s){ $j = 0; last; }
	}
	last if $j;
}
print " done.\n";
$selves{'rank'} = $ranks[$i];
$selves{'rank_id'} = $i;
$selves{'rank_taxid'} = $s;
$selves{'rank_name'} = $ranksdb{$s};
print "  All input genomes belong to $selves{'rank'} $selves{'rank_name'}.\n";

foreach (split (/,/, $selfRank)){
	if (/\d/){
		%h = ('rank',$ranks[$selves{'rank_id'}+$_],'id',$selves{'rank_id'}+$_);
	}else{
		%h = ('rank',$_);
		for ($j=0; $j<=$#ranks; $j++){
			if ($ranks[$j] eq $_){ $h{'id'} = $j; last; }
		}
	}
	foreach my $key (keys %selfinfo){ # information of this level
		$h{'taxid'} = $taxadb{$selfinfo{$key}{'taxid'}}{$h{'rank'}};
		$h{'name'} = $ranksdb{$h{'taxid'}};
		$i = $h{'id'}+1;
		if ($i <= $#ranks){
			$h{'parentRank'} = $ranks[$i];
			$h{'parentTaxid'} = $taxadb{$selfinfo{$key}{'taxid'}}{$ranks[$i]};
			$h{'parentName'} = $ranksdb{$h{'parentTaxid'}};
		}
		last;
	}
	
	if (@selfGroup){ # user-defined self group
		$h{'rank'} = "(user defined self group)"; $h{'id'} = -1; $h{'taxid'} = join (",", @selfGroup);
		@a = ();
		foreach (@selfGroup){
			if (exists $ranksdb{$_}){ push @a, $ranksdb{$_}; }
			elsif (exists $taxadb{$_}) { push @a, $taxadb{$_}{'name'}; }
			else { push @a, "unknown"; }
		}
		$h{'name'} = join (",", @a);
	}
	if (@closeGroup){ # user-defined close group
		$h{'parentRank'} = "(user defined close group)"; $h{'parentTaxid'} = join (",", @closeGroup);
		@a = ();
		foreach (@closeGroup){
			if (exists $ranksdb{$_}){ push @a, $ranksdb{$_}; }
			elsif (exists $taxadb{$_}) { push @a, $taxadb{$_}{'name'}; }
			else { push @a, "unknown"; }
		}
		$h{'parentName'} = join (",", @a);
	}	
	push @lvs, {%h}; # merge data
	
	my %ids = ();
	unless (exists $lvList{$h{'rank'}}){ # list of taxids of same level
		foreach my $key (keys %taxadb){
			next unless $taxadb{$key}{'rank'};
			foreach (split (/,/, $h{'taxid'})){
				$ids{$key} = 1 if ($taxadb{$key}{'rank'} =~ /\/$_$/ or $taxadb{$key}{'rank'} =~ /\/$_\//);
			}
			# next unless $taxadb{$key}{$h{'rank'}};
			# $ids{$key} = 1 if ($taxadb{$key}{$h{'rank'}} eq $h{'taxid'});
		}
		$lvList{$h{'rank'}} = {%ids};
	}
	%ids = ();
	unless (exists $lvList{$h{'parentRank'}}){ # list of taxids of parent level
		foreach my $key (keys %taxadb){
			next unless $taxadb{$key}{'rank'};
			foreach (split (/,/, $h{'parentTaxid'})){
				$ids{$key} = 1 if ($taxadb{$key}{'rank'} =~ /\/$_$/ or $taxadb{$key}{'rank'} =~ /\/$_\//);
			}
			# next unless $taxadb{$key}{$h{'parentRank'}};
			# $ids{$key} = 1 if ($taxadb{$key}{$h{'parentRank'}} eq $h{'parentTaxid'});
		}
		$lvList{$h{'parentRank'}} = {%ids};
	}
}
print "Analysis will work on the following taxonomic ranks:\n";
for ($i=0; $i<=$#lvs; $i++){
	print "  Self: $lvs[$i]{'rank'} $lvs[$i]{'name'} (taxid: $lvs[$i]{'taxid'}) (".(keys %{$lvList{$lvs[$i]{'rank'}}})." members),\n";
	print "  Close: $lvs[$i]{'parentRank'} $lvs[$i]{'parentName'} (taxid: $lvs[$i]{'parentTaxid'}) (".(keys %{$lvList{$lvs[$i]{'parentRank'}}})." members),\n";
	print "  Distal: all other organisms.\n";
}
if ($interactive){
	print "Press Enter to continue, or press Ctrl+C to exit:";
	$s = <STDIN>;
}



# Test whether levels are successfully detected.
# open OUT, ">self.txt"; foreach (sort keys %{$lvList{$lvs[0]{'rank'}}}){ print OUT $_."\t".$taxadb{$_}{'name'}."\n"; } close OUT;
# open OUT, ">parent.txt"; foreach (sort keys %{$lvList{$lvs[0]{'parentRank'}}}){ print OUT $_."\t".$taxadb{$_}{'name'}."\n";} close OUT;
# exit 1;


## Identify gene orthology ##

if ($deHypo == 2){
	if ($smOrthology){
		unless (-e $smOrthology){
			print "Orthology scheme file $smOrthology not accessible.\n";
			$smOrthology = "$wkDir/taxonomy/orthology.db";
		}
	}else{ $smOrthology = "$wkDir/taxonomy/orthology.db"; }
	unless (-e $smOrthology){ system "perl scripts/orthologer.pl $wkDir"; }
	unless (-e $smOrthology){ print "Error: Identification of orthology failed.\n" and exit 1; }
	open IN, "<$smOrthology";
	while (<IN>){
		s/\s+$//; next if /^#/; next unless $_;
		@a = split (/\t/); next unless $#a;
		$a[0] =~ s/^\d+\|//;
		@b = split (/[\s\/,]/, $a[$#a]);
		foreach (@b){
			$proteins{$_} = $a[0] unless exists $proteins{$_};
		}
	}
	close IN;
	$deHypo = 1 unless %proteins # if no name is defined by the scheme file, then use original protein name.
}


## Read protein sets ##

print "Reading protein sets...";
opendir (DIR, "$wkDir/blast");
@a = readdir(DIR);
close DIR;
foreach (@a){
	next if (/^\./);
	push @sets, $_ if -d "$wkDir/blast/$_";
}
print "done. ";
die "No genome detected.\n" unless @sets;
print @sets." sets detected.\n";


# Summarize blast reports ##

print "Analyzing Blast results...\n";
print "0-------------25-------------50------------75------------100%\n";

foreach my $set (@sets){
	if (@inSets){ $i = 0; foreach (@inSets){ if ($set eq $_){ $i = 1; last; } } next unless $i; }
	if (@exSets){ $i = 0; foreach (@exSets){ if ($set eq $_){ $i = 1; last; } } next if $i; }
	opendir (DIR, "$wkDir/blast/$set");
	@files = grep(/\.bla$/,readdir(DIR));
	close DIR;
	print "No protein found in $set\n" and next unless @files;

	## varibles to show a progress bar
	my $iProtein = 0;
	my $iProgress = 0;
	my $nProtein = $#files+1;
	print "$set has $nProtein proteins. Analyzing...\n";

	## varibles to create phyletic patterns
	my %selfNumbers = (); my %closeNumbers = (); my %distalNumbers = ();
	my %selfScores = (); my %closeScores = (); my %distalScores = ();
	
	unless ($unite){
		@a = ('0','1','2'); @b = ();
		# @a = ('0','1','2','01','12','012'); @b = ();
		$fpN{$set}{$_}{'data'} = [@b] for (@a);
		# $fpS{$set}{$_}{'data'} = [@b] for (@a);
	}	
	
	## Information of self
	my %self = ();

	foreach my $file (@files){
		$iProtein ++;
		
		my %result = ();				# a record to store everything about this blast search, including:
										# query, length, product, n,
										# details for sGenus, sMiddle, sGroup, nGenus, nMiddle, nGroup, sMiddle_nGenus, sGroup_nMiddle, all
										## each of which contains:
										## n, mean, median, min, max, q1, q3, stdev, mad.
										# names of top hit in nGenus, nMiddle, nGroup
										
		my %scores = ();				# scores of each hit by category, as a buffer for computing the statistics above
		my @hits;						# parameters of the hits. one hit contains:
										#accn, organism, group, taxid, genus, score
		$file =~ /(.+)\.[^.]+$/;
		$result{'query'} = $1;

		if ($deHypo == 2 and exists $proteins{$result{'query'}}){
			$s = $proteins{$result{'query'}};
			next if ($s =~ /hypothetical/ or $s =~ /hypotethical/ or $s =~ /hypothetcial/);
		}

		my $nHits = 0;					# total number of hits. just for convenience
		my $nScore = 0;					# total blast score
		my $lastHit = "";				# store the last hit in the organism table, for the identification of duplicated taxa
		my $selfScore = 0;


		# Read hit table #
		open IN, "<$wkDir/blast/$set/$file" or next;
		my $reading = 0;
		while (<IN>) {
			s/\s+$//;
			if (/^BEGIN QUERY/){ $reading = "query"; next; }
			if (/^BEGIN ORGANISM/){ $reading = "organism"; next; }
			if (/^BEGIN DATA/){ $reading = "data"; next; }
			if (/^END;/){ $reading = 0; next; }
			if ($reading eq "query"){ # read query (self)
				if (/^\tAccession=(.+);$/){
					$result{'accn'} = $1;
					$result{'accn'} =~ s/\.[\d]+$//;
				}
				$result{'gi'} = $1 if /^\tLength=(\d+);$/;
				$result{'length'} = $1 if /^\tLength=(.+);$/;
				$result{'product'} = $1 if /^\tProduct=(.+)\s*;$/;
				$result{'organism'} = $1 if /^\tOrganism=(.+)\s*;$/;
			}
			if ($reading eq "organism"){ # read organisms
				next if /^;/;
				@a = split (/\t/);
				next if ($#a >= 6 and $a[6] eq "x");
				next if ($evalue and $a[5] > $evalue);
				my %hit = ();
				$hit{'accns'} = $a[0];
				$hit{'organism'} = $a[1];
				$hit{'group'} = $a[2];
				$hit{'taxid'} = $a[3];
				$hit{'evalue'} = $a[5];
				if ($useDistance){
					if ($#a >= 6){ $hit{'score'} = 1 - $a[6]; } # use distance instead of bit score
					elsif (!@hits){ $hit{'score'} = 0; }
					else { $hit{'score'} = $hits[$#hits]{'score'}; }
				}else{
					$hit{'score'} = $a[4];
				}
				@a = split(/\//, $a[0]);
				$hit{'accn'} = $a[0];
				push @hits, {%hit};
			}
			last if ($maxHits and $#hits >= $maxHits-1);
		}
		close IN;
		next unless @hits;
		# next if ($minSize and ($result{'length'} < $minSize));
		print "Problems reading $file of $set.\n" unless (exists $result{'query'} and exists $result{'length'} and exists $result{'product'});

		if ($deHypo == 1 and $result{'product'}){
			$s = $result{'product'};
			next if ($s =~ /hypothetical/ or $s =~ /hypotethical/ or $s =~ /hypothetcial/);
		}


		## Intepret hit table ##

		# total number of hits #
		$result{'n'} = $#hits+1;

		# sort by score or distance #
		@hits = sort {$b->{'score'} <=> $a->{'score'}} @hits;

		# identify self (query) information #
		for ($i=0; $i<=$#hits; $i++){
			@a = split(/\//, $hits[$i]{'accns'});
			foreach (@a){
				if ($result{'accn'} eq $_){
					$result{'id'} = $i;
					$result{'taxid'} = $hits[$i]{'taxid'};
					$result{'score'} = $hits[$i]{'score'};
					$result{'organism'} = $hits[$i]{'organism'};
					last;
				}
			}
			last if exists $result{'id'};
		}
		unless (exists $result{'id'}){
			$result{'id'} = 0;
			$result{'taxid'} = $hits[0]{'taxid'};
			$result{'score'} = $hits[0]{'score'};
			$result{'organism'} = $hits[0]{'organism'};		
		}

		# Use absolute or relative bit scores #
		if ($normalize and not $useDistance){
			for ($i=0; $i<=$#hits; $i++){
				$hits[$i]{'score'} = sprintf("%.3f", $hits[$i]{'score'}/$result{'score'});
			}
		}

		# initialize values of prediction results
		
		$result{'in'} = "";					# whether incoming HGT or origination took place within the group
		$result{'loss'} = "";				# gene loss event
		$result{'origin'} = "";				# gene origination event
		$result{'income'} = "";				# incoming HGT event
		$result{'outcome'} = "";			# outcoming HGT event

		# Summarize numbers and scores #
		## 0 - self group, 1 - close groups, 2 - distal groups
		## N - number, S - scores
		## hit1 - first close hit, hit 2 - first distal hit
		
		$result{'N0'} = 0; $result{'N1'} = 0; $result{'N2'} = 0;
		@a = (); $result{'S0'} = [@a]; $result{'S1'} = [@a]; $result{'S2'} = [@a];
		for ($i=0; $i<=$#hits; $i++){
			if (exists $lvList{$lvs[0]{'rank'}}{$hits[$i]{'taxid'}}){
				if ($useWeight){ $result{'N0'} += $hits[$i]{'score'}; }
				else { $result{'N0'} ++; }
				push @{$result{'S0'}}, $hits[$i]{'score'};
			}elsif (exists $lvList{$lvs[0]{'parentRank'}}{$hits[$i]{'taxid'}}){
				if ($useWeight){ $result{'N1'} += $hits[$i]{'score'}; }
				else { $result{'N1'} ++; }
				push @{$result{'S1'}}, $hits[$i]{'score'};
				$result{'hit1'} = $hits[$i]{'taxid'} unless exists $result{'hit1'};
				$result{'BBH'} = "" unless exists $result{'BBH'};
			}else{
				if ($useWeight){ $result{'N2'} += $hits[$i]{'score'}; }
				else { $result{'N2'} ++; }
				push @{$result{'S2'}}, $hits[$i]{'score'};
				unless (exists $result{'hit2'}){
					$result{'hit2'} = $hits[$i]{'taxid'};
					$result{'hit2accn'} = $hits[$i]{'accn'};
				}
				$result{'BBH'} = 1 unless exists $result{'BBH'};
			}
		}
		
		$result{'N0'} = sprintf("%.3f", $result{'N0'});
		$result{'N1'} = sprintf("%.3f", $result{'N1'});
		$result{'N2'} = sprintf("%.3f", $result{'N2'});
#		$result{'N01'} = $result{'N0'} + $result{'N1'};
#		$result{'N12'} = $result{'N1'} + $result{'N2'};
#		$result{'N012'} = $result{'N0'} + $result{'N1'} + $result{'N2'};

		# proteins without non-self hits are considered as de novo originated and not considered for statistics

		$result{'origin'} = 1 if $result{'N1'}+$result{'N2'} < 0.000001;

		# proteins with hits lower than minimum cutoff are considered as de novo originated and not considered for statistics

		if ($minHits and ($result{'n'} < $minHits)){ $result{'origin'} = 1; $result{'BBH'} = ""; }

		# Put this record into the master record
		push @{$results{$set}}, {%result};
		
		# Record fingerprint #
		if ($unite){ $s = '0'; }
		else { $s = $set; }

		if ($result{'N0'} and not $result{'origin'}){ # skip in case distance is used and no distance is available
			push @{$fpN{$s}{'0'}{'data'}}, $result{'N0'}; # push @{$fpS{$s}{'0'}{'data'}}, @{$result{'S0'}};
			push @{$fpN{$s}{'1'}{'data'}}, $result{'N1'}; # push @{$fpS{$s}{'1'}{'data'}}, @{$result{'S1'}};
			push @{$fpN{$s}{'2'}{'data'}}, $result{'N2'}; # push @{$fpS{$s}{'2'}{'data'}}, @{$result{'S2'}};
#			push @{$fpN{$s}{'01'}{'data'}}, $result{'N01'}; # push @{$fpS{$s}{'01'}{'data'}}, (@{$result{'S0'}}, @{$result{'S1'}});
#			push @{$fpN{$s}{'12'}{'data'}}, $result{'N12'}; # push @{$fpS{$s}{'12'}{'data'}}, (@{$result{'S1'}}, @{$result{'S2'}});
#			push @{$fpN{$s}{'012'}{'data'}}, $result{'N012'}; # push @{$fpS{$s}{'012'}{'data'}}, (@{$result{'S0'}}, @{$result{'S1'}}, @{$result{'S2'}});
		}
		
		# Show progress
		print "." and $iProgress++ if ($iProtein/$nProtein >= $iProgress/60);
	}
	print "\n";
}
print " done.\n";

## create folder to contain results ##

mkdir "$wkDir/result" unless -d "$wkDir/result";
mkdir "$wkDir/result/statistics" unless -d "$wkDir/result/statistics";

## output raw data for further statistical analysis ##

if ($outRaw){
	open OUT, ">$wkDir/result/statistics/rawdata.txt";
	print OUT "Query\tSet\tLength\tHits\tSelf\tClose\tDistal\n";
	foreach my $set (sort keys %results){
		$n = @{$results{$set}};
		for ($i=0; $i<$n; $i++){
			my %res = %{$results{$set}[$i]};
			next if $res{'origin'};
			next unless $res{'N0'};
			print OUT "$res{'accn'}\t$set\t$res{'length'}\t$res{'n'}\t$res{'N0'}\t$res{'N1'}\t$res{'N2'}\n";
		}
	}
	close OUT;
	print "Raw data are saved in result/statistics/rawdata.txt.\n";
	print "You may conduct further analyses on these data.\n";
	if ($interactive){
		print "Press Enter to continue, or press Ctrl+C to exit:";
		$s = <STDIN>;
	}
}

## graph fingerprints with R ##

if ($graphFp){
	print "\nGraphing fingerprints with R...";
	$R->startR;
	print "R cannot be started. Make sure it is properly installed in the system.\n" and exit 1 unless $R->is_started();
	if ($plot3D){
		$R->send("library('rgl')");
		print "\n  You chose to display interactive 3D scatter plots. They will be displayed sequentially. Use mouse to rotate plots. Press Enter in the terminal to move to the next plot.\n";
	}
	foreach my $set (sort keys %fpN){
		my $fpre = "$wkDir/result/statistics/".("$set." x ($set ne "0")); # prefix of filename
		my $tpost = " of $set" x ($set ne "0"); # postfix of title
		for (0..2){
			@b = @{$fpN{$set}{$_}{'data'}};
			$_ = sprintf("%.3f", $_) for (@b);
			$R->send("x$_<-c(".join (",", @b).")");
		}
		if ($boxPlot){
			$R->send("pdf(\"$fpre"."box.pdf\")");
			$R->send("boxplot(x0,x1,x2,names=c(\"self\",\"close\",\"distal\"), main='Box plot$tpost', xlab='Group', ylab='Weight')");
			$R->send("dev.off()");
		}
		if ($histogram){
			$R->send("pdf(\"$fpre"."hist.pdf\", width=21, height=8)");
			$R->send("par(mfrow=c(1,3), oma=c(0,0,3,0))");
			for (0..2){
				if ($_ == 0){ $s = "Self"; }elsif ($_ == 1){ $s = "Close"; }else{ $s = "Distal"; }
				$R->send("hist(x$_, breaks=$nBin, freq=F, col='lightgrey', xlab='Weight', ylab='Probability density', main='$s')");
				$R->send("lines(density(x$_".(",bw=bw.nrd0(x$_)*$bwF" x ($bwF and $bwF != 1))."), lwd=2)") if $histDensity;
			}
			$R->send("title(\"Histogram".(" and density function" x $histDensity)."$tpost\",outer=T)");
			$R->send("dev.off()");
		}
		if ($densityPlot){
			if ($merge3Groups){
				$R->send("pdf(\"$fpre"."density.pdf\")");
				$R->send("plot(density(x0".(",bw=bw.nrd0(x0)*$bwF" x ($bwF and $bwF != 1))."), xlab='Weight', ylab='Probability density', main='Density plot$tpost')"); # xlim=range(0:1),
				$R->send("lines(density(x1".(",bw=bw.nrd0(x1)*$bwF" x ($bwF and $bwF != 1))."), col=2)");
				$R->send("lines(density(x2".(",bw=bw.nrd0(x2)*$bwF" x ($bwF and $bwF != 1))."), col=3)");
				$R->send("legend(\"topright\",legend=c(\"self\",\"close\",\"distal\"),col=(1:3),lwd=2,lty=1)");
				$R->send("dev.off()");
			}elsif (not ($histogram and $histDensity)){
				$R->send("pdf(\"$fpre"."density.pdf\", width=21, height=8)");
				$R->send("par(mfrow=c(1,3), oma=c(0,0,3,0))");
				for (0..2){
					if ($_ == 0){ $s = "Self"; }elsif ($_ == 1){ $s = "Close"; }else{ $s = "Distal"; }
					$R->send("plot(density(x$_".(",bw=bw.nrd0(x$_)*$bwF" x ($bwF and $bwF != 1))."), lwd=2, xlab='Weight', ylab='Probability density', main='$s')");
				}
				$R->send("title(\"Density plot$tpost\",outer=T)");
				$R->send("dev.off()");
			}
		}
		if ($scatterPlot){
			$R->send("pdf(\"$fpre"."scatter.pdf\", width=21, height=8)");
			$R->send("par(mfrow=c(1,3), oma=c(0,0,3,0))");
			$R->send("plot(x0, x1, xlab='Self weight', ylab='Close weight')");
			$R->send("plot(x1, x2, xlab='Close weight', ylab='Distal weight')");
			$R->send("plot(x2, x0, xlab='Distal weight', ylab='Self weight')");
			$R->send("title(\"Scatter plot$tpost\",outer=T)");
			$R->send("dev.off()");
		}
		if ($plot3D){
			# if ($^O=~/Win/){ $R->send("windows()"); }elsif ($^O=~/Mac/){ $R->send("quartz()"); }else{ $R->send("x11()"); }
			$R->send("plot3d(x0, x1, x2, xlab='Self weight', ylab='Close weight', zlab='Distal weight', title=\"3D scatter plot$tpost\")");
			print "Displaying 3D plot$tpost";
			$s = <STDIN>;
		}
	}
	$R->stopR();
	print " done.\n";
	print "Graphs are saved in result/statistics/.\n";
	if ($interactive){
		print "Press Enter to continue, or press Ctrl+C to exit:";
		$s = <STDIN>;
	}
}


## Compute the statistics for the whole genome(s), i.e., phyletic pattern, or "fingerprint" ##

print "\nComputing statistics...";
if (($howCO == 4 and ($toolKDE or $toolExtrema)) or ($howCO == 5) or $dipTest){
	$R->startR;
	print "R cannot be started. Make sure it is properly installed in the system.\n" and exit 1 unless $R->is_started();
	$R->send("library(pastecs)") if ($howCO == 4 and $toolExtrema);
	$R->send("library(diptest)") if $dipTest;
}

foreach my $set (keys %fpN){
	if ($set eq "0"){ print "\n  All protein sets:\n"; }else{ print "\n  Protein set $set:\n"; }
	foreach my $key ('0','1','2'){ #(keys %{$fpN{$set}}){
		next unless @{$fpN{$set}{$key}{'data'}}; # and @{$fpS{$set}{$key}{'data'}});
		
		if ($key eq "0"){ print "  Self group:\n"; }elsif ($key eq "1"){ print "  Close group:\n"; }elsif ($key eq "2"){ print "  Distal group:\n"; }
		
		@a = sort {$a<=>$b} @{$fpN{$set}{$key}{'data'}}; # sort low to high
		# @{$fpN{$set}{$key}{'data'}} = sort {$b<=>$a} @{$fpN{$set}{$key}{'data'}};
		# @a = @{$fpN{$set}{$key}{'data'}};
		# next unless @a;

		my $global_cutoff;
		my $computed_cutoff;
		my $use_global = 0;

		# compute basic statistical parameters
		
		$fpN{$set}{$key}{'n'} = @a;
		$s = 0; $s += $_ for @a;
		$fpN{$set}{$key}{'mean'} = $s/$fpN{$set}{$key}{'n'};
		$s = 0; $s += ($fpN{$set}{$key}{'mean'} - $_)**2 for @a;
		$fpN{$set}{$key}{'stdev'} = sqrt($s/($fpN{$set}{$key}{'n'}-1));
		$fpN{$set}{$key}{'min'} = $a[0];
		$fpN{$set}{$key}{'max'} = $a[$#a];
		$fpN{$set}{$key}{'median'} = median(@a);
		$fpN{$set}{$key}{'mad'} = mad(@a);
		($fpN{$set}{$key}{'q1'}, $fpN{$set}{$key}{'q3'}) = quantiles(@a);
		
		# determine cutoff using global cutoff

		$i = $fpN{$set}{$key}{'n'}*$globalCO;
		if (int($i) == $i){
			$global_cutoff = ($a[$#a-$i+1]+$a[$#a-$i])/2;
		}else{
			if ($i-int($i) <= int($i)-$i+1){
				$global_cutoff = $a[$#a-int($i)+1];
			}else{
				$global_cutoff = $a[$#a-int($i)];
			}
		}
		print "  Global cutoff ($globalCO) = $global_cutoff.\n";

		# override individual cutoff with global cutoff

		$use_global = 1 if ((($key eq '0') and ($selfCO eq 'G')) or (($key eq '1') and ($closeCO eq 'G')) or (($key eq '2') and ($distalCO eq 'G')));
		
		# exclude outliers

		if ($exOutlier){
			@c = boxplot(@a) if ($exOutlier == 1);
			@c = adjusted_boxplot(@a) if ($exOutlier == 2);
			@c = modified_z(@a) if ($exOutlier == 3);
			if ($a[$#a] > $c[1]){
				for ($i=$#a; $i>=0; $i--){
					if ($a[$i] <= $c[1]){
						@a = @a[0..$i];
						last;
					}
				}
			}
		}

		# perform Hartigan's dip test to assess non-unimodality

		if ($dipTest){
			print "  Performing Hartigan's dip test...";
			$R->send("x<-c(".join (",", @a).")");
			$R->send("dip.test(x)");
			$s = $R->read;
			if ($s =~ /D = \S+, p-value [<=] (\S+)\n/){
				print "  done.\n  $&";
				if ($dipSig){
					if ($1 >= $dipSig){
						print "  The weight distribution is NOT significantly non-unimodal.\n";
						if ($howCO >= 3){
							if ($interactive){
								print "  Proceed with statistical analysis anyway (yes) or use global cutoff ($globalCO) instead (no)? ";
								while (<STDIN>){
									chomp; last unless $_;
									last if (/^y$/i or /^yes$/i);
									if (/^n$/i or /^no$/i){ $use_global = 1; last; }
								}
							}else{ $use_global = 1; }
						}
					}else{ print "  The weight distribution is significantly non-unimodal.\n"; }
				}else{ print "  The weight distribution is ". ("NOT" x ($1 >= 0.05)). " significantly non-unimodal.\n"; }
			}else{ print "  failed.\n"; }
		}

		# determine cutoff using histogram

		my @freqs = (0)x$nBin;
		my $interval = ($a[$#a]-$a[0])/$nBin;
		my $cid = 0;
		for ($j=1; $j<$nBin; $j++){
			my $high_bound = $interval*$j;
			for ($i=$cid; $i<=$#a; $i++){
				if ($a[$i] < $high_bound){
					$freqs[$j-1] ++;
				}else{
					$cid = $i;
					last;
				}
			}
		}
		$freqs[$nBin-1] = $#a-$cid+1;
		my $local_min = 0; # index of the lowest bar from left
		for ($i=1; $i<$nBin-1; $i++){
			if ($freqs[$i]<$freqs[$i-1] and $freqs[$i]<=$freqs[$i+1]){
				$computed_cutoff = $interval*$i;
				$local_min = $i;
				last;
			}
		}
		
		# draw histogram
		
		if ($plotHist){
			print "  Histogram:";
			@c = sort {$b<=>$a} @freqs;
			$s = 50/$c[0];
			my @widths = (0)x$nBin;
			my @labels = (0)x$nBin;
			for ($i=0; $i<$nBin; $i++){
				$widths[$i] = int($freqs[$i]*$s);
				$labels[$i] = sprintf("%.2f", $interval*$i)."-".sprintf("%.2f", $interval*($i+1));
			}
			@c = sort {length($b)<=>length($a)} @labels;
			$s = length($c[0]);
			print "\n";
			for ($i=0; $i<$nBin; $i++){
				print " "x($s-length($labels[$i]))."$labels[$i] ".("*"x$widths[$i])."$freqs[$i]".(" (local minimum)"x($local_min and $local_min==$i))."\n";
			}
		}

		# determine cutoff using kernel density estimation (KDE)

		if ($howCO == 4 and not $use_global){
		
			# perform kernel density estimation (KDE)
		
			print "  Performing kernel density estimation...";
			my (@dx, @dy); # x- and y-coornidates of density function
			
			# perform KDE using basic R command "density"
			
			if ($toolKDE){
				$R->send("x<-c(".join (",", @a).")");
				if ($bwF and $bwF != 1){
					$R->send("bwx<-bw.nrd0(x)*$bwF");
					$R->send("d<-density(x,bw=bwx)");
				}else{
					$R->send("d<-density(x)");
				}
				@dx = @{$R->get('d$x')};
				@dy = @{$R->get('d$y')};
			}
			
			# perform KDE using self-written Perl code
			
			else{
				my $n = scalar @a;
				my $mean; $mean += $_ for @a; $mean = $mean/$n;
				my $stdev; $stdev += ($mean-$_)**2 for @a; $stdev = sqrt($stdev/($n-1));
				my @Q = quantiles(@a); my $iqr = $Q[1]-$Q[0];
				my $bw;
				if ($stdev == 0 and $iqr == 0){ $bw = 1; }
				elsif ($stdev == 0){ $bw = $iqr/1.34; }
				elsif ($iqr == 0){ $bw = $stdev; }
				elsif ($stdev <= $iqr/1.34){ $bw = $stdev; }
				else{ $bw = $iqr/1.34; }
				$bw = 0.9*$bw*$n**(-1/5); # select bandwidth by Silverman's ¡°rule of thumb¡± (1986)
				$bw = $bw*$bwF if ($bwF and $bwF != 1);
				my ($min, $max) = ($a[0]-3*$bw, $a[$#a]+3*$bw); # cut = 3
				print "\n  N = $n, bandwidth = ".sprintf("%.3f", $bw).".\n";
				for (my $x=$min; $x<=$max; $x+=($max-$min)/511) { # 512 points
					my $e = 0; $e += exp(-(($x-$_)/$bw)**2/2)/sqrt(2*3.1415926536) for @a; # Gaussian kernel
					push @dx, $x;
					push @dy, 1/$n/$bw*$e;
				}
			}
			
			print " done.\n";
			
			# plot density function on screen

			if ($plotKDE){
				print "  Density function:\n";
				my $k_x = int(@dx/100);
				$k_x = int(@dx/$plotKDE) if ($plotKDE > 1);
				my $k_y = 0;
				foreach (@dy){
					$k_y = $_ if ($_ > $k_y);
				}
				$k_y = 64/$k_y;
				for ($i=0; $i<=$#dy; $i+=$k_x){
					print " "x(7-length(sprintf("%.2f", $dx[$i]))).sprintf("%.2f", $dx[$i]);
					print " "x (int($dy[$i]*$k_y)+1)."*\n";
				}
			}

			my $peak1st;
			my $failed;
			my ($peak_i, $peak_x, $peak_y);
			my ($pit_i, $pit_x, $pit_y);

			# identify local extrema using R package "pastecs"
			
			if ($toolExtrema){
				if ($toolKDE < 2){
					$R->send("d<-data.frame(x=c(".join (",", @dx)."),y=c(".join (",", @dy)."))");
				}
				$R->send("tp<-turnpoints(ts(d\$y))");
				if ($R->get('tp$firstispeak') eq "TRUE"){
					if ($R->get('tp$nturns') == 1){ $failed = 1; }
					else{ $peak1st = 1; }
				}else{ $peak1st = 0; }
				# $R->send("summary(tp)");
				# $s = $R->read;
				# if ($s =~ /nbr turning points: 1 (first point is a peak)/){
				#	$failed = 1;
				# }else{
				#	if ($s =~ /first point is a peak/){
				#		$peak1st = 1;
				#	}elsif ($s =~ /first point is a pit/){
				#		$peak1st = 0;
				#	}else{
				#		$failed = 1;
				#	}
				#}

				unless ($failed){
					my @tpx = @{$R->get("d\$x[tp\$tppos]")};
					my @tpy = @{$R->get("d\$y[tp\$tppos]")};
					if ($peak1st){
						($peak_x, $peak_y) = ($tpx[0], $tpy[0]);
						($pit_x, $pit_y) = ($tpx[1], $tpy[1]);
					}else{
						($peak_x, $peak_y) = ($dx[0], $dy[0]);
						($pit_x, $pit_y) = ($tpx[0], $tpy[0]);
					}
					# $R->send("d\$x[tp\$tppos]");
					# $s = $R->read;
					# $s =~ s/\[\d+?\]//g;
					# $s =~ s/^\s+//;
					# @b = split(/\s+/, $s);
					# if ($peak1st){ # first extremum is a peak
					#	if ($stringency){ $computed_cutoff = ($b[0]+$b[1])/2; } # conservative: midpoint between peak and pit
					#	else{ $computed_cutoff = $b[1]; } # liberal: pit
					# }else{ # first extremum is a pit
					#	if ($stringency){ $computed_cutoff = $b[0]/2; } # conservative: midpoint between pit and zero
					#	else{ $computed_cutoff = $b[0]; } # liberal: pit
					# }
				}
			}
			
			# identify local extrema using self-written Perl code
			
			else{
				for ($i=0; $i<$#dy; $i++){
					if ($dy[$i] < $dy[$i+1]){ $peak1st = 1; last; }
					elsif ($dy[$i] > $dy[$i+1]){ $peak1st = 0; last; }
					else{ next; }
				}
				if ($peak1st){
					for ($i=1; $i<$#dy; $i++){
						if (($dy[$i-1] <= $dy[$i]) and ($dy[$i] > $dy[$i+1])){
							($peak_i, $peak_x, $peak_y) = ($i, $dx[$i], $dy[$i]);
						}
						if (($dy[$i-1] > $dy[$i]) and ($dy[$i] <= $dy[$i+1])){
							($pit_i, $pit_x, $pit_y) = ($i, $dx[$i], $dy[$i]);
						}
						last if ($peak_i and $pit_i);
					}
					$failed = 1 unless $pit_i;
				}else{
					($peak_i, $peak_x, $peak_y) = (0, $dx[0], $dy[0]);
					for (my $i=1; $i<$#dx; $i++){
						if (($dy[$i-1] > $dy[$i]) and ($dy[$i] <= $dy[$i+1])){
							($pit_i, $pit_x, $pit_y) = ($i, $dx[$i], $dy[$i]);
						}
						last if $pit_i;
					}
				}
			}
			
			# locate cutoff point
			
			if ($modKCO == 0){ # pit
				$computed_cutoff = $pit_x;
			}elsif ($modKCO == 1){ # horizontal midpoint
				$computed_cutoff = ($pit_x+$peak_x)/2;
			}elsif ($modKCO == 2){ # horizontal quantile
				$computed_cutoff = $pit_x-($pit_x-$peak_x)*$qKCO;
			}else{ # vertical quantile
				my $vCO = 0;
				$vCO = $pit_y+($peak_y-$pit_y)*$qKCO;
				for (my $i=$peak_i; $i<$pit_i; $i++){
					if (($dy[$i] >= $vCO) and ($vCO > $dy[$i+1])){
						$computed_cutoff = $dx[$i]+($dx[$i+1]-$dx[$i])*($dy[$i]-$vCO)/($dy[$i]-$dy[$i+1]);
						last;
					}
				}
			}

		}
		
		# determine cutoff using hierarchical clustering

		if ($howCO == 5 and not $use_global){
			print "  Performing hierarchical clustering...";
			$R->send("x<-c(".join (",", @a).")");
			$R->send("d<-dist(x,method=\"euclidean\")");
			$R->send("fit<-hclust(d,method=\"ward\")");
			print "  done.\n";
			my $nCluster = 1;
			while ($nCluster++){
				$R->send("c<-cutree(fit,k=$nCluster)");
				my $clusters = $R->get('c'); # tip: return vector from R
				@c = (()) x $nCluster; # data of clusters
				for ($i=0; $i<=@{$clusters}-1; $i++){
					foreach (1..$nCluster){
						if (@{$clusters}[$i] == $_){
							push @{$c[$_-1]}, $a[$i];
							last;
						}
					}
				}
				for ($i=1; $i<=$#c; $i++){
					@{$c[$i]} = sort {$a <=> $b} @{$c[$i]};
				}
				@b = sort {$a->[0] <=> $b->[0]} @c;
				print "  With $nCluster clusters, cutoff is $b[1][0]. Accept? (yes/no) ";
				my $user_okay = 1;
				while (<STDIN>){
					chomp; last unless $_;
					last if (/^y$/i or /^yes$/i);
					if (/^n$/i or /^no$/i){ $user_okay = 0; last; }
				}
				last if $user_okay;
			}
			$computed_cutoff = $b[1][0];
		}
		
		# report cutoff to user
		
		if ($howCO == 0 or ($howCO >= 3 and $use_global)){ # user-defined global cutoff
			$fpN{$set}{$key}{'cutoff'} = $global_cutoff;
			print "  Cutoff is $global_cutoff (determined by global cutoff $globalCO.\n";
		}
		if ($howCO == 1 and $unite){ # user-defined individual cutoff
			if ($key eq '0'){ $s = $selfCO; }
			elsif ($key eq '1'){ $s = $closeCO; }
			elsif ($key eq '2'){ $s = $distalCO; }
			if ($s){
				$fpN{$set}{$key}{'cutoff'} = $s;
				print "  Cutoff is $s (user-defined).\n";
			}else{
				$fpN{$set}{$key}{'cutoff'} = $global_cutoff;
				print "  User-defined cutoff is not available. Use global cutoff $global_cutoff instead.\n";
			}
		}
		if ($howCO >= 3 and not $use_global){ # computed cutoff
			$s = 0; if ($howCO == 3){ $s = "histogram"; }
			elsif ($howCO == 4){ $s = "kernel density estimation"; }
			elsif ($howCO == 5){ $s = "clustering analysis"; }
			if ($computed_cutoff){
				$fpN{$set}{$key}{'cutoff'} = $computed_cutoff;
				print "  Cutoff is ".sprintf("%.3f", $computed_cutoff)." (determined by $s).\n";
			}else{
				$fpN{$set}{$key}{'cutoff'} = $global_cutoff;
				print "  ". ucfirst($s) ." failed to identify a cutoff. Use global cutoff $global_cutoff instead.\n";
			}
		}
		
		# ask user to enter cutoff
		
		my $user_enter;
		if ($interactive and $fpN{$set}{$key}{'cutoff'}){
			print "  Accept? (yes/no) ";
			while (<STDIN>){
				chomp;
				last unless $_;
				last if (/^y$/i or /^yes$/i);
				if (/^n$/i or /^no$/i){ $user_enter = 1; last; }
			}
		}
		if ($user_enter or $howCO == 2){
			print "  Enter user-specified cutoff: ";
			while (<STDIN>){
				chomp;
				last if /^\d+\.?\d*$/;
				print "  Invalid cutoff value. Re-enter: ";
			}
			$fpN{$set}{$key}{'cutoff'} = $_;
		}

		# @{$fpS{$set}{$key}{'data'}} = sort {$b<=>$a} @{$fpS{$set}{$key}{'data'}};
		# @a = @{$fpS{$set}{$key}{'data'}};
		# next unless @a;
		# $fpS{$set}{$key}{'n'} = @a;
		# $s = 0; $s += $_ for @a;
		# $fpS{$set}{$key}{'mean'} = sprintf("%.3f",$s/$fpS{$set}{$key}{'n'});
		# $s = 0; $s += ($fpS{$set}{$key}{'mean'} - $_)**2 for @a;
		# $fpS{$set}{$key}{'stdev'} = sprintf("%.3f",sqrt($s/($fpS{$set}{$key}{'n'}-1)));
		# $fpS{$set}{$key}{'min'} = $a[$#a];
		# $fpS{$set}{$key}{'max'} = $a[0];
		# $fpS{$set}{$key}{'median'} = median(@a);
		# $fpS{$set}{$key}{'mad'} = sprintf("%.3f",mad(@a));
		# @b = quantiles(@a);
		# $fpS{$set}{$key}{'q1'} = $b[1];
		# $fpS{$set}{$key}{'q3'} = $b[0];
		# $i = $fpS{$set}{$key}{'n'}*$globalCO;
		# if (int($i) == $i){
		# 	$fpS{$set}{$key}{'cutoff'} = ($a[$#a-$i+1]+$a[$#a-$i])/2;
		# }else{
		# 	if ($i-int($i) <= int($i)-$i+1){ $fpS{$set}{$key}{'cutoff'} = $a[$#a-int($i)+1]; }
		# 	else{ $fpS{$set}{$key}{'cutoff'} = $a[$#a-int($i)]; }
	}
}
if (($howCO == 4 and ($toolKDE or $toolExtrema)) or ($howCO == 5) or $dipTest){
	$R->stopR;
}
print " done.\n";


## output fingerprint ##

if ($outFp){
	open OUT, ">$wkDir/result/statistics/fingerprint.txt";
	print OUT "#NEXUS\nBEGIN STATISTICS;\n";
	print OUT "\tGroup\tNumber\tMean\tSD\tMin\tMax\tMedian\tMAD\tQ1\tQ3\tCutoff\n";
	foreach my $set (sort keys %fpN){
		foreach ('0','1','2'){
			%h = %{$fpN{$set}{$_}};
			if ($_ eq '0'){ $s = "self"; }elsif ($_ eq '1'){ $s = "close"; }elsif ($_ eq '2'){ $s = "distal"; }
			if ($set eq '0'){ print OUT "all"; }else{ print OUT $set; }
			print OUT "\t$s\t$h{'n'}\t".sprintf("%.2f", $h{'mean'})."\t".sprintf("%.2f", $h{'stdev'})."\t".sprintf("%.2f", $h{'min'})."\t".sprintf("%.2f", $h{'max'})."\t".sprintf("%.2f", $h{'median'})."\t".sprintf("%.2f", $h{'mad'})."\t".sprintf("%.2f", $h{'q1'})."\t".sprintf("%.2f", $h{'q3'})."\t".sprintf("%.2f", $h{'cutoff'})."\n";
		}
	}
	print OUT "END;\n";
	close OUT;
	print "Result is saved in result/statistics/fingerprint.txt.\n";
}

if ($interactive){
	print "Press Enter to proceed with prediction, or press Ctrl+C to exit:";
	$s = <STDIN>;
}


## conduct bidirectional best hit (BBH) search ##

if ($BBH){
	print "Conducting bidirectional best hit (BBH) search...";
	unless (-e "$wkDir/result/bbh.txt"){
		open OUT, ">$wkDir/result/bbh_input.txt";
		foreach my $set (keys %results){
			for ($i=0; $i<@{$results{$set}}; $i++){
				my %res = %{$results{$set}[$i]};
				if (exists $res{'BBH'} and $res{'BBH'}){
					# set id query_accn query_taxid subject_accn subject_taxid
					print OUT $set."\t".$i."\t".$res{'accn'}."\t".$selfinfo{$set}{'taxid'}."\t".$res{'hit2accn'}."\t".$res{'hit2'}."\n";
				}
			}
		}
		system "perl scripts/bbh.pl $wkDir";
	}
	open IN, "<$wkDir/result/bbh.txt";
	while (<IN>){
		s/\s+$//;
		@a = split (/\t/);
		$results{$a[0]}[$a[1]]{'BBH'} = "" if ($a[$#a] ne '1');
	}
	close IN;
	print " done.\n";
}


################################################
## Analyze the whole master record and output ##
################################################

mkdir "$wkDir/result/detail" unless -d "$wkDir/result/detail";

print "Predicting...";
foreach my $set (keys %results){
	$n = @{$results{$set}};
	for ($i=0; $i<$n; $i++){
		my %res = %{$results{$set}[$i]};
		next if ($minHits && ($res{'n'} < $minHits));
		next if ($minSize && ($res{'length'} < $minSize));

		my $fp = $set; 	$fp = '0' if $unite;

		# principle: fewer hits suggests income or loss
		# check if self hits are low (indicating the gene is not prevalent within the group)
		if ($res{'N0'} < $fpN{$fp}{'0'}{'cutoff'}){
			$res{'in'} = 1;
		}

		## predict incoming HGT events ##
		# hits from close sister groups are low, indicating there's no vertical ancestor
		if ($res{'N1'} < $fpN{$fp}{'1'}{'cutoff'}){
			# overall non-self hits are normal, indicating it's not an origination event
			# if ($res{'N1'}+$res{'N2'} >= $fpN{$set}{'12'}{'cutoff'}){ ##########################################!!!!!!!!!
			if ($res{'N2'} >= $fpN{$fp}{'2'}{'cutoff'}){
				$res{'income'} = 1;
				$res{'income'} = "" if ($selfLow and ($res{'N0'} >= $fpN{$fp}{'0'}{'cutoff'}));
			}
		}
		
		## predict gene loss events ##
		# hits from close sister groups are normal,
		elsif ($res{'in'}){
			# self hits are low, suggesting the gene is absent in some lineages
			$res{'loss'} = 1;
		}

		if (0){ # not available in this version
			## predict gene origination event ##
			# detection won't work for saturated blast report.
			if ($maxHits and ($res{'N0'}+$res{'N1'}+$res{'N2'}) < $maxHits){	
				# if every hit is self, then must be origination.
				unless ($res{'N1'}+$res{'N2'}){	
					$res{'origin'} = 1;							
				}else{
					# if non-self hits are significantly weak.
					$s = 0;
					$s = $res{'S1'}[0] if $res{'N1'};
					$s = $res{'S2'}[0] if ($res{'N2'} and $res{'S2'}[0] > $s);
					# if ($s < $fpS{$fp}{'1'}{'mean'}-2*$fpS{$fp}{'1'}{'stdev'}){
					# 	$res{'origin'} = 1;
					# }
				}
			}
		}

		# predicting gene outcoming events
		# self hits must be normal, excluding any paralogs
		unless ($res{'in'}){
			if ($res{'N1'}+$res{'N2'} and exists($res{'sGenus'}{'q1'})){
				if ($res{'nGenus'}{'max'} > $res{'sGenus'}{'q1'}){		# the best hit from other genera fall within the range of self genus hits
					@a = @{$res{'nGenus'}{'scores'}};
					$s = 0;
					foreach (@a){
						if ($_ < $res{'sGenus'}{'min'}){ last;}
						else{ $s = $_;}
					}
					if ($s > $res{'sGenus'}{'q1'}){
						$res{'genus_outcome'} = $res{'nGenusHit'};
					}
				}
			}
		}
		$results{$set}[$i] = {%res};
	}

	###### Output report ######

	open (OUT, ">$wkDir/result/detail/$set.txt");
	@a = ('Query','Length','Product','Hits','Self','Close','Distal','HGT');
	unless ($BBH){
		push (@a, "Loss") if $loss;
		push (@a, "POE") if $POE;
	}
	push @a, "Match";
	print OUT "HGTector result of $set\n".join("\t",@a)."\n";
	for ($i=0; $i<$n; $i++){
		%h = %{$results{$set}[$i]};
		print OUT $h{'query'}."\t".$h{'length'}."\t".$h{'product'}."\t".$h{'n'}."\t";
		print OUT "\n" and next if ($minHits and ($h{'n'} < $minHits)) or ($minSize and ($h{'length'} < $minSize)); #####??????#######
		print OUT sprintf("%.2f", $h{'N0'})."\t".sprintf("%.2f", $h{'N1'})."\t".sprintf("%.2f", $h{'N2'})."\t";
		unless ($BBH){
			print OUT $h{'income'};
			print OUT "\t".$h{'loss'} if $loss;
			print OUT "\t".$h{'origin'} if $POE;
		}else{
			print OUT $h{'BBH'} if exists $h{'BBH'} and $h{'BBH'};
		}
		print OUT "\t";
		if (exists $h{'hit2'} and exists $taxadb{$h{'hit2'}}){
			print OUT $h{'hit2'}." (".$taxadb{$h{'hit2'}}{'name'}.")";
		}
		print OUT "\n";
	}
	close OUT;
}
print " done.\n";
print "Prediction results are saved in result/detail/.\n";
exit 0;



##################
## sub routines ##
##################

sub median (@){
	my $mid = int(@_/2);
	return $_[0] unless $mid;
	if (@_ % 2) {
		return $_[$mid];
	}else{
		return ($_[$mid-1] + $_[$mid])/2;
	}
}

sub mad (@){
	my @absdev;
	my $median = median(@_);
	foreach my $x(@_){
		push @absdev, abs($x - $median);
	}
	return median(sort{$a<=>$b}@absdev);
}

sub quantiles (@){
	my $Q1, my $Q3;
	my $mid = int(@_/2);
	return ($_[0],$_[0]) unless $mid;
	if (@_ % 2) {
		$Q1 = median(@_[0..$mid]);
		$Q3 = median(@_[$mid..$#_]);
	}else{
		$Q1 = median(@_[0..($mid-1)]);
		$Q3 = median(@_[$mid..$#_]);
	}
	return ($Q1,$Q3);
}


# compute Z-score (Z = (xi-x^)/s)

sub z_scores(@){
	my $mean = 0;
	$mean += $_ for @_;
	$mean = $mean / @_;
	my $stdev = 0;
	$stdev += ($mean-$_)**2 for @_;
	$stdev = sqrt($stdev/(@_-1));
	my @z = ();
	push (@z, ($_-$mean)/$stdev) for @_;
	return @z;
}


# Z-score test for outliers (|Z| > 3)

sub z_test(@){
	my @data = sort{$a<=>$b}@_;
	my @z = z_scores(@data);
	my $lower_fence = $data[0];
	my $upper_fence = $data[$#data];
	for (my $i=0; $i<=$#data; $i++){
		if (abs($z[$i]) <= 3){
			$lower_fence = $data[$i];
			last;
		}
	}
	for (my $i=$#data; $i>=0; $i--){
		if (abs($z[$i]) <= 3){
			$upper_fence = $data[$i];
			last;
		}
	}
	return ($lower_fence, $upper_fence);
}

# modified Z-score test for outliers (|modified_Z| > 3.5) (Iglewicz and Hoaglin, 1993)

sub modified_z(@){
	my @data = sort{$a<=>$b}@_;
	my $lower_fence = $data[0];
	my $upper_fence = $data[$#data];
	my $median = median(@data);
	my $mad = mad(@data);
	return ($data[0],$data[$#data]) unless $mad;
	for (my $i=0; $i<=$#data; $i++){
		if (abs(0.6745*($data[$i]-$median)/$mad) <= 3.5){
			$lower_fence = $data[$i];
			last;
		}
	}
	for (my $i=$#data; $i>=0; $i--){
		if (abs(0.6745*($data[$i]-$median)/$mad) <= 3.5){
			$upper_fence = $data[$i];
			last;
		}
	}
	return ($lower_fence, $upper_fence);
}


# boxplot test for outliers

sub boxplot(@){
	my $lower_fence, my $upper_fence;
	my @data = sort{$a<=>$b}@_;
	my @Q = quantiles(@data);
	my $iqr = $Q[1]-$Q[0];
	my $f = 3*exp(10/@data);
	$lower_fence = $Q[0]-$f*$iqr;
	$upper_fence = $Q[1]+$f*$iqr;
	return ($lower_fence, $upper_fence);
}


# adjusted boxplot test for outliers

sub adjusted_boxplot(@){
	my $lower_fence, my $upper_fence;
	my @data = sort{$a<=>$b}@_;
	my $median = median(@data);
	my @lower; my @upper;
	foreach (@data){
		push (@lower, $_) if ($_ <= $median);
		push (@upper, $_) if ($_ >= $median);
	}
	my @kernel;
	foreach my $i (@lower){
		foreach my $j (@upper){
			next if ($i == $j);
			push @kernel, (($j-$median)-($median-$i))/($j-$i);
		}
	}
	my $mc = median(sort{$a<=>$b}@kernel);
	my @Q = quantiles(@data);
	my $iqr = $Q[1]-$Q[0];
	my $f = 1.5; # *exp(10/@data);
	if ($mc >= 0){
		$lower_fence = $Q[0]-$f*exp(-3.5*$mc)*$iqr;
		$upper_fence = $Q[1]+$f*exp(4*$mc)*$iqr;
	}else{
		$lower_fence = $Q[0]-$f*exp(-4*$mc)*$iqr;
		$upper_fence = $Q[1]+$f*exp(3.5*$mc)*$iqr;
	}
	return ($lower_fence, $upper_fence);
}

sub recurse_deOutlier(@){		# assume data are high to low
	my @data = @_;
	my @fences = adjusted_boxplot @data;
	for ($i=0; $i<=$#data; $i++){
		if ($data[$i] < $fences[0]){
			$i--;
			last;
		}
	}
	if ($i < $#data){
		@data = @data[0..$i];
		@data = recurse_deOutlier @data;
	}
	return @data;
}

sub recurse_Z(@){
	my @data = @_;
	if (modified_z(@data) > 3.5){
		pop @data;
		@data = recurse_Z(@data);
	}
	return @data;
}

# Mann-Whitney U test
# use Statistics::Test::WilcoxonRankSum;
# my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();


	## Compute the statistics for each subset of hits
	#foreach my $key (keys %scores){
	#	$result{$key}{'n'} = 0;
	#	@a = @{$scores{$key}};
	#	next unless @a;
	#	$result{$key}{'n'} = @a;
	#	next unless ($#a);
	#	$s = 0; $s += $_ for @a;
	#	$result{$key}{'mean'} = $s/$result{$key}{'n'};
	#	$s = 0; $s += ($result{$key}{'mean'} - $_)**2 for @a;
	#	$result{$key}{'stdev'} = sqrt($s/($result{$key}{'n'}-1));
	#	$result{$key}{'min'} = $a[$#a];
	#	$result{$key}{'max'} = $a[0];
	#
	#	$result{$key}{'median'} = median(@a);
	#	$result{$key}{'mad'} = mad(@a);
	#	@b = quantiles(@a);
	#	$result{$key}{'q1'} = $b[1];
	#	$result{$key}{'q3'} = $b[0];
	#	@{$result{$key}{'scores'}} = @a;
	#}









