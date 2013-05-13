#!/usr/bin/perl

use warnings;
use strict;
$| = 1;

print "

This script names each cluster of genes/proteins by the name of its member with the shortest name.

Usage:
	perl cog_namer.pl <working directory> <clusterring scheme file> <optional: output file>

Output:
	output.txt

" and exit unless @ARGV;

print "  Naming orthologous groups...";

## Program parameters ##

my $wkDir = $ARGV[0];								# working directory
my $infile = $ARGV[1];								# clusterring scheme file
my $outfile;										# output file
if ($#ARGV >= 2){ $outfile = $ARGV[2]; }
else{ $outfile = "output.txt"; }

my %proteins = ();									# all input proteins (accn -> id)


my @sets;											# protein sets (genomes)

if (-e "$wkDir/config.txt"){
	open IN, "<$wkDir/config.txt";
	while (<IN>){
		s/#.*$//; s/\s+$//g; s/^\s+//g; next unless $_;
		@sets = split (/,/, $1) if /^inSets=(.+)$/;
	}
}

unless (@sets){
	if (-d "$wkDir/result/detail"){
		opendir (DIR, "$wkDir/result/detail");
		@sets = grep(/\.txt$/,readdir(DIR));
		s/\.txt$// for @sets;
		close DIR;
	}else{
		opendir (DIR, "$wkDir/blast");
		@sets = readdir (DIR);
		close DIR;
	}
}

my %setAccns = (); # protein accns specific to set
foreach my $set (@sets){
	$setAccns{$set} = [];
	if (-e "$wkDir/result/detail/$set.txt"){
		open IN, "<$wkDir/result/detail/$set.txt";
		while (<IN>){
			next if /^HGTector/; next if /^Query/;
			s/\s+$//; next if /^#/; next unless $_;
			my @a = split (/\t/);
			$proteins{$a[0]} = $a[2] if $#a >= 2;
			push @{$setAccns{$set}}, $a[0];
		}
		close IN;
	}else{
		next unless -d "$wkDir/blast/$set";
		opendir (DIR, "$wkDir/blast/$set");
		my @blasts = grep(/\.bla$/,readdir(DIR));
		close DIR;
		foreach (@blasts){
			/^(.+)\.bla$/;
			my $accn = $1;
			push @{$setAccns{$set}}, $1;
			open IN, "<$wkDir/blast/$set/$_" or next;
			while (<IN>){
				s/\s+$//;
				last if /^;/;
				if (/^\tProduct=(.+);$/){
					$proteins{$accn} = $1;
					$proteins{$accn} =~ s/\s+$//;
					last;
				}
			}
			close IN;
		}
	}
	# print $set.": ".(scalar @{$setAccns{$set}})."\n";
}

my $i = 0;
open OUT, ">$outfile";
open IN, "<$infile";
while (<IN>){
	my @a = split (/\s+/); # blastclust output
	next unless @a;
	shift @a if $a[0] =~ /:$/; # orthomcl output
	next unless @a;
	$i ++;
	my $name;
	my @names = ();
	foreach (@a){
		s/^.+\|//;
		if (exists $proteins{$_} and $proteins{$_}){
			push (@names, $proteins{$_});
			delete $proteins{$_};
		}
	}
	if (@names){
		@names = sort {length($a) <=> length($b)}(@names);
		my @noHypo = ();
		foreach (@names){
			push (@noHypo, $_) unless (/hypothetical/ or /hypotethical/ or /hypothetcial/);
		}
		if (@noHypo){
			$name = "$i|$noHypo[0]";
		}else{
			$name = "$i|$names[0]";
		}
	}else{
		$name = "$i|unknown";
	}
	print OUT $name."\t".join ("\/", @a)."\n";
}
close IN;

foreach my $set (@sets){
	foreach my $accn (@{$setAccns{$set}}){
		if (exists $proteins{$accn}){
			print OUT (++$i)."|$proteins{$accn}\t$accn\n";
		}
	}
}

close OUT;

print " done.\n";










