#!/usr/bin/perl

# HGTector (version 0.1.6): Genome-wide detection of horizontal gene transfer based on Blast hit distribution statistics
# Copyright (C) 2013, Qiyun Zhu, Katharina Dittmar. All rights reserved.
# Licensed under BSD 2-clause license.

use warnings;
use strict;
$| = 1;

print "

|---------------------|
|   HGTector v0.1.6   |
|---------------------|

";


## program parameters ##

my $s;

my $wkDir = "working";
my $tasks = 0;										# whether user specified any tasks

my $blaster = 0;									# Step 1 - batch BLASTp searches
my $taxonomer = 0;									# Step 2 - identify taxonomic information of returning matches
my $purifier = 0;									# Step 3 - post-blast modifications of matches
my $analyzer = 0;									# Step 4 - predict HGT events as well as gene loss and origination events
my $summarizer = 0;									# Step 6 - summarize results and generate reports

my $blastMode = 0;									# BLAST mode (0: http, 1: local, 2: remote)
my $nBlaster = 1;									# run multiple Blasters simultaneously (default: 1) (0: number of sets)

my $nProtein = 0;									# total number of proteins to Blast

## working directory ##

if (@ARGV){
	foreach (@ARGV)	{
		if (/^-/){
			if (/^-b/i){ $blaster = 1; }
			elsif (/^-t/i){ $taxonomer = 1; }
			elsif (/^-p/i){ $purifier = 1; }
			elsif (/^-a/i){ $analyzer = 1; }
			elsif (/^-s/i){ $summarizer = 1; }
		}else{
			$wkDir = $_;
		}
	}
}else{
	print "No arguments specified. Take working/ as working directory.\n";
}
die "Error: invalid working directory $wkDir.\n" unless -d $wkDir;


## read configuration ##

if (-e "$wkDir/config.txt"){
	print "Reading configuration...";
	open IN, "<$wkDir/config.txt";
	while (<IN>){
		s/#.*$//; s/\s+$//; s/^\s+//; next unless $_;
		$blaster = 1 if /^blaster=1$/;
		$taxonomer = 1 if /^taxonomer=1$/;
		$purifier = 1 if /^purifier=1$/;
		$analyzer = 1 if /^analyzer=1$/;
		$summarizer = 1 if /^summarizer=1$/;
		$blastMode = $1 if /^blastMode=(\d)$/;
		$nBlaster = $1 if /^nBlaster=(\d+)$/;
	}
	print " done.\n";
}else{
	print "config.txt not found in working directory. Use default settings.\n";
}


## check status of tasks ##

$tasks = $blaster + $taxonomer + $purifier + $analyzer + $summarizer;

unless ($tasks){
	if (-d "$wkDir/blast"){
		print "Blast results seem to already exist. Running Blaster will resume a halted batch Blast task. If you want to start a new batch Blast task and overwrite the old results, you should manually delete the \"blast\" folder beforing running this program.\n";
		print "Type 1 to proceed with Blaster, or type 0 to skip:";
		$s = <STDIN>; chomp $s; $s = 0 unless $s;
		if ($s){
			$blaster = 1;
		}else{
			$s = "$wkDir/taxonomy";
			if (-d $s and -e "$s/taxa.db" and -e "$s/ranks.db" and -e "$s/self.info"){
				print "Taxonomic information seems to already exist. Running Taxonomer will overwrite the information.\n";
				print "Type 1 to re-run Taxonomer, or type 0 to skip:";
				$s = <STDIN>; chomp $s; $s = 0 unless $s;
				if ($s){
					$taxonomer = 1;
				}else{
					print "Purifier is an optional step, which will further process Blast results.\n";
					print "Type 1 to run Purifier, or type 0 to skip:";
					$s = <STDIN>; chomp $s; $s = 0 unless $s;
					if ($s){
						$purifier = 1;
					}else{
						if (-d "$wkDir/result" and -d "$wkDir/result/detail"){
							print "Prediction result seems to already exist. Running Analyzer will overwrite the result.\n";
							print "Type 1 to run Analyzer, or type 0 to skip:";
							$s = <STDIN>; chomp $s; $s = 0 unless $s;
							if ($s){
								$analyzer = 1;
							}else{
								if (-d "$wkDir/result/HGT"){
									print "Summary of result seems to already exist. Running Summarizer will overwrite the report.\n";
									print "Type 1 to run Summarizer, or type 0 to exit:";
									$s = <STDIN>; chomp $s; $s = 0 unless $s;
									if ($s){
										$summarizer = 1;
									}else{
										print "Nothing is done.\n";
										exit 0;
									}
								}else{
									$summarizer = 1;
								}
							}
						}else{
							$analyzer = 1;
						}
					}
				}
			}else{
				$taxonomer = 1;
			}
		}
	}else{
		$blaster = 1;
	}
}

## execute blaster in single or parallel mode ##

if ($blaster){
	print "Reading input...";
	my @sets; # sets of proteins (genes), aka. genomes
	opendir (DIR, "$wkDir/input");
	my @files = readdir(DIR);
	close DIR;
	while (1){
		foreach (@files){
			next if (/^\./);
			open IN, "<$wkDir/input/$_" or next;
			push @sets, $_;
			if (/\.gbk?$/i){ # GenBank file
				while(<IN>){ s/\s+$//; $nProtein ++ if /^\s+\/protein_id="(.+)"$/; }
			}else{ # protein summary
				while(<IN>){ s/\s+$//; next unless $_; next if /^\d/; $nProtein ++; }
			}
			close IN;
		}
		print " done. $nProtein proteins from ".@sets." set(s) to Blast.\n";
		die unless @sets;
	
		if ($nBlaster == 1 or $blastMode != 0){
			foreach (@sets){
				system "perl scripts/blaster.pl $wkDir $_";
			}
		}else{
			print "Running multiple Blasters simultaneously...\n";
			print "The status will be updated every ten minutes.\n";
			my %runs = (); # Blaster instances, key = set name, 0 - done, 1 - ongoing, not exist - not started yet
			while (1){
				if (@sets - keys (%runs)){ # start new Blaster instances
					foreach my $set (@sets){
						last if ($nBlaster > 1) and ((grep { $_ == 1} values %runs) >= $nBlaster); # do nothing if the max instances are reached.
						next if exists $runs{$set};
						my $command = "perl scripts/blaster.pl $wkDir $set";
						$set =~ s/\.[^.]+$//;
						print "Blaster of $set started.\n";
						$runs{$set} = 1;
						sleep 1;
						if (!fork()) {
							close(STDOUT);
							exec $command;
						}
					}
				}
				last unless grep { $_ == 1} values %runs; # all Blaster instances are completed.
				sleep 600;
				print "Blaster(s) running:";
				my $nDone = 0; my $nFailed = 0;
				foreach my $run (keys %runs){ # check if Blaster instances are completed.
					# next unless $runs{$run};
					open LOG, "<$wkDir/blast/$run.log";
					while (<LOG>){
						if (/\t\.\s+$/){ $nDone ++; }
						elsif (/\tfailed\.\s+$/){ $nFailed ++; }
						if (eof and /^Program ended/ and $runs{$run}){
							print " $run (completed)";
							$runs{$run} = 0;
						}
					}
					close LOG;
					print " $run" if $runs{$run};
				}
				print ".\nProteins: $nDone completed, $nFailed failed, ".($nProtein-$nDone-$nFailed)." left.\n";
			}
			print "\nAll Blasters completed.\n";
		}
		unless ($tasks){
			if (-d "$wkDir/blast"){
				print "Type 1 to proceed with Taxonomer, type 2 to re-run Blaster to fill any incomplete ones, or type 0 to exit.";
				$s = <STDIN>; chomp $s; $s = 0 unless $s;
				if ($s == 1){ $taxonomer = 1; last; }
				elsif ($s == 2){ next; }
				elsif ($s == 0){ exit 0; }
			}
		}else{
			last;
		}
	}
}


## execute other steps ##

if ($taxonomer){
	system "perl scripts/taxonomer.pl $wkDir";
	unless ($tasks){
		if (-d "$wkDir/taxonomy"){
			if (-e "$wkDir/taxonomy/taxa.db" and -e "$wkDir/taxonomy/ranks.db" and -e "$wkDir/taxonomy/self.info"){
				print "Type 1 to proceed with Purifier, type 2 to skip Purifier and proceed with Analyzer, or type 0 to exit.";
				$s = <STDIN>; chomp $s; $s = 0 unless $s;
				if ($s == 1){ $purifier = 1; }
				elsif ($s == 2){ $analyzer = 1; }
				elsif ($s == 0){ exit 0; }
			}
		}
	}
}

if ($purifier){
	system "perl scripts/purifier.pl $wkDir";
	unless ($tasks){
		print "Type 1 to proceed with Analyzer, or type 0 to exit.";
		$s = <STDIN>; chomp $s; $s = 0 unless $s;
		if ($s == 1){ $analyzer = 1; }
		elsif ($s == 0){ exit 0; }
	}
}

if ($analyzer){
	system "perl scripts/analyzer.pl $wkDir";
	unless ($tasks){
		if (-d "$wkDir/result" and -d "$wkDir/result/detail"){
			print "Type 1 to proceed with Summarizer, or type 0 to exit.";
			$s = <STDIN>; chomp $s; $s = 0 unless $s;
			if ($s == 1){ $summarizer = 1; }
			elsif ($s == 0){ exit 0; }
		}
	}
}

if ($summarizer){
	system "perl scripts/summarizer.pl $wkDir";
}

print "All steps completed.\n";

exit 0;


