#!/usr/bin/perl -w
$|++; #---turn on the auto flush for the progress bar
use strict;
use Math::Random; #---For generating normal distribution random numbers
#use File::Path; #---for removing tmp GNUPLOT dat files
#use List::Util 'shuffle';#--- for shuffling an array
#use POSIX qw(log10);
use Time::HiRes qw( time );
#use Getopt::Long;

######################################################################################################################################################
#
#	Description
#		This is a perl script to pick up reads (in a fastq file) which has a poly nucleotide track, e.g. polyA.
#
#	Input
#		--mode=						
#			fastGrab				all parameters will be ignored expect strndSpcfc and minTrgtBaseTrckLen, just grab and print the seq with a strench of target base;
#			multiTrackPerRead		will detect multiple target tracks within one read, and will extract the trimmed sequences as each track
#		--fastqPath=				path of the fastq file
#		--inDirPath=				path of the input directory;
#		--outDir=					output directory;
#		--trgtBase=					target base, default = A;
#		--minTrgtBaseTrckLen=		the minimum length of the target base track, default = 20;
#		--minNonTargetLen=			the minimum number of Non target base track bases in the output seq; default = 10;
#		--strndSpcfc=				"yes" or "no"; strand specific or not; if no, it will detected the poly base track in reverse complement too; 
#		--printNonRedudant=			"yes" or "no"; WARNING: may cause memory crash if the output sequence file is big; print non-redundant read or not as an extra-file; default = "no"
#		--checkTrgtbaseTrack		"yes" or "no"; check the quality of the target base track or not; default = no;
#		--checkLeftPriorTrgtbaseQual"yes" or "no"; check the quality of the base track immediately before the target base track;  default = no;
#		--minTrgtbaseTrackHiQScore	integer; minimum score of a base in the target base track that is regarded as high quality;
#		--minTrgtbaseTrackHiQLen	integer; minimum number of bases in the target base track that is regarded as high quality;
#		--minPriorTrgtbaseHiQScore	integer; minimum score of a base in the bases immediately before the target base track that is regarded as high quality;
#		--checkPriorTrgtbaseHiQLen 	integer; the number of bases immediately before the target base track to be as high quality;
#		--stdOut=					"yes" or "no"; ["no"]; if "yes", will print the result to standard out 
#
#	Output
#
#	Usage
#		
#		perl fastqQPolyBaseTrackPicker_v0.3.pl --fastqPath=/Volumes/CCHON1TBA/NGSData/fastq/run.Mar.2011/lib2ABC/lib2C_HM1_Rep3_pair_s8.extendedFrags.fastq --outDir=/Volumes/CCHON1TBA/NGSData/fastq/run.Mar.2011/polyA/lib2C_HM1_Rep3_pair_s8/ --strndSpcfc=no -outReadTag=HM1_Rep3
#		
#	Assumption
#		the quality values if referring to the sequence immediately prior to this quality string
#
#	Version history
#
#		Previous versions as fastqIDLengthNCleaner 		
#
#		v0.1
#			-debut;
#
#		v0.2
#			-minLength option added;
#			-will generate histogram and cumulative curve
#
#		v0.3
#			-name changed as fastqIDAndLengthTrimmer
#			-changed to read 4 lines at a time rather than putting into a seqID hash, hope it runs faster
#			-printEdit and printReject and printIDTable options added.
#			-a log file will be produced;
#
#		v0.4
#			-occur option added;
#
#		v0.5
#			-renamed to fastqIDLengthNCleaner
#			-added the randRepN option;
#			-added effective length statistics; effective length is defined as number of non-B qual flag base;
#
#		v0.6
#			-added option to filtered out sequences containing certain pattern;
#			-header option added;
#			-lowQFlag option added;
#
#		v0.7
#			-multiple line progress indicator changed to single-line;
#			-printfastGrabModefastGrabModeTrim5SortSeq option added;
#
#		v0.8
#			-fix a bug for wrong taken the fastGrabModefastGrabModeTrim5Seq from the trimmed but not the original sequence
#
#		v0.9
#			-added function to count AT content
#
#		v0.10
#			--trim3= option added
#			--use rejectCountHsh to count separately the lenght and HiQ rejections;
#
#		v0.11
#			--added outPolyARead option
#
#		Debut version of fastqQPolyBaseTrackPicker
#	
#		v0.1	debut
#
#		V0.2	
#			 -added strndSpcfcty option
#
#		V0.3
#			-added gapped track search, by allowing maxGapSize and maxGapNum	
#
#		v0.4
#			-added the maxNonTargetHiQLen option, aimed to pick up the 27nt with polyA tail
#			-added minNonTargetLen maxNonTargetLen;
#			-added simpleStrndSpcfcTailMode and printNonMatch options;
#			-added fastGrabMode option;
#			-added dirPath option, will loop over all fastqs in the directory, and will override the fastqPath option;
#
#		v0.5
#			-added the multiTrackPerRead mode; instead of getting the sequence from either side of on the leftmost target track per sequence, it'll get the flanking sequence from EVERY target track; e.g. ACTAGAGAGGACTAGCTAGAAAAAAAAAAAAACTAGCATCGACTACGATCACGTAAAAAAAAAAAAAAAA, both ACTAGAGAGGACTAGCTAG and ACTAGAGAGGACTAGCTAGAAAAAAAAAAAAACTAGCATCGACTACGATCACGT will be taken.
#			-with outTrimLen option added;
#			-added printNonRedudant option;
#			-in the multiTrackPerRead mode, a more thorough check for both strand will be done; also read contains > 50% homopolyer will not be taken;
#
#		v0.6
#			-added the fastGrabModeTrim5 option
#			-added autoCheck header function
#		
#		v0.7
#			-added autodetect quality score scale
#			-lowQFlag became obsolete
#			-mode simpleStrndSpcfcTail and gapExtending became obsolete
#			-all right sequence output were removed
#			
#		v0.8
#			-added --stdOut= option, allow printing the fastq directly to stdOut and suppress all other verbose output;
#
#
#####################################################################################################################################################

#==========================================================Main body starts==========================================================================#
use vars qw ($mode $fastqPathHsh_ref $outDir $trgtBase $minTrgtBaseTrckLen $minNonTargetLen $checkTrgtbaseTrack $checkLeftPriorTrgtbaseQual $minTrgtbaseTrackHiQScore $minTrgtbaseTrackHiQLen $minPriorTrgtbaseHiQScore $checkPriorTrgtbaseHiQLen $strndSpcfc $printNonRedudant $stdOut);
($mode, $fastqPathHsh_ref, $outDir, $trgtBase, $minTrgtBaseTrckLen, $minNonTargetLen, $checkTrgtbaseTrack, $checkLeftPriorTrgtbaseQual, $minTrgtbaseTrackHiQScore, $minTrgtbaseTrackHiQLen, $minPriorTrgtbaseHiQScore, $checkPriorTrgtbaseHiQLen, $strndSpcfc, $printNonRedudant, $stdOut) = readParameters();

printCMDLogOrFinishMessage("CMDLog");

my %fastqPathHsh = %{$fastqPathHsh_ref};
my $fastqNum = keys %fastqPathHsh;
die "No fastq input. Quitting" if ($fastqNum == 0);

print "$fastqNum files found\n" if $stdOut eq 'no';
my %tempFileNameHsh;
foreach my $fastqPath (sort {$a cmp $b} keys %fastqPathHsh) {
	my $fileName = $fastqPathHsh{$fastqPath};
	print $fileName."\n" if $stdOut eq 'no';
	die "$fileName is duplicated, please change filename before proceed\n" if (exists $tempFileNameHsh{$fileName});
	$tempFileNameHsh{$fileName}++;
}
print "\n" if $stdOut eq 'no';

#---define qualityScale
my($qualScaleHsh_ref, $qualPhred33Or64Hsh_ref, $scalePriorityHsh_ref) = defineQualityScale();

my $procNum = 0;
foreach my $fastqPath (sort {$a cmp $b} keys %fastqPathHsh) {
	$procNum++;
	print "Processing $procNum of $fastqNum, $fastqPathHsh{$fastqPath}.\n" if $stdOut eq 'no';
	my $printRight = "no";
	my ($header, $checkedPct) = autoCheckHeader($fastqPath, 3, 10000, 60, $stdOut);
	my ($qualCharHsh_ref) = getQualityCharFromFastq($fastqPath, 10000, $header, $stdOut);
	my ($validQualScale) = autoCheckQualScale($qualCharHsh_ref, $qualScaleHsh_ref, $qualPhred33Or64Hsh_ref, $scalePriorityHsh_ref, $stdOut);
	printEditedFileOnTheFly($fastqPath, $printRight, $validQualScale, $qualScaleHsh_ref, $header, $stdOut);
}

printCMDLogOrFinishMessage("finishMessage") if $stdOut eq 'no';

########################################################################## readParameters
sub readParameters {
	
	$mode = "multiTrackPerRead";
	my $fastqPath = "default";
	my $inDirPath = "default";
	$outDir = "default";
	$trgtBase = "A";
	$minTrgtBaseTrckLen = 5;
	$minNonTargetLen = 18;
	$checkTrgtbaseTrack = "no";
	$checkLeftPriorTrgtbaseQual = "no";
	$minTrgtbaseTrackHiQScore = 10;
	$minTrgtbaseTrackHiQLen = 5;
	$minPriorTrgtbaseHiQScore = 20;
	$checkPriorTrgtbaseHiQLen = 2;
	$strndSpcfc = "no";
	$printNonRedudant = "no";
	$stdOut = 'no';
	
	my (%fastqPathHsh, %inDirPathHsh);
	
	foreach my $param (@ARGV) {
		if ($param =~ m/--mode=/) {$mode = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--fastqPath=/) {$fastqPath = substr ($param, index ($param, "=")+1); $fastqPathHsh{$fastqPath}++;}
		elsif ($param =~ m/--inDirPath=/) {$inDirPath = substr ($param, index ($param, "=")+1); $inDirPath =~ s/\/$//; $inDirPathHsh{$inDirPath}++;}
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--trgtBase=/) {$trgtBase = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--minTrgtBaseTrckLen=/) {$minTrgtBaseTrckLen = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--minNonTargetLen=/) {$minNonTargetLen = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--checkTrgtbaseTrack=/) {$checkTrgtbaseTrack = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--checkLeftPriorTrgtbaseQual=/) {$checkLeftPriorTrgtbaseQual = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--minTrgtbaseTrackHiQScore=/) {$minTrgtbaseTrackHiQScore = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--minTrgtbaseTrackHiQLen=/) {$minTrgtbaseTrackHiQLen = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--minPriorTrgtbaseHiQScore=/) {$minPriorTrgtbaseHiQScore = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--checkPriorTrgtbaseHiQLen=/) {$checkPriorTrgtbaseHiQLen = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--strndSpcfc=/) {$strndSpcfc = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--printNonRedudant=/) {$printNonRedudant = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--stdOut=/) {$stdOut = substr ($param, index ($param, "=")+1);}
	}
	
	system ("mkdir -p -m 777 $outDir");
	system ("mkdir -p -m 777 $outDir/trimLeftTotal/");
	system ("mkdir -p -m 777 $outDir/trimLeftNR/");
	system ("mkdir -p -m 777 $outDir/fastGrab/");
	system ("mkdir -p -m 777 $outDir/log/");
	
	#--check the input
	if ($fastqPath ne "default") {
		foreach my $chkFastqPath (keys %fastqPathHsh) {
			my @chkFastqPathAry = split /\//, $chkFastqPath;
			my $fileName = $chkFastqPathAry[-1];
			open (TEST, "$chkFastqPath") || die "Can't open $chkFastqPath\n"; close TEST;
			$fastqPathHsh{$chkFastqPath} = $fileName;
		}
	}
		
	if ($inDirPath ne "default") {

		foreach my $chkInDirPath (keys %inDirPathHsh) {
			opendir (DIR, $chkInDirPath);
			my @dirContentAry = grep /\.fastq/, readdir DIR;
			closedir (DIR);
			
			foreach my $fileName (@dirContentAry) {
				my $chkFastqPath = $chkInDirPath."/".$fileName;
				open (TEST, "$chkFastqPath") || die "Can't open $chkFastqPath\n"; close TEST;
				$fastqPathHsh{$chkFastqPath} = $fileName;
			}
		}
	}
	
	return ($mode, \%fastqPathHsh, $outDir, $trgtBase, $minTrgtBaseTrckLen, $minNonTargetLen, $checkTrgtbaseTrack, $checkLeftPriorTrgtbaseQual, $minTrgtbaseTrackHiQScore, $minTrgtbaseTrackHiQLen, $minPriorTrgtbaseHiQScore, $checkPriorTrgtbaseHiQLen, $strndSpcfc, $printNonRedudant, $stdOut);
}
########################################################################## printEditedFileOnTheFly
sub printEditedFileOnTheFly {
	
	my $fastqPath = $_[0];
	my $printRight = $_[1];
	my $validQualScale = $_[2];
	my %qualScaleHsh = %{$_[3]};
	my $header = $_[4];
	my $stdOut = $_[5];
	
	#---get score info
	my $phred = ${$qualScaleHsh{$validQualScale}}{"phred"};
	
	my @fastqPathSplt = split /\//, $fastqPath;
	my $fastqName = $fastqPathSplt[-1];
	my @fastqNameSplit = split /\./, $fastqName;
	my $outReadTag = $fastqNameSplit[0];

	#--complement the target base for fastGrab
	my $compTrgtBase = $trgtBase;
	$compTrgtBase =~ tr/ACGTacgt/TGCAtgca/;

	if ($stdOut eq 'no') {
		if ($mode eq "fastGrab") {
			open (FASTGRABSEQ, ">$outDir/fastGrab/$fastqName.poly$trgtBase.min$minTrgtBaseTrckLen.fastGrab.fastq");
		} else {
			if ($printNonRedudant eq "yes") {
				open (NRTRIMMEDLEFTFQ, ">$outDir/trimLeftNR/$fastqName.poly$trgtBase.min$minTrgtBaseTrckLen.trimLeftNR.fastq");
			} else {
				open (TRIMMEDLEFTFQ, ">$outDir/trimLeftTotal/$fastqName.poly$trgtBase.min$minTrgtBaseTrckLen.trimLeftTotal.fastq");
			}
		}
	}
	
	my (%NRTrimLeftSeqHsh);
	
	#---get the total line number and determine the interval size

	#---define the start time and counters
	my $intervalStart = time();
	my $lineProc = my $progCount = my $acceptSeqNum = my $lenAndQRejectSeqNum = my $totalCount = my $filterSeqNum = 0;
	
	my $totalReadNum = my $readWithTrgtTrckNum =  my $leftPositiveReadNum = my $leftDiscardReadNum = my $trgtBaseDiscardTrackNum = my $trgtBaseHiQTrackNum = 0;

	my %polyBaseInfoHsh;
	
	my $uniqueID = 0;
	
	if ($fastqPath =~ m/\.gz$/) {
		open (INFILE, "gzip -dc $fastqPath |");
	} else {
		open (INFILE, "$fastqPath");
	}
	while (my $theLine = <INFILE>) {

		#---get the seq and ID
		if ($theLine =~ m/^\@$header/) {#---sequence header, assume 4 lines consists a read
			
			my $readWithTrgtTrck = "no";
			$totalReadNum++;

			my $seqHeader = $theLine; chomp $seqHeader; $seqHeader =~ s/^\@//;
			my $seq = <INFILE>; chomp $seq;
			
			my $seqLen = length $seq;
			
			my $qualHeader = <INFILE>; chomp $qualHeader; $qualHeader =~ s/^\+//;
			my $qual = <INFILE>; chomp $qual;
			
			#---check 4-lines-per-read format
			die "This fastq file doesnt seem to be in 4-lines-per-read format. Program terminated.\n" if ($qualHeader =~ m/^\+/);
			
			if ($mode eq "fastGrab") {

				#--will check both strand in any cases
				my %seqQualRevCompInfoHsh; 
				@{$seqQualRevCompInfoHsh{"oriDrtn"}} = ($seq, $qual);
				
				if ($strndSpcfc eq "no") {
					my $revCompSeq = reverse($seq);
					my $revCompQual = reverse($qual);
					$revCompSeq =~ tr/ACGTacgt/TGCAtgca/;
					@{$seqQualRevCompInfoHsh{"revComp"}} = ($revCompSeq, $revCompQual);
				}

				foreach my $directionToChk (keys %seqQualRevCompInfoHsh) {
					my ($seqToChk, $qualToChk) = @{$seqQualRevCompInfoHsh{$directionToChk}};
				
					if ($seqToChk =~ m/([$trgtBase]{$minTrgtBaseTrckLen,})/) {
						$readWithTrgtTrck = "yes";
						if ($stdOut eq 'no') {
							print FASTGRABSEQ "\@".$seqHeader."_".$directionToChk."\n";
							print FASTGRABSEQ $seqToChk."\n";
							print FASTGRABSEQ "\+\n";
							print FASTGRABSEQ $qualToChk."\n";
						} elsif ($stdOut eq 'yes') {
							print "\@".$seqHeader."_".$directionToChk."\n";
							print $seqToChk."\n";
							print "\+\n";
							print $qualToChk."\n";
						} else {
							die "stdOut option unknown\n";
						}
					}
				}
			
			} elsif ($mode eq "multiTrackPerRead") {

				#--will check both strand in any cases
				my %seqQualRevCompInfoHsh; 
				@{$seqQualRevCompInfoHsh{"oriDrtn"}} = ($seq, $qual);
				
				if ($strndSpcfc eq "no") {
					my $revCompSeq = reverse($seq);
					my $revCompQual = reverse($qual);
					$revCompSeq =~ tr/ACGTacgt/TGCAtgca/;
					@{$seqQualRevCompInfoHsh{"revComp"}} = ($revCompSeq, $revCompQual);
				}

				foreach my $directionToChk (keys %seqQualRevCompInfoHsh) {

					my ($seqToChk, $qualToChk) = @{$seqQualRevCompInfoHsh{$directionToChk}};
				
					if ($seqToChk =~ m/([$trgtBase]{$minTrgtBaseTrckLen,})/) {
						
						$readWithTrgtTrck = "yes";
						
						my $choppedSeq = $seqToChk;
						
						my (%trgtBaseTrackLenHsh, %trgtBaseTrackAbbSeqHsh);
						my $choppedPos = 0;
						
						#---collect the start pos and length of the target tracks
						while ($choppedSeq =~ m/([$trgtBase]{$minTrgtBaseTrckLen,})/) {
							
							my $trgtBaseTrackSeq = $1;
							my $trgtBaseTrackLen = length $trgtBaseTrackSeq;
							my $trgtBaseTrackStartOnChoppedSeq = index $choppedSeq, $trgtBaseTrackSeq;
							my $trgtBaseTrackStartOnSeq = $choppedPos + $trgtBaseTrackStartOnChoppedSeq;
	
							my $trgtBaseTrackLowQ = 0;
							if ($checkTrgtbaseTrack eq "yes") {
								my $trgtBaseTrackQual = substr $qualToChk, $trgtBaseTrackStartOnSeq, $trgtBaseTrackLen;
								my @trgtBaseTrackUnscaleScoreAry = unpack("C*", $trgtBaseTrackQual); #---returns ASCII code directly
								foreach my $unscaleScore (@trgtBaseTrackUnscaleScoreAry) {
									$trgtBaseTrackLowQ++ if (($unscaleScore - $phred) < $minTrgtbaseTrackHiQScore);
								}
							}
							
							my $trgtBaseTrackHiQ = $trgtBaseTrackLen - $trgtBaseTrackLowQ;
	
							if ($trgtBaseTrackHiQ >= $minTrgtbaseTrackHiQLen) {
								$trgtBaseHiQTrackNum++;
								my $abbTrgtBaseTrackSeq = abbreviateRepSeq($trgtBaseTrackSeq);
								$trgtBaseTrackLenHsh{$trgtBaseTrackStartOnSeq} = $trgtBaseTrackLen;
								$trgtBaseTrackAbbSeqHsh{$trgtBaseTrackStartOnSeq} = $abbTrgtBaseTrackSeq;
							} else {
								$trgtBaseDiscardTrackNum++;
							}
							
							$choppedPos = $trgtBaseTrackStartOnSeq+$trgtBaseTrackLen;
							$choppedSeq = substr $seqToChk, $choppedPos;
						}
						
						#----go through each target track
						foreach my $trgtBaseTrackStart (sort {$a <=> $b} keys %trgtBaseTrackLenHsh) {
							
							$uniqueID++;
							
							my $trgtBaseTrackLen = $trgtBaseTrackLenHsh{$trgtBaseTrackStart};
							my $abbTrgtBaseTrackSeq = $trgtBaseTrackAbbSeqHsh{$trgtBaseTrackStart};
							my $leftNonTrgtSeq = substr $seqToChk, 0, $trgtBaseTrackStart;
							my $leftNonTrgtQual = substr $qualToChk, 0, $trgtBaseTrackStart;
							my $leftNonTrgtLen = length $leftNonTrgtSeq;
							next if $leftNonTrgtLen < $minNonTargetLen;
							
							#----check the quality of the immediate upstream bases
							my $leftPriorTrgtbaseHiQ = "yes";
							
							if ($checkLeftPriorTrgtbaseQual eq "yes") {
								my $leftPriorTrgtbaseQualStartPos = $leftNonTrgtLen - $checkPriorTrgtbaseHiQLen;
								my $leftPriorTrgtbaseQualStr = substr $leftNonTrgtQual, $leftPriorTrgtbaseQualStartPos;
								my @leftPriorTrgtbaseUnscaleScoreAry = unpack("C*", $leftPriorTrgtbaseQualStr); #---returns ASCII code directly
								foreach my $unscaleScore (@leftPriorTrgtbaseUnscaleScoreAry) {
									if (($unscaleScore - $phred) < $minPriorTrgtbaseHiQScore) {
										$leftPriorTrgtbaseHiQ = "no";
										last;
									}
								}
							}
								
							if ($leftPriorTrgtbaseHiQ eq "yes") {
								
								$leftPositiveReadNum++;
								
								my $newSeqHeader = "HWI_".$directionToChk."_".$outReadTag."_".$uniqueID."_L_".$abbTrgtBaseTrackSeq;
				
								if ($printNonRedudant eq "yes") {
									@{$NRTrimLeftSeqHsh{$leftNonTrgtSeq}} = ($newSeqHeader, $leftNonTrgtQual);
								} else {
								
									if ($stdOut eq 'no') {
										print TRIMMEDLEFTFQ "\@".$newSeqHeader."\n";
										print TRIMMEDLEFTFQ $leftNonTrgtSeq."\n";
										print TRIMMEDLEFTFQ "\+\n";
										print TRIMMEDLEFTFQ $leftNonTrgtQual."\n";
									} elsif ($stdOut eq 'yes') {
										print "\@".$newSeqHeader."\n";
										print $leftNonTrgtSeq."\n";
										print "\+\n";
										print $leftNonTrgtQual."\n";
									} else {
										die "stdOut option unknown\n";
									}
								}
							} else {
								
								$leftDiscardReadNum++;
							}
						}
					}
				}
			
			} else {
			
				die "Unknown mode of search\n.";
			}
			
			$readWithTrgtTrckNum++ if $readWithTrgtTrck eq "yes";
		}
	}
	close INFILE;

	print "\n\n";

	if ($printNonRedudant eq "yes") {
		print "Print non-Redundant trimmed reads\n";
		foreach my $leftNonTrgtSeq (keys %NRTrimLeftSeqHsh) {
			my ($newSeqHeader, $leftNonTrgtQual) = @{$NRTrimLeftSeqHsh{$leftNonTrgtSeq}};
			if ($stdOut eq 'no') {
				print NRTRIMMEDLEFTFQ "\@".$newSeqHeader."\n";
				print NRTRIMMEDLEFTFQ $leftNonTrgtSeq."\n";
				print NRTRIMMEDLEFTFQ "\+\n";
				print NRTRIMMEDLEFTFQ $leftNonTrgtQual."\n";
			} elsif ($stdOut eq 'yes') {
				print "\@".$newSeqHeader."\n";
				print $leftNonTrgtSeq."\n";
				print "\+\n";
				print $leftNonTrgtQual."\n";
			} else {
				die "stdOut option unknown\n";
			}
		}
	}

	#----print the log file;
	open (LOGFILE, ">$outDir/log/$fastqName.fastqQPolyBaseTrackPicker.log.txt");
	print LOGFILE $0."\n\n";
	print LOGFILE "totalReadNum = $totalReadNum\n";
	print LOGFILE "readWithTrgtTrckNum = $readWithTrgtTrckNum\n";
	print LOGFILE "trgtBaseDiscardTrackNum = $trgtBaseDiscardTrackNum\n";
	print LOGFILE "trgtBaseHiQTrackNum = $trgtBaseHiQTrackNum\n";
	print LOGFILE "leftPositiveReadNum = $leftPositiveReadNum\n";
	print LOGFILE "leftDiscardReadNum = $leftDiscardReadNum\n";

	foreach my $stat (sort {$a cmp $b} keys %polyBaseInfoHsh) {
		print LOGFILE "\n".$stat."\n";
		no warnings 'numeric';
		foreach my $value (sort {$a <=> $b} keys %{$polyBaseInfoHsh{$stat}}) {
			my $count = ${$polyBaseInfoHsh{$stat}}{$value};
			print LOGFILE $value."\t".$count."\n";
		}
	}
	close LOGFILE;
}
########################################################################## checkFastqReadNumAndDefineIntervalSize
sub checkFastqReadNumAndDefineIntervalSize {
    
    my $fileToCheckPath = $_[0];
    my $linesToSample = $_[1];
	my $header = $_[2];
	my $stdOut = $_[3];
	
	#---make sure $linesToSample is a non-zero number, if not set to 10000
	my $linesToSampleInt = int $linesToSample;
	$linesToSampleInt = 100000 if (($linesToSampleInt != $linesToSample) or ($linesToSampleInt == 0));
	
	#---get the filename from the path
	my @fileToCheckPathSplt = split /\//, $fileToCheckPath;
	my $fileToCheckName = $fileToCheckPathSplt[-1];

	print "Estimating the number of lines in $fileToCheckName.\n";
	
	#---estimate the number of lines in the file
	open (INFILE, $fileToCheckPath) || die "Can't open $fileToCheckPath.\n";
	my $tmpFilePath = $fileToCheckPath."_tmp.txt";
	system "head -$linesToSampleInt $fileToCheckPath >$tmpFilePath";
	my $fileToCheckSize = -s "$fileToCheckPath";
	my $tmpFileSize = -s "$tmpFilePath";
	my $seqHeaderInTmp = `grep -e \@$header -c $tmpFilePath`;
	die "$fileToCheckName contains no fastq header. Quitting\n" if ((not (defined $seqHeaderInTmp)) or ($seqHeaderInTmp == 0));
	system "rm $tmpFilePath";
	my $fileToCheckSizeTotalReadNum = int (($fileToCheckSize/$tmpFileSize)*$seqHeaderInTmp);
	print "Estimated to have ".$fileToCheckSizeTotalReadNum." reads in $fileToCheckName.\n";

	my $intervalSize = int ($fileToCheckSizeTotalReadNum/100); #---define as
	$intervalSize = 1000000 if ($intervalSize > 1000000);
	$intervalSize = 1 if ($intervalSize == 0);
	return ($fileToCheckSizeTotalReadNum, $intervalSize);

}
########################################################################## reportProgress
sub reportProgress {

	my $progCount = $_[0];
	my $itemProc = $_[1];
	my $intervalSize = $_[2];
	my $fileTotalItemNum = $_[3];
	my $intervalStart = $_[4];
	my $item = $_[5]; #---reads or lines or alignments

	$progCount=0;
	my $intervalEnd = time();
	my $timeElapsed = $intervalEnd - $intervalStart;
	$timeElapsed = sprintf ("%.2f", $timeElapsed);
	my $estimatedEnd = (($fileTotalItemNum - $itemProc)*$timeElapsed)/$intervalSize;
	$estimatedEnd = sprintf ("%.2f", $estimatedEnd/60);
	print "$itemProc $item processed. Last $intervalSize $item:".$timeElapsed." sec. Estimated end: ".$estimatedEnd." mins.          \r";
	$intervalStart = time();
	
	return ($progCount, $intervalStart);
		
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
########################################################################## abbreviateRepSeq
sub abbreviateRepSeq {

	my $seq = $_[0];
	
	my @seqSplt = split "", $seq;
	
	my $abbStr = "";
	my $chrCount = 1;
	my $abbToAdd = "";

	if (length $seq == 1) {
	
		$abbStr = "[".$seq."x1]";
	
	} else {
		
		for (my $i = 0; $i < $#seqSplt; $i++) {
			my $curntChr = $seqSplt[$i];
			my $nextChr = $seqSplt[$i+1];
			if ($curntChr eq $nextChr) {
				$chrCount++;
				$abbToAdd = "[".$curntChr."x".$chrCount."]";
			} else {
				if ($chrCount>1) {
					$abbStr .=  $abbToAdd;
				} else {
					$abbStr .= $curntChr;
				}
				$chrCount = 1;
				$abbToAdd = $nextChr;
			}
		}
		
		$abbStr .=  $abbToAdd;
	}
	
	return $abbStr;

}
########################################################################## autoCheckHeader
sub autoCheckHeader {
	
	#---subroutine dependency: 
	#---in/out: my ($checkedHeader, $checkedPct) = autoCheckHeader($fastqPath, $checkWordSize, $checkHeaderNum, $minPct);
	
	my $fastqPath = $_[0];
	my $checkWordSize = $_[1];
	my $checkHeaderNum = $_[2];
	my $minPct = $_[3];
	my $stdOut = $_[4];
	
	my %checkHeaderHsh;
	my $sampleCount = 0;
	print "Auto-checking the header of fastq\n" if $stdOut eq 'no';

	if ($fastqPath =~ m/\.gz$/) {
		open (INFILE, "gzip -dc $fastqPath | ");
	} else {
		open (INFILE, "$fastqPath");
	}
	while (my $theLine = <INFILE>) {
		if ($theLine =~ m/^\@/) {#---sequence header, assume 4 lines consists a read
			my $header = substr $theLine, 1, $checkWordSize;
			$checkHeaderHsh{$header}++;
			$sampleCount++;
			last if ($sampleCount > $checkHeaderNum);
		}
	}
	close INFILE;
	
	my $checkedHeader = "";
	my $checkedPct;
	foreach my $header (sort {$checkHeaderHsh{$b} <=> $checkHeaderHsh{$a}} keys %checkHeaderHsh) {
		my $headerCount = $checkHeaderHsh{$header};
		$checkedPct = 100*$headerCount/$sampleCount;
		$checkedHeader = $header if ($checkedPct >= $minPct);
		last;
	}
	print "Header $checkedHeader sampled at $checkedPct%.\n" if $stdOut eq 'no';
	
	return ($checkedHeader, $checkedPct);
}
########################################################################## autoCheckHeader
sub getQualityCharFromFastq {
	
	#---subroutine dependency: 
	#---in/out: my ($qualCharHsh_ref) = getQualityCharFromFastq($fastqPath, $checkReadNum, $header);
	
	my $fastqPath = $_[0];
	my $checkReadNum = $_[1];
	my $header = $_[2];

	my $sampleCount = 0;
	my %qualCharHsh;
	
	if ($fastqPath =~ m/\.gz$/) {
		open (INFILE, "gzip -dc $fastqPath | ");
	} else {
		open (INFILE, "$fastqPath");
	}
	print "Getting quality character from fastq.\n" if ($stdOut eq 'no');
	while (my $theLine = <INFILE>) {

		chomp $theLine;
		#---get the seq and ID
		if ($theLine =~ m/^\@$header/) {#---sequence header, assume 4 lines consists a read
			
			$sampleCount++;
			
			my $seqHeader = $theLine;
			my $seq = <INFILE>;
			my $qualHeader = <INFILE>;
			die "fastq seems not to be in 4 lines per read format. Quitting\n" if not ($qualHeader =~ m/^\+/);
			my $qual = <INFILE>; chomp $qual;
			my @qualSplt = split //, $qual;
			foreach my $qualChar (@qualSplt) {
				$qualCharHsh{$qualChar}++;
			}
			last if $sampleCount >= $checkReadNum;
		}
	}
	close INFILE;

	return (\%qualCharHsh);
}
########################################################################## autoCheckHeader
sub autoCheckQualScale {
	
	#---subroutine dependency: 
	#---in/out: my ($validScale) = autoCheckQualScale($qualCharHsh_ref, $qualScaleHsh_ref, $qualPhred33Or64Hsh_ref, $scalePriorityHsh_ref);
	
	my %qualCharHsh = %{$_[0]};
	my %qualScaleHsh = %{$_[1]};
	my %qualPhred33Or64Hsh = %{$_[2]};
	my %scalePriorityHsh = %{$_[3]};
	my $stdOut = $_[4];

	my %phredScoreRngHsh;
	foreach my $phredScale ((64, 33)) {
		${$phredScoreRngHsh{$phredScale}}{"max"} = -999999;
		${$phredScoreRngHsh{$phredScale}}{"min"} = 999999;
	}

	foreach my $qualChar (keys %qualCharHsh) {
		my $ASCII = ord($qualChar);
		foreach my $phredScale ((64, 33)) {#---just record the max and min
			my $score = $ASCII - $phredScale;
			${$phredScoreRngHsh{$phredScale}}{"max"} = $score if $score > ${$phredScoreRngHsh{$phredScale}}{"max"};
			${$phredScoreRngHsh{$phredScale}}{"min"} = $score if $score < ${$phredScoreRngHsh{$phredScale}}{"min"};
		}
	}

	my %validScaleHsh;
	foreach my $phredScale ((64, 33)) {
		foreach my $scale (@{$qualPhred33Or64Hsh{$phredScale}}) {
			$validScaleHsh{$scale}++;
		}
	}
	
	foreach my $phredScale ((64, 33)) {
		my $phredMax = ${$phredScoreRngHsh{$phredScale}}{"max"};
		my $phredMin = ${$phredScoreRngHsh{$phredScale}}{"min"};
		foreach my $scale (@{$qualPhred33Or64Hsh{$phredScale}}) {
			my $min = ${$qualScaleHsh{$scale}}{"min"};
			my $max = ${$qualScaleHsh{$scale}}{"max"};
			if (($phredMax > $max) or ($phredMin < $min)) {
				delete $validScaleHsh{$scale};
			}
		}
	}
	
	my $validScaleNum = keys %validScaleHsh;
	die "qualityScale seems to be incorrect. Qutting\n" if $validScaleNum < 1;
	foreach my $scale (keys %validScaleHsh) {
		print $scale." is valid.\n" if ($stdOut eq 'no');
	}

	my $validQualScale;
	foreach my $scale (sort {$scalePriorityHsh{$a} <=> $scalePriorityHsh{$b}} keys %scalePriorityHsh) {
		if (exists $validScaleHsh{$scale}) {
			$validQualScale = $scale;
			print $scale." has been chosen.\n" if ($stdOut eq 'no');
			last;
		}
	}

	return ($validQualScale);
}
########################################################################## defineQualityScale
sub defineQualityScale {

#	my($qualScaleHsh_ref, $qualPhred33Or64Hsh_ref, $scalePriorityHsh_ref) = defineQualityScale();
#
#  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
#  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
#  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
#  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
#  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
#  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
#  |                         |    |        |                              |                     |
# 33                        59   64       73                            104                   126#
#
# S - Sanger        Phred+33,  raw reads typically (0, 40)
# X - Solexa        Solexa+64, raw reads typically (-5, 40)
# I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
# J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
#    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
#    (Note: See discussion above).
# L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
#
	my %qualScaleHsh;
	my %qualPhred33Or64Hsh; 
	my %scalePriorityHsh; #---priority of choosing a scale if more than one matched;
	
	my $scale = "sanger";
	$scalePriorityHsh{$scale} = 4;
	${$qualScaleHsh{$scale}}{"min"} = 0;
	${$qualScaleHsh{$scale}}{"max"} = 40;
	${$qualScaleHsh{$scale}}{"phred"} = 33;
	${$qualScaleHsh{$scale}}{"lowQFlag"} = "!";#---score = 0
	push @{$qualPhred33Or64Hsh{33}}, $scale;

	$scale = "illumina18";
	$scalePriorityHsh{$scale} = 3;
	${$qualScaleHsh{$scale}}{"min"} = 0;
	${$qualScaleHsh{$scale}}{"max"} = 41;
	${$qualScaleHsh{$scale}}{"phred"} = 33;
	${$qualScaleHsh{$scale}}{"lowQFlag"} = "!";#---score = 0
	push @{$qualPhred33Or64Hsh{33}}, $scale;

	$scale = "solexa";
	$scalePriorityHsh{$scale} = 5;
	${$qualScaleHsh{$scale}}{"min"} = -5;
	${$qualScaleHsh{$scale}}{"max"} = 40;
	${$qualScaleHsh{$scale}}{"phred"} = 64;
	${$qualScaleHsh{$scale}}{"lowQFlag"} = ";";#---score = -5
	push @{$qualPhred33Or64Hsh{64}}, $scale;
	
	$scale = "illumina13";
	$scalePriorityHsh{$scale} = 2;
	${$qualScaleHsh{$scale}}{"min"} = 0;
	${$qualScaleHsh{$scale}}{"max"} = 40;
	${$qualScaleHsh{$scale}}{"phred"} = 64;
	${$qualScaleHsh{$scale}}{"lowQFlag"} = "@";#---score = 0
	push @{$qualPhred33Or64Hsh{64}}, $scale;

	$scale = "illumina15";
	$scalePriorityHsh{$scale} = 1;
	${$qualScaleHsh{$scale}}{"min"} = 2;#---2 is B, flag quality;
	${$qualScaleHsh{$scale}}{"max"} = 41;
	${$qualScaleHsh{$scale}}{"phred"} = 64;
	${$qualScaleHsh{$scale}}{"lowQFlag"} = "B";#---score = 2
	push @{$qualPhred33Or64Hsh{64}}, $scale;
	
	return \%qualScaleHsh, \%qualPhred33Or64Hsh, \%scalePriorityHsh;
}
