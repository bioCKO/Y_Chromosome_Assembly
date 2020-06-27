#Purpose: Identify BACs that are concordant with the Y chromosome scaffold and merge sequence at gaps\n";
#Use:  identify_concordant_BACs_and_merge.pl <nucmer coordinates output file> <BAC fasta file> <raw Y chromosome fasta file from 3D-DNA> <Assembly from 3D-DNA> <Assembled Y chromosome fasta file from 3D-DNA>\n";

#use strict;
use warnings;
use Data::Dumper;

open (NUC, "$ARGV[0]") || die "Unable to open nucmer file; $!\n";
open (FASTA, "$ARGV[1]") || die "Unable to open YBAC sequence: $!\n";
open (CONTIGFA, "$ARGV[2]") || die "Unable to open raw chrom fasta file: $!\n";
open (ASSEMBLY, "$ARGV[3]") || die "Unable to open assembly file\n";
open (ASSEMBLYFA, "$ARGV[4]") || die "Unable to open contig fasta file: $!\n";

my %fasta;
my %assembly;
my %contigsToName;
my %nameToContigs;
my $contigLength;
my %contigFasta;
my @assemblyPositions;
my %usedBACs;
my %alignedPositions;
my %alignedGaps;
my @contigBreakpoints;
my $hicGap;
#Use these to determine which groups of bacs are aligned through quickmerge
my @orderedAssembly;
my %assemblyOrientation;

#separate out contigs from fasta file
my $header;
while(<FASTA>) {
	if ($_ !~ /^>/) {
		chomp($_);
		$yBACFasta{$header} .= $_;
	}
	else {
	 	chomp($_);
	 	my @subHeader = split(/\s+/, $_);
	 	$header = substr($subHeader[0], 1);
	 	$yBACFasta{$header} = "";
	}
}
while(<CONTIGFA>) {
	if ($_ !~ /^>/) {
		chomp($_);
		$contigFasta{$header} .= $_;
	}
	else {
	 	chomp($_);
	 	my @subHeader = split(/\s+/, $_);
	 	$header = substr($subHeader[0], 1);
	 	$contigFasta{$header} = "";
	}
}

while(<ASSEMBLYFA>) {
	if ($_ !~ /^>/) {
		chomp($_);
		$assembly{$header} .= $_;
	}
	else {
	 	chomp($_);
	 	my @subHeader = split(/\s+/, $_);
	 	$header = substr($subHeader[0], 1);
	 	$assembly{$header} = "";
	}
}
while(<ASSEMBLY>) {
	if ($_ =~ /^>/) {
		chomp($_);
		my @subHeader = split(/\s+/, $_);
	 	my $contigNumber = $subHeader[1];
	 	my $contigSize = $subHeader[2];
	 	my $contigName = substr($subHeader[0], 1);
	 	if ($contigName =~ /^hic_gap/) {
	 		$hicGap = $contigName;
	 	}
		$contigsToName{$contigNumber} = $contigName;
		$nameToContigs{$contigName} = $contigNumber;
		$contigLength{$contigNumber} = $contigSize;
	}
	else {
	 	chomp($_);
	 	my @splitLine = split(/\s+/, $_);
	 	my %tempAssembly;
	 	my @tempBreakpoints = (-499);
	 	my @tempContigAssembly;
	 	my $j = 1;
	 	foreach my $contig (@splitLine) {
	 		if ($contig =~ /^-/) {
	 			$contig = substr($contig, 1);
	 			$assemblyOrientation{$contig} = "-";
	 		}
	 		else {
	 			$assemblyOrientation{$contig} = "+"; 
	 		}
	 		for (my $i = 0; $i < $contigLength{$contig}; $i++) {
	 			$tempAssembly{$j} = $contig;
	 			$j++
	 		}
	 		next if ($contig == $nameToContigs{$hicGap});
	 		push(@tempContigAssembly, $contig);
	 		my $start = $tempBreakpoints[-1] + 500;
	 		my $end = $start + $contigLength{$contig};
	 		push(@tempBreakpoints, $start);
	 		push(@tempBreakpoints, $end);
	 	}
	 	push(@assemblyPositions, \%tempAssembly);
	 	push(@contigBreakpoints, \@tempBreakpoints);
	 	push(@orderedAssembly, \@tempContigAssembly);
	}
}

my %alignedBACStart;
my %alignedBACEnd;
my %bacStart;
my %bacEnd;
my %bacEndNotAligned;
my %bacLength;
my %bacAlignLength;
my %bacOverallPercentIdentity;
my %bacAlignOrientation;
my %bacsThatSpanGaps;
my %bacsThatExtendContigsStart;
my %bacsThatExtendContigsEnd;
my %bacsThatEncompassContigs;
my %bacAlignPositions;

while (<NUC>) {
	#determine how many BAC contigs are aligned at beginning and end (± 1000bp)
	#determine how many BAC contigs have concordant alignments based on overall known length
	#determine overall percent identity in alignments
	chomp($_);
	my @splitLine = split(/\t/, $_);
	my $bac = $splitLine[10];
	my $scaffold = $splitLine[9];
	my $lenQ = $splitLine[8];
	my $individualIdentity = $splitLine[6];
	my $alignLength = $splitLine[5];
	next if ($individualIdentity < 99);
	next if ($scaffold !~ /HiC_scaffold_1$/);
	next if ($alignLength < 10000);
	#new BAC alignment
	if (!exists($bacStart{$bac})) {
		#first alignment will include start position of BAC.
		my $s2 = $splitLine[2];
		my $e2  = $splitLine[3];
		for (my $i = $splitLine[0]; $i <= $splitLine[1]; $i++) {
			$bacAlignPositions{$bac}{$assemblyPositions[0]{$i}} = "";
		}
		if ($s2 < $e2) {
			$bacStart{$bac} = $s2;
			$bacEnd{$bac} = $e2;
			$alignedBACStart{$bac} = $splitLine[0];
			$alignedBACEnd{$bac} = $splitLine[1];
			$bacAlignOrientation{$bac} = "+";
		} 
		else {
			$bacStart{$bac} = $e2;
			$bacEnd{$bac} = $s2;
			$alignedBACStart{$bac} = $splitLine[1];
			$alignedBACEnd{$bac} = $splitLine[0];
			$bacAlignOrientation{$bac} = "-";
		}
		$bacLength{$bac} = $lenQ;
	}
	#remainder of the alignments. Update with new end positions if there are additional alignments
	else {
		my $s2 = $splitLine[2];
		my $e2  = $splitLine[3];
		for (my $i = $splitLine[0]; $i <= $splitLine[1]; $i++) {
			$bacAlignPositions{$bac}{$assemblyPositions[0]{$i}} = "";
		}
		if ($s2 < $e2) {
			$bacEnd{$bac} = $e2;
			if ($splitLine[0] < $alignedBACStart{$bac}) {
				$alignedBACStart{$bac} = $splitLine[0];
			}
			if ($splitLine[1] > $alignedBACEnd{$bac}) {
				$alignedBACEnd{$bac} = $splitLine[1];
			}
			if ($bacAlignOrientation{$bac} !~ /\+/) {
				$bacAlignOrientation{$bac} = 'mixed';
			}
		}
		else {
			$bacEnd{$bac} = $s2;
			if ($splitLine[1] > $alignedBACStart{$bac}) {
				$alignedBACStart{$bac} = $splitLine[1];
			}
			if ($splitLine[0] < $alignedBACEnd{$bac}) {
				$alignedBACEnd{$bac} = $splitLine[0];
			}
			if ($bacAlignOrientation{$bac} !~ /\-/) {
				$bacAlignOrientation{$bac} = 'mixed';
			}
		}
	}
	my $percentIdentity = $individualIdentity * $alignLength;
	if (!exists($bacOverallPercentIdentity{$bac})) {
		$bacOverallPercentIdentity{$bac} = $percentIdentity;
		$bacAlignLength{$bac} = $alignLength;
	}
	else {
		$bacOverallPercentIdentity{$bac} += $percentIdentity;
		$bacAlignLength{$bac} += $alignLength;
	}
}

#open (OUT, ">BAC_aligned.bed") || die "Unable to write BAC alignment bed file\n";
#print OUT "chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\tcomment\n";
my %alignedPacBioContigs;
my $concordantBAC = 0;
my $notAlignedBAC = 0;
my %notConcordantBoundaryBAC;
my %notCocordantEncompassBAC;
my $spansConcordantLarge = 0;
my %spansConcordantExtends;
my $rearrangements = 0;
my $duplications = 0;
my $notConcordant = 0;
my %spanPercentIdentity;
my %extendPercentIdentityStart;
my %extendPercentIdentityEnd;

# print Dumper(\%bacAlignPositions);
# exit(0);

open (BED, ">concordantBACs.bedTEST") || die "Unable to write output file:$!\n";
open (SUMMARY, ">manualAssembly.summary.txtTEST") || die "Unable to write output file: $!\n";
open (CONTIGLIST, ">contiglist.summary.txtTEST") || die "Unable to write contig list file: $!\n";

foreach my $bac (keys(%bacStart)) {
	my $startBoundary = "Start position NOT by boundary";
	my $endBoundary = "End position NOT by boundary";
	my $averagePercentIdentity = $bacOverallPercentIdentity{$bac} / $bacAlignLength{$bac};
	my $alignLength = abs($alignedBACEnd{$bac} - $alignedBACStart{$bac});
	my $bacAlignLength = abs($bacEnd{$bac} - $bacStart{$bac});
	my $concordant = 'Not concordant';
	$alignedPacBioContigs{$assemblyPositions[0]{$alignedBACEnd{$bac}}} = "";
	$alignedPacBioContigs{$assemblyPositions[0]{$alignedBACStart{$bac}}} = "";

	#determine all the contigs each BAC aligns to. BAC names are in assemblyPositions
	my @alignedContigs = (215);
	if ($alignedBACStart{$bac} < $alignedBACEnd{$bac}) {
		for (my $i = $alignedBACStart{$bac}; $i < $alignedBACEnd{$bac}; $i++) {
			my $tempAlignedContig = $assemblyPositions[0]{$i};
			if ($alignedContigs[-1] != $tempAlignedContig) {
				push(@alignedContigs, $tempAlignedContig);
			}
		}
	}
	else {
		for (my $i = $alignedBACEnd{$bac}; $i < $alignedBACStart{$bac}; $i++) {
			my $tempAlignedContig = $assemblyPositions[0]{$i};
			if ($alignedContigs[-1] != $tempAlignedContig) {
				push(@alignedContigs, $tempAlignedContig);
			}
		}
	}
	#remove 215
	print "Aligned Bacs: " . join(" ", keys(%{ $bacAlignPositions{$bac} })) . "\n";
	shift(@alignedContigs);
	
	if ($bacAlignOrientation{$bac} =~ /mixed/) {
		$concordant = 'Not concordant-mixed';
		print "$bac\t Mixed orientations\n";
		$rearrangements++;
		next;
	}
	#each BAC end must be almost completely aligned (within 1 kb). If multiple BACs span gaps at that contig, keep the one with the highest percent identity
	#if the BAC spans multiple contigs, check to see if it aligns at all to contigs in between
	if (($bacLength{$bac} - $alignLength) > 0 && $bacStart{$bac} == 1 && ($bacLength{$bac} - $bacEnd{$bac}) == 0 && @alignedContigs > 1 && @alignedContigs == (keys(%{ $bacAlignPositions{$bac} }) + (keys(%{ $bacAlignPositions{$bac} }) - 1))) {
		$concordant = 'Concordant-spans gap';
		$spansConcordantLarge++;
		if (@alignedContigs == 3) {
			my $startContig = $alignedContigs[0];
			if (exists($spanPercentIdentity{$startContig})) {
				if ($spanPercentIdentity{$startContig} < $averagePercentIdentity) {
					$bacsThatSpanGaps{$startContig} = $bac;
					$spanPercentIdentity{$startContig} = $averagePercentIdentity;
				}
			}
			else {
				$bacsThatSpanGaps{$startContig} = $bac;
				$spanPercentIdentity{$startContig} = $averagePercentIdentity;
			}
		}
		if (@alignedContigs > 3) {
			for (my $i = 0; $i < @alignedContigs - 1; $i += 2) {
				my $startContig = $alignedContigs[$i];
				if (exists($spanPercentIdentity{$startContig})) {
					if ($spanPercentIdentity{$startContig} < $averagePercentIdentity) {
						$bacsThatSpanGaps{$startContig} = $bac;
						$spanPercentIdentity{$startContig} = $averagePercentIdentity;
					}
				}
				else {
					$bacsThatSpanGaps{$startContig} = $bac;
					$spanPercentIdentity{$startContig} = $averagePercentIdentity;
				}
			}
		}
	}

	if (($bacLength{$bac} - $alignLength) > 0 && @alignedContigs == (keys(%{ $bacAlignPositions{$bac} }) + (keys(%{ $bacAlignPositions{$bac} }) - 1))) {
		#determine if BAC end positions fall close to a contig boundary. The other end of the BAC must be almost completely aligned (within 1000 bp)
		#If multiple BACs extend from contigs, keep the one with the highest percent identity
		foreach my $breakpoint (@{ $contigBreakpoints[0] }) {
			next if ($breakpoint < 0);
			if (abs($alignedBACStart{$bac} - $breakpoint) == 0 && abs($bacEnd{$bac} - $bacLength{$bac}) == 0 && $alignLength < $bacLength{$bac}) {
				$startBoundary = "Start position IS by boundary";
				my $startContig = $alignedContigs[0];
				$notConcordantBoundaryBAC{$bac} = "";
				if (exists($extendPercentIdentityStart{$startContig})) {
					#print "Previous percent identity: $extendPercentIdentityStart{$startContig}\n";
					#print "New percent identity: $averagePercentIdentity\n";
					if ($extendPercentIdentityStart{$startContig} < $averagePercentIdentity) {
						$bacsThatExtendContigsStart{$startContig} = $bac;
						$extendPercentIdentityStart{$startContig} = $averagePercentIdentity;
					}
				}
				else {
					$bacsThatExtendContigsStart{$startContig} = $bac;
					$notConcordantBoundaryBAC{$bac} = "";
					$extendPercentIdentityStart{$startContig} = $averagePercentIdentity;
				}
				if ($concordant =~ /Concordant-spans gap/) {
					$concordant = 'Concordant-spans gap-extends';
					$spansConcordantExtends{$bac} = "";
				}
				else {
					$concordant = 'Concordant-extends';
				}
			}
			elsif (abs($alignedBACEnd{$bac} - $breakpoint) == 0 && $bacStart{$bac} == 1 && $alignLength < $bacLength{$bac}) {
				$endBoundary = "End position IS by boundary";
				my $startContig = $alignedContigs[-1];
				$notConcordantBoundaryBAC{$bac} = "";
				if (exists($extendPercentIdentityEnd{$startContig})) {
					if ($extendPercentIdentityEnd{$startContig} < $averagePercentIdentity) {
						$bacsThatExtendContigsEnd{$startContig} = $bac;
						$extendPercentIdentityEnd{$startContig} = $averagePercentIdentity;
					}
				}
				else {
					$bacsThatExtendContigsEnd{$startContig} = $bac;
					$notConcordantBoundaryBAC{$bac} = "";
					$extendPercentIdentityEnd{$startContig} = $averagePercentIdentity;
				}
				if ($concordant =~ /Concordant-spans gap/) {
					$concordant = 'Concordant-spans gap-extends';
					$spansConcordantExtends{$bac} = "";
				}
				else {
					$concordant = 'Concordant-extends';
				}
			}
			elsif (abs($alignedBACEnd{$bac} - $breakpoint) == 0 && ($alignedBACStart{$bac} - $breakpoint) == 0 && $alignLength < $bacLength{$bac}) {
				print "ENTERED encompass\n";
				exit(0);
				$startBoundary = "Start position IS by boundary";
				$endBoundary = "End position IS by boundary";
				$notConcordantEncompassBAC{$bac} = "";
				if ($concordant =~ /Concordant-spans gap/) {
					$concordant = 'Concordant-spans gap-encompass';
				}
				else {
					$concordant = 'Concordant-encompass';
				}
			}
		}
	}
	if (abs($alignLength - $bacLength{$bac}) <= 10000 && $concordant =~ /Not concordant/ && @alignedContigs == (keys(%{ $bacAlignPositions{$bac} }) + (keys(%{ $bacAlignPositions{$bac} }) - 1))) {
		$concordant = 'Concordant';
		$concordantBAC++;
	}
	
	print "$bac\t$bacStart{$bac}\t$bacEnd{$bac}\t$bacLength{$bac}\t$alignedBACStart{$bac}\t$alignedBACEnd{$bac}\t$alignLength\t$averagePercentIdentity\t$concordant\t";
	print "$startBoundary\t$endBoundary\t@alignedContigs\n";
	
	print SUMMARY "$bac\t$bacStart{$bac}\t$bacEnd{$bac}\t$bacLength{$bac}\t$alignedBACStart{$bac}\t$alignedBACEnd{$bac}\t$alignLength\t$averagePercentIdentity\t$concordant\t";
	print SUMMARY "$startBoundary\t$endBoundary\t@alignedContigs\n";
	if ($concordant !~ /Not concordant/) {
		foreach my $contig (@alignedContigs) {
			next if $contigsToName{$contig} =~ /$hicGap/;
			print CONTIGLIST "$contigsToName{$contig}\n";
		}
	}
	if ($concordant =~ /^Concordant/) {
		print BED "assembly\t$alignedBACStart{$bac}\t$alignedBACEnd{$bac}\t$bac\_$concordant\n";
		if (@alignedContigs > 1) {
			my @sortedContigs = sort {$b <=> $a} @alignedContigs;
			my $temp = join("\t", @sortedContigs);
			$alignedGaps{$temp} = '';
		}
		if ($alignedBACStart{$bac} < $alignedBACEnd{$bac}) {
			for (my $i = $alignedBACStart{$bac}; $i <= $alignedBACEnd{$bac}; $i++) {
				$alignedPositions{$i} = '';
			}
		}
		else {
			for (my $i = $alignedBACEnd{$bac}; $i <= $alignedBACStart{$bac}; $i++) {
				$alignedPositions{$i} = '';
			}
		}
	}
}
foreach my $bac (keys(%yBACFasta)) {
	if (!exists($bacStart{$bac})) {
		print "Not aligned: $bac\n";
		$notAlignedBAC++;
	}
}
my $notConcordantBoundaryBACCount = scalar(keys(%notConcordantBoundaryBAC));
my $notConcordantEncompassBACCount = scalar(keys(%notConcordantEncompassBAC));
my $spansConcordantExtendsCount = scalar(keys(%spansConcordantExtends));

print "Concordant: $concordantBAC\nNot concordant Boundary: $notConcordantBoundaryBACCount\nNot concordant Encompass: $notConcordantEncompassBACCount\nSpans Large Gap concordant: $spansConcordantLarge\nSpans Gap and Extends: $spansConcordantExtendsCount\nRearrangement: $rearrangements\nDuplication: $duplications\nNot Aligned: $notAlignedBAC\n";
my $total = $concordantBAC + $notConcordantBoundaryBACCount + $notConcordantEncompassBACCount + $spansConcordantLarge - $spansConcordantExtendsCount + $rearrangements + $duplications + $notAlignedBAC;
print "Total: $total\n";

my $alignedPositions = scalar(keys(%alignedPositions));
foreach my $gap (keys(%alignedGaps)) {
	my @split = split(/\t/, $gap);
	print "@split\n";
	foreach my $contig (@split) {
		if ($contig == 406) {
			print "GAP\n";
			$alignedPositions -= 500;
		}
	}
}

print "Total aligned positions: $alignedPositions\n";
 
#iterate through assembly and use quickmerge to merge bacs with 
#nucmer --mum -l 100 --prefix test arrowCombined.fa fastaSplit.259I16_seg2.fa
#delta-filter -r -q -l 10000 test.delta > test.filtered.delta
#~/quickmerge/quickmerge -d test.filtered.delta -q fastaSplit.259I16_seg2.fa -r arrowCombined.fa -hco 5.0 -c 1.5 -l 25000 -ml 10000 -p quickmerge.test

#first merge BACs that extend from contigs
my %extendedPacBio;
for (my $i = 0; $i < @{ $orderedAssembly[0] }; $i++) {
	if (exists($bacsThatExtendContigsStart{$orderedAssembly[0][$i]})) {
		print "START EXTEND\n";
		my $mergedFasta = $contigFasta{$contigsToName{$orderedAssembly[0][$i]}};
		#open nucmer files
		open (BACNUCMER, ">nucmer.bac.fa") || die "Unable to write nucmer Bac file: $!\n";
		open (PACBIONUCMER, ">nucmer.pacbio.fa") || die "Unable to write pacbio file: $!\n";
		print PACBIONUCMER ">$contigsToName{$orderedAssembly[0][$i]}.$orderedAssembly[0][$i]\n$mergedFasta\n";
		print BACNUCMER ">$bacsThatExtendContigsStart{$orderedAssembly[0][$i]}\n$yBACFasta{$bacsThatExtendContigsStart{$orderedAssembly[0][$i]}}\n";
		print "Before, bac $bacsThatExtendContigsStart{$orderedAssembly[0][$i]}: " . length($yBACFasta{$bacsThatExtendContigsStart{$orderedAssembly[0][$i]}}) . "\n";
		print "Before, pacbio $contigsToName{$orderedAssembly[0][$i]}.$orderedAssembly[0][$i]: " . length($mergedFasta) . "\n";
		
		#run nucmer, delta-filter, and quickmerge
		system("nucmer -l 100 --prefix temp nucmer.pacbio.fa nucmer.bac.fa");
		system("delta-filter -r -q -l 10000 -i 99 temp.delta > temp.filtered.delta");
		system("show-coords -lTH temp.filtered.delta > temp.filtered.delta.coords");
		
		#determine if they aligned. If not, keep the pacbio sequence
		open (COORDS, "temp.filtered.delta.coords") || die "Unable to open show-coords output file: $!\n";
		
		my $bacLength;
		my $pacbioLength;
		my $orientation;
		my @bacAlignPositions;
		my @pacbioAlignPositions;
		
		while (<COORDS>) {
			my @splitLine = split(/\t/, $_);
			print "@splitLine\n";
			$bacLength = $splitLine[8];
			$pacbioLength = $splitLine[7];
			if(!defined($orientation) && $splitLine[2] < $splitLine[3]) {
				$orientation = 'forward';
			}
			if(!defined($orientation) && $splitLine[2] > $splitLine[3]) {
				$orientation = 'reverse';
			}
			push(@bacAlignPositions, $splitLine[2]);
			push(@bacAlignPositions, $splitLine[3]);
			push(@pacbioAlignPositions, $splitLine[0]);
			push(@pacbioAlignPositions, $splitLine[1]);
		}	
		my @sortedBACAlignPositions = sort {$a <=> $b} @bacAlignPositions;
		my @sortedPacbioAlignPositions = sort {$a <=> $b} @pacbioAlignPositions;
		
		my $tempBacSequence;
		if ($bacLength == ($sortedBACAlignPositions[-1] - $sortedBACAlignPositions[0] + 1)) {
			print "Entire BAC aligned: NEXT!\n";
			next;
		}
		if (abs($sortedBACAlignPositions[-1] - $bacLength) == 0 && $orientation =~ /reverse/) {
			$tempBACSequence = substr($yBACFasta{$bacsThatExtendContigsStart{$orderedAssembly[0][$i]}}, 0, $sortedBACAlignPositions[0]);
			$tempBACSequence = reverse($tempBACSequence);
			$tempBACSequence =~ tr/ACGTacgt/TGCAtgca/;
		}
		elsif (abs($sortedBACAlignPositions[-1] - $bacLength) == 0 && $orientation =~ /forward/) {
			$tempBACSequence = substr($yBACFasta{$bacsThatExtendContigsStart{$orderedAssembly[0][$i]}}, 0, $sortedBACAlignPositions[0]);
		}
		elsif ($sortedBACAlignPositions[0] == 1 && $orientation =~ /reverse/) {
			$tempBACSequence = substr($yBACFasta{$bacsThatExtendContigsEnd{$orderedAssembly[0][$i]}}, $sortedBACAlignPositions[-1]);
			$tempBACSequence = reverse($tempBACSequence);
			$tempBACSequence =~ tr/ACGTacgt/TGCAtgca/;
		}
		elsif ($sortedBACAlignPositions[0] == 1 && $orientation =~ /forward/) {
			$tempBACSequence = substr($yBACFasta{$bacsThatExtendContigsEnd{$orderedAssembly[0][$i]}}, $sortedBACAlignPositions[-1]);
		}
		#if at start of BAC
		if ($sortedPacbioAlignPositions[0] == 1) {
			$mergedFasta = $tempBACSequence . $contigFasta{$contigsToName{$orderedAssembly[0][$i]}};
		}
		#if at end of BAC
		elsif ($sortedPacbioAlignPositions[-1] == $pacbioLength) {
			$mergedFasta = $contigFasta{$contigsToName{$orderedAssembly[0][$i]}} . $tempBACSequence;
		}
		else {
			print "alignment not at pacbio end\n";
			next;
		}
		my $length = length($mergedFasta);
		print "After merged: $length\n";
		
		unlink('nucmer.pacbio.fa', 'nucmer.bac.fa', 'temp.filtered.delta', 'temp.delta', 'temp.filtered.delta.coords');
		$extendedPacBio{$contigsToName{$orderedAssembly[0][$i]}} = $mergedFasta;
		
	}
	if (exists($bacsThatExtendContigsEnd{$orderedAssembly[0][$i]})) {
		my $mergedFasta;
		print "END EXTEND\n";
		if (exists($extendedPacBio{$contigsToName{$orderedAssembly[0][$i]}})) {
			$mergedFasta = $extendedPacBio{$contigsToName{$orderedAssembly[0][$i]}};
		}
		else {
			$mergedFasta = $contigFasta{$contigsToName{$orderedAssembly[0][$i]}};
		}
		
		#open nucmer files
		open (BACNUCMER, ">nucmer.bac.fa") || die "Unable to write nucmer Bac file: $!\n";
		open (PACBIONUCMER, ">nucmer.pacbio.fa") || die "Unable to write pacbio file: $!\n";
		print PACBIONUCMER ">$contigsToName{$orderedAssembly[0][$i]}.$orderedAssembly[0][$i]\n$mergedFasta\n";
		print BACNUCMER ">$bacsThatExtendContigsEnd{$orderedAssembly[0][$i]}\n$yBACFasta{$bacsThatExtendContigsEnd{$orderedAssembly[0][$i]}}\n";
		print "Before, bac $bacsThatExtendContigsEnd{$orderedAssembly[0][$i]}: " . length($yBACFasta{$bacsThatExtendContigsEnd{$orderedAssembly[0][$i]}}) . "\n";
		print "Before, pacbio $contigsToName{$orderedAssembly[0][$i]}.$orderedAssembly[0][$i]: " . length($mergedFasta) . "\n";
		
		#run nucmer, delta-filter, and quickmerge
		system("nucmer -l 100 --prefix temp nucmer.pacbio.fa nucmer.bac.fa");
		system("delta-filter -r -q -l 10000 -i 99 temp.delta > temp.filtered.delta");
		system("show-coords -lTH temp.filtered.delta > temp.filtered.delta.coords");
		
		#determine if they aligned. If not, keep the pacbio sequence
		open (COORDS, "temp.filtered.delta.coords") || die "Unable to open show-coords output file: $!\n";
		
		my $bacLength;
		my $pacbioLength;
		my $orientation;
		my @bacAlignPositions;
		my @pacbioAlignPositions;
		
		while (<COORDS>) {
			my @splitLine = split(/\t/, $_);
			print "@splitLine\n";
			$bacLength = $splitLine[8];
			$pacbioLength = $splitLine[7];
			if(!defined($orientation) && $splitLine[2] < $splitLine[3]) {
				$orientation = 'forward';
			}
			if(!defined($orientation) && $splitLine[2] > $splitLine[3]) {
				$orientation = 'reverse';
			}
			push(@bacAlignPositions, $splitLine[2]);
			push(@bacAlignPositions, $splitLine[3]);
			push(@pacbioAlignPositions, $splitLine[0]);
			push(@pacbioAlignPositions, $splitLine[1]);
		}	
		my @sortedBACAlignPositions = sort {$a <=> $b} @bacAlignPositions;
		my @sortedPacbioAlignPositions = sort {$a <=> $b} @pacbioAlignPositions;
		
		my $tempBacSequence;
		if ($bacLength == ($sortedBACAlignPositions[-1] - $sortedBACAlignPositions[0] + 1)) {
			print "Entire BAC aligned: NEXT!\n";
			next;
		}
		if (abs($sortedBACAlignPositions[-1] - $bacLength) == 0 && $orientation =~ /reverse/) {
			$tempBACSequence = substr($yBACFasta{$bacsThatExtendContigsStart{$orderedAssembly[0][$i]}}, 0, $sortedBACAlignPositions[0]);
			$tempBACSequence = reverse($tempBACSequence);
			$tempBACSequence =~ tr/ACGTacgt/TGCAtgca/;
		}
		elsif (abs($sortedBACAlignPositions[-1] - $bacLength) == 0 && $orientation =~ /forward/) {
			$tempBACSequence = substr($yBACFasta{$bacsThatExtendContigsStart{$orderedAssembly[0][$i]}}, 0, $sortedBACAlignPositions[0]);
		}
		elsif ($sortedBACAlignPositions[0] == 1 && $orientation =~ /reverse/) {
			$tempBACSequence = substr($yBACFasta{$bacsThatExtendContigsEnd{$orderedAssembly[0][$i]}}, $sortedBACAlignPositions[-1]);
			$tempBACSequence = reverse($tempBACSequence);
			$tempBACSequence =~ tr/ACGTacgt/TGCAtgca/;
		}
		elsif ($sortedBACAlignPositions[0] == 1 && $orientation =~ /forward/) {
			$tempBACSequence = substr($yBACFasta{$bacsThatExtendContigsEnd{$orderedAssembly[0][$i]}}, $sortedBACAlignPositions[-1]);
		}
		#if at start of contig
		if ($sortedPacbioAlignPositions[0] == 1) {
			$mergedFasta = $tempBACSequence . $mergedFasta;
		}
		#if at end of contig
		elsif ($sortedPacbioAlignPositions[-1] == $pacbioLength) {
			$mergedFasta = $mergedFasta . $tempBACSequence;
		}
		else {
			print "alignment not at pacbio end\n";
			next;
		}
		my $length = length($mergedFasta);
		print "After merged: $length\n";
		
		unlink('nucmer.pacbio.fa', 'nucmer.bac.fa', 'temp.filtered.delta', 'temp.delta', 'temp.filtered.delta.coords');
		$extendedPacBio{$contigsToName{$orderedAssembly[0][$i]}} = $mergedFasta;
	}
	if (!exists($extendedPacBio{$contigsToName{$orderedAssembly[0][$i]}})) {
		$extendedPacBio{$contigsToName{$orderedAssembly[0][$i]}} = $contigFasta{$contigsToName{$orderedAssembly[0][$i]}};
	}
}

#join contigs that are spanned by a bac
my %joinedContigs;
my %joinedContigsComposition;
my @finalAssembly;
for (my $i = 0; $i < @{ $orderedAssembly[0] }; $i++) {
	my $j = $i + 1;
	my %tempAssembliesSequence;
	my %tempAssembliesCounts;
	if (exists($bacsThatSpanGaps{$orderedAssembly[0][$i]})) {
		#open nucmer files
		open (BACNUCMER, ">nucmer.bac.fa") || die "Unable to write nucmer Bac file: $!\n";
		open (PACBIONUCMER, ">nucmer.pacbio.fa") || die "Unable to write pacbio file: $!\n";
		my $pacbio1Sequence;
		my $pacbio2Sequence;
		if (exists($joinedContigs{$contigsToName{$orderedAssembly[0][$i]}})) {
			$pacbio1Sequence = $joinedContigs{$contigsToName{$orderedAssembly[0][$i]}};
			print "previous pacbio length: " . length($pacbio1Sequence) . "\n";
			print PACBIONUCMER ">$contigsToName{$orderedAssembly[0][$i]}.$orderedAssembly[0][$i]\n$pacbio1Sequence\n";
		}
		else {
			$pacbio1Sequence = $extendedPacBio{$contigsToName{$orderedAssembly[0][$i]}};
			if ($assemblyOrientation{$orderedAssembly[0][$i]} =~ /-/) {
				$pacbio1Sequence = reverse($pacbio1Sequence);
				$pacbio1Sequence =~ tr/ACGTacgt/TGCAtgca/;
			}
			print "new pacbio length $contigsToName{$orderedAssembly[0][$i]}: " . length($pacbio1Sequence) . "\n";
			print PACBIONUCMER ">$contigsToName{$orderedAssembly[0][$i]}.$orderedAssembly[0][$i]\n$pacbio1Sequence\n";
		}
		print BACNUCMER ">$bacsThatSpanGaps{$orderedAssembly[0][$i]}\n$yBACFasta{$bacsThatSpanGaps{$orderedAssembly[0][$i]}}\n";
		print "bac length: $bacsThatSpanGaps{$orderedAssembly[0][$i]}: " . length($yBACFasta{$bacsThatSpanGaps{$orderedAssembly[0][$i]}}) . "\n";
		$pacbio2Sequence = $extendedPacBio{$contigsToName{$orderedAssembly[0][$j]}};
		if ($assemblyOrientation{$orderedAssembly[0][$j]} =~ /-/) {
				$pacbio2Sequence = reverse($pacbio2Sequence);
				$pacbio2Sequence =~ tr/ACGTacgt/TGCAtgca/;
		}
		print PACBIONUCMER ">$contigsToName{$orderedAssembly[0][$j]}.$orderedAssembly[0][$j]\n$pacbio2Sequence\n";
		print "connecting pacbio length $contigsToName{$orderedAssembly[0][$j]}: " . length($pacbio2Sequence) . "\n";

		#run nucmer, delta-filter, and quickmerge
		system("nucmer -l 100 --prefix temp nucmer.pacbio.fa nucmer.bac.fa");
		system("delta-filter -r -q -l 10000 -i 99 temp.delta > temp.filtered.delta");
		system("show-coords -lTH temp.filtered.delta > temp.filtered.delta.coords");
		
		#determine if they aligned. If not, keep the pacbio sequence
		open (COORDS, "temp.filtered.delta.coords") || die "Unable to open show-coords output file: $!\n";
		
		my $bacLength;
		my $pacbioLength1;
		my $pacbioLength2;
		my $orientation1;
		my $orientation2;
		my @bac1AlignPositions;
		my @bac2AlignPositions;
		my @pacbio1AlignPositions;
		my @pacbio2AlignPositions;
		
		while (<COORDS>) {
			my @splitLine = split(/\t/, $_);
			chomp(@splitLine);
			print "@splitLine\n";
			my $pacbio = $splitLine[9];
			my $bac = $splitLine[10];
			if ($pacbio =~ /$contigsToName{$orderedAssembly[0][$i]}.$orderedAssembly[0][$i]/) {
				push(@bac1AlignPositions, $splitLine[2]);
				push(@bac1AlignPositions, $splitLine[3]);
				push(@pacbio1AlignPositions, $splitLine[0]);
				push(@pacbio1AlignPositions, $splitLine[1]);
				$bacLength = $splitLine[8];
				$pacbioLength1 = $splitLine[7];
				if(!defined($orientation1) && $splitLine[2] < $splitLine[3]) {
					$orientation1 = 'forward';
				}
				if(!defined($orientation1) && $splitLine[2] > $splitLine[3]) {
					$orientation1 = 'reverse';
				}
			}
			if ($pacbio =~ /$contigsToName{$orderedAssembly[0][$j]}.$orderedAssembly[0][$j]/) {
				push(@bac2AlignPositions, $splitLine[2]);
				push(@bac2AlignPositions, $splitLine[3]);
				push(@pacbio2AlignPositions, $splitLine[0]);
				push(@pacbio2AlignPositions, $splitLine[1]);
				$bacLength = $splitLine[8];
				$pacbioLength2 = $splitLine[7];
				if(!defined($orientation2) && $splitLine[2] < $splitLine[3]) {
					$orientation2 = 'forward';
				}
				if(!defined($orientation2) && $splitLine[2] > $splitLine[3]) {
					$orientation2 = 'reverse';
				}
			}
		}
		my @sortedBAC1AlignPositions = sort {$a <=> $b} @bac1AlignPositions;
		my @sortedBAC2AlignPositions = sort {$a <=> $b} @bac2AlignPositions;
		my @sortedPacbio1AlignPositions = sort {$a <=> $b} @pacbio1AlignPositions;
		my @sortedPacbio2AlignPositions = sort {$a <=> $b} @pacbio2AlignPositions;
		my $pacbioOrientation1 = 'forward';
		my $pacbioOrientation2 = 'forward';
		
		my $mergedFasta;
		
		if ($orientation1 =~ /forward/ && $orientation2 =~ /forward/) {
			print "Span: forward forward\n";
			if ($sortedPacbio1AlignPositions[-1] != $pacbioLength1 && $sortedPacbio2AlignPositions[0] != 1) {
				print "alignment does not reach end of pacbio segments\n";
				if ($sortedPacbio1AlignPositions[-1] != $pacbioLength1) {
					$pacbio1Sequence = substr($pacbio1Sequence, 0, $sortedPacbio1AlignPositions[-1]);
					print "test length 1: " . length($pacbio1Sequence) . "\n";
				}
				if ($sortedPacbio2AlignPositions[0] != 1) {
					$pacbio2Sequence = substr($pacbio2Sequence, $sortedPacbio2AlignPositions[0] - 1);
					print "test length 2: " . length($pacbio2Sequence) . "\n";
				}
			}
			my $bac1AlignLength = $sortedBAC1AlignPositions[-1] - $sortedBAC1AlignPositions[0];
			my $bac2AlignLength = $sortedBAC2AlignPositions[-1] - $sortedBAC2AlignPositions[0];
			my $tempBACSequence = substr($yBACFasta{$bacsThatSpanGaps{$orderedAssembly[0][$i]}}, $bac1AlignLength - 1);
			substr($tempBACSequence, 0 - $bac2AlignLength - 1) = "";

			print "length of bac inserted: " . length($tempBACSequence) . "\n";
			$mergedFasta = $pacbio1Sequence . $tempBACSequence . $pacbio2Sequence;
		}
		elsif ($orientation1 =~ /reverse/ && $orientation2 =~ /reverse/) {
			print "Span: reverse reverse\n";
			if ($sortedPacbio1AlignPositions[-1] != $pacbioLength1 && $sortedPacbio2AlignPositions[-1] != $pacbioLength2) {
				print "alignment does not reach end of pacbio segments\n";
				if ($sortedPacbio1AlignPositions[-1] != $pacbioLength1) {
					$pacbio1Sequence = substr($pacbio1Sequence, 0, $sortedPacbio1AlignPositions[-1]);
					print "test length 1: " . length($pacbio1Sequence) . "\n";
				}
				if ($sortedPacbio2AlignPositions[0] != 1) {
					$pacbio2Sequence = substr($pacbio2Sequence, $sortedPacbio2AlignPositions[0] - 1);
					print "test length 2: " . length($pacbio2Sequence) . "\n";
				}
			}
			my $bac1AlignLength = $sortedBAC1AlignPositions[-1] - $sortedBAC1AlignPositions[0];
			my $bac2AlignLength = $sortedBAC2AlignPositions[-1] - $sortedBAC2AlignPositions[0];
			my $tempBACSequence = substr($yBACFasta{$bacsThatSpanGaps{$orderedAssembly[0][$i]}}, $bac2AlignLength - 1);
			substr($tempBACSequence, 0 - $bac1AlignLength - 1) = "";

			$tempBACSequence = reverse($tempBACSequence);
			$tempBACSequence =~ tr/ACGTacgt/TGCAtgca/;
			print "length of bac inserted: " . length($tempBACSequence) . "\n";
			$mergedFasta = $pacbio1Sequence . $tempBACSequence . $pacbio2Sequence;
		}
		elsif ($orientation1 =~ /forward/ && $orientation2 =~ /reverse/) {
			print "Span: forward reverse\nImproper pacbio contig orientation\n";
			exit(0);
		}
		
		elsif ($orientation1 =~ /reverse/ && $orientation2 =~ /forward/) {
			print "Span: reverse forward\nImproper pacbio contig orientation\n";
			exit(0);
		}
		my $length = length($mergedFasta);
		print "After merged: $length\n";
		
		$joinedContigs{$contigsToName{$orderedAssembly[0][$j]}} = $mergedFasta;
		my @contigs = ($orderedAssembly[0][$j]);
		if (exists($joinedContigsComposition{$orderedAssembly[0][$i]})) {
			unshift(@contigs, @{ $joinedContigsComposition{$orderedAssembly[0][$i]} });
			delete($joinedContigsComposition{$orderedAssembly[0][$i]});
		}
		else {
			unshift(@contigs, $orderedAssembly[0][$i]);
		}
		print "@contigs\n";
		$joinedContigsComposition{$orderedAssembly[0][$j]} = \@contigs;
	}
}


#print out revised assembly
my $finalContig = 'false';
my @finalSequence;
my %printedContigs;
open (OUT, ">separatedContigs.faTEST") || die "Unable to write output file: $!\n";
for (my $i = scalar(@{ $orderedAssembly[0] }) - 1; $i >= 0; $i--) {
	#skip the contigs that should not be in the assembly
	if ($orderedAssembly[0][$i] == 77) {
		last;
	}
	#skip contigs that are already in the assembly
	next if exists($printedContigs{$orderedAssembly[0][$i]});

	if (exists($joinedContigsComposition{$orderedAssembly[0][$i]})) {
		my $tempSequence = $joinedContigs{$contigsToName{$orderedAssembly[0][$i]}};

		unshift(@finalSequence, $tempSequence);
		foreach $contig (@{ $joinedContigsComposition{$orderedAssembly[0][$i]} }) {
			$printedContigs{$contig} = "";
		}
		print OUT ">$orderedAssembly[0][$i]\n$tempSequence\n";
	}
	else {
		#determine orientation
		my $tempSequence = $extendedPacBio{$contigsToName{$orderedAssembly[0][$i]}};
		if ($assemblyOrientation{$orderedAssembly[0][$i]} =~ /-/) {
 			$tempSequence = reverse($tempSequence);
  			$tempSequence =~ tr/ACGTacgt/TGCAtgca/;
 		}
		unshift(@finalSequence, $tempSequence);
		print OUT ">$orderedAssembly[0][$i]\n$tempSequence\n";
		$printedContigs{$orderedAssembly[0][$i]} = "";
	}
}
print "Initial contig count: " . scalar(@{ $orderedAssembly[0] }) . "\n";
print "Final contig count: " . scalar(@finalSequence) . "\n";
open (OUT2, ">mergedFasta_BAC.faTEST") || die "Unable to write output fasta file: $!\n";
print OUT2 ">joinedYChr\n" . join('N' x 500, @finalSequence) . "\n";