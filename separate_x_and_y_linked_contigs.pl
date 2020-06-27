#Purpose: Separate putative X and Y-linked contigs based on nucmer alignment
#Use:  separate_x_and_y_linked_contigs.pl <nucmer coordinates output file> <Canu contigs fasta file>

#use strict;
use warnings;
use Data::Dumper;

open (NUC, "$ARGV[0]") || die "Unable to open nucmer file: $!\n";
open (FASTA, "$ARGV[1]") || die "Unable to open Canu contigs: $!\n";

my %fasta;

#separate out contigs from fasta file
my $header;
while(<FASTA>) {
	if ($_ !~ /^>/) {
		chomp($_);
		$fasta{$header} .= $_;
	}
	else {
	 	chomp($_);
	 	my @subHeader = split(/\s+/, $_);
	 	$header = substr($subHeader[0], 1);
	 	$fasta{$header} = "";
	}
}

my %contigs;
my %contigLength;
my %contigIdentity;
my %maxIdentity;
while (<NUC>) {
	my @splitLine = split(/\s+/, $_);
	my $contig = $splitLine[10];
	my $alignLength = $splitLine[5];
	my $percentIdentity = $splitLine[6] * $alignLength;
	my $chr = $splitLine[9];
	my $length = $splitLine[8];
	if (!exists($contigLength{$contig})) {
		$contigLength{$contig} = $length;
	}
	if (!exists($contigs{$contig}{$chr})) {
		$contigs{$contig}{$chr} = $alignLength;
		$contigIdentity{$contig}{$chr} = $percentIdentity;
	}
	else {
		$contigs{$contig}{$chr} += $alignLength;
		$contigIdentity{$contig}{$chr} += $percentIdentity;
	}
}
my @chrXIXContigList;
my @autoContigList;
my @noAlignContigList;

foreach my $contig (keys(%contigLength)) {
	my $percentAlign = 0.25;
	my $maxChr = 'none';
	foreach my $chr (keys(%{ $contigs{$contig} })) {
		my $proportionAlign = $contigs{$contig}{$chr} / $contigLength{$contig};
		if ($proportionAlign > $percentAlign) {
			$percentAlign = $proportionAlign;
			$maxChr = $chr;
			$maxIdentity{$contig} = $contigIdentity{$contig}{$chr} / $contigs{$contig}{$chr};
		}
	}
	if ($maxChr =~ /chrXIX/) {
		push(@chrXIXContigList, $contig);
	}
	elsif ($maxChr =~ /none/) {
		push(@noAlignContigList, $contig);
	}
	else {
		push(@autoContigList, $contig);
	}
}
#print out contigs. Determine contig alignment with greatest alignment length
my $counter = 1;
my $chrXIXCount = 0;
my $autoCount = 0;
my $noneWithAlignCount = 0;
my $noneWithoutAlignCount = 0;
my $yPossibleCount = 0;
my $autoLength = 0;
my $chrXIXLength = 0;
my $yPossibleLength = 0;
my $noneWithAlignLength = 0;
my $noneWithoutAlignLength = 0;
open (CHRXIX, ">$ARGV[0].chrXIX.txt") || die "Unable to write output file\n";
open (CHRXIXYFASTA, ">$ARGV[0].chrXIX_YPossible.fa") || die "Unable to write output file\n";
open (CHRXIXFASTA, ">$ARGV[0].chrXIX.fa") || die "Unable to write output file\n";
foreach my $contig (@chrXIXContigList) {
	 print CHRXIX "$maxIdentity{$contig}\t$contigLength{$contig}\n";
	if ($maxIdentity{$contig} <= 96) {
		print CHRXIXYFASTA ">$contig\n$fasta{$contig}\n";
		$yPossibleCount++;
		$yPossibleLength += $contigLength{$contig};
	}
	else {
		print CHRXIXFASTA ">$contig\n$fasta{$contig}\n";
		$chrXIXCount++;
		$chrXIXLength  += $contigLength{$contig};
	}
}
open (AUTO, ">$ARGV[0].auto.txt") || die "Unable to write output file\n";
open (AUTOFASTA, ">$ARGV[0].auto.fa") || die "Unable to write output file\n";
foreach my $contig (@autoContigList) {
	print AUTO "$contig\t$maxIdentity{$contig}\t$contigLength{$contig}\n";
 	print AUTOFASTA ">$contig\n$fasta{$contig}\n";
	$autoCount++;
	$autoLength += $contigLength{$contig};
}
open (NONEFASTA, ">$ARGV[0].noneWithAlign.fa") || die "Unable to write output file\n";
open (NONEWITHOUT, ">$ARGV[0].noneWithoutAlign.txt") || die "Unable to write output file\n";
open (NONEWITHOUTFASTA, ">$ARGV[0].noneWithoutAlign.fa") || die "Unable to write output file\n";
foreach my $contig (@noAlignContigList) {
	print NONEFASTA ">$contig\n$fasta{$contig}\n";
	$noneWithAlignCount++;
	$noneWithAlignLength += $contigLength{$contig};
}

#determine remaining contigs that did not align at all to the genome
foreach my $totalContig (keys(%fasta)) {
	if (!exists($contigLength{$totalContig})) {
		print NONEWITHOUT "$totalContig\t" . length($fasta{$totalContig}) . "\n";
		print NONEWITHOUTFASTA ">$totalContig\n$fasta{$totalContig}\n";
		$noneWithoutAlignCount++;
		$noneWithoutAlignLength += length($fasta{$totalContig});
		
	}
}

print "Auto align: $autoCount\tAuto length: $autoLength\n";
print "None with align: $noneWithAlignCount\tNone length: $noneWithAlignLength\n";
print "None without align: $noneWithoutAlignCount\tNone length: $noneWithoutAlignLength\n";
print "chrXIX align: $chrXIXCount\tchrXIX length: $chrXIXLength\n";
print "Y possible: $yPossibleCount\tY possible length: $yPossibleLength\n";
