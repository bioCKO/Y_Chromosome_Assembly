#print "***Purpose: Compute mean pairwise identity between centromere repeats\n";
#print "***Use:  centromere_pairwise_identity.pl <blast PSL file>\n";

#use strict;
use warnings;
use Data::Dumper;

sub median {
    my (@data) = sort { $a <=> $b } @_;
    if ( scalar(@data) % 2 ) {
        return ( $data[ @data / 2 ] );
    } else {
        my ( $upper, $lower );
        $lower = $data[ @data / 2 ];
        $upper = $data[ @data / 2 - 1 ];
        return ( mean( $lower, $upper ) );
    }
}
sub mean {
    my (@data) = @_;
    my $sum;
    foreach (@data) {
        $sum += $_;
    }
    return ( $sum / @data );
}

#format of PSL needs to be -outfmt '6 qaccver saccver pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore'
open (PSL, $ARGV[0]) || die "Unable to open PSL file: $!\n";
open (OUT, ">centromerePairwiseComparison.txt") || die "Unable to write output file: $!\n";

my %positions;
my %psl;
my @pairwiseIdentity;

while (<PSL>) {
	my @splitLine = split(/\t/, $_);
	my @start1 = split(/\./, $splitLine[0]);
	my @start2 = split(/\./, $splitLine[1]);
	my $percentIdentity = $splitLine[2];
	next if ($splitLine[3] < 177 || $splitLine[3] > 197);
	if (!exists($psl{$start1[0]}{$start2[0]})) {
		$psl{$start1[0]}{$start2[0]} = $percentIdentity;
	}
	else {
		print "DUPLICATE\n";
	}
	$positions{$start1[0]} = "";
	$positions{$start2[0]} = "";
}

my @keys = keys(%positions);
my @sortedKeys = sort {$a <=> $b} @keys;
print OUT "\t" . join("\t", @sortedKeys) . "\n";

for (my $i = 0; $i < @sortedKeys; $i++) {
	print OUT "$sortedKeys[$i]\t";
	for (my $j = 0; $j < @sortedKeys - 1; $j++) {
		if (exists($psl{$sortedKeys[$i]}{$sortedKeys[$j]})) {
			print OUT "$psl{$sortedKeys[$i]}{$sortedKeys[$j]}\t";
			if ($sortedKeys[$i] != $sortedKeys[$j]) {
				push(@pairwiseIdentity, $psl{$sortedKeys[$i]}{$sortedKeys[$j]});
			}
		}
		else {
			print OUT "NA\t";
		}
	}
	if (exists($psl{$sortedKeys[$i]}{$sortedKeys[-1]})) {
		print OUT "$psl{$sortedKeys[$i]}{$sortedKeys[-1]}\n";
		if ($sortedKeys[$i] != $sortedKeys[-1]) {
			push(@pairwiseIdentity, $psl{$sortedKeys[$i]}{$sortedKeys[-1]});
		}
	}
	else {
		print OUT "NA\n";
	}
}

print "Mean pairwise identity: " . mean(@pairwiseIdentity) . "\n";