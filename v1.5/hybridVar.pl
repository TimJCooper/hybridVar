#!/usr/bin/env perl
#Version: 1.5

####################################################################################################################
# Author(s): T.J.Cooper
# Updated: 7/9/2016
# Extracts high-quality, filtered variants from VCF hybrid-spore data and constructs a variant reference genome
# Requires non-core module (Bio::SeqIO)
# GATK Example: java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R S288c.fa -I sample_sorted.bam -o sample_out.vcf
####################################################################################################################

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename qw(basename);
use List::Util qw(all);
my @files = glob("*.vcf");
my $chk   = scalar(@files);
print "\nFailed to detect any .VCF files within the current directory.\n\n"
  if $chk == 0;
exit if $chk == 0;    #Stop script if no .vcf files are found
my ( $fasta, $tRDfilt, $vRDfilt, $freqlow, $freqhigh );
my $scriptname = basename($0);    #Obtain script-name
my $usage =
"Usage: $scriptname -r <ReferenceFASTA> -lf <CallFreqLowerLim> -uf <CallFreqUpperLim> -trd <MinReadDepth> -vd <MinVarDepth>"; #Error/usage message
GetOptions(
    'r=s'   => \$fasta,   #Command-line arguments
    'lf=f'  => \$freqlow,
    'uf=f'  => \$freqhigh,
    'trd=i' => \$tRDfilt,
    'vd=f'  => \$vRDfilt
) or die("\n$usage\n");
die(
"\nError: Arguments or -flags are missing and/or incorrectly specified.\n\n$usage\n\n"
) unless all { defined } $fasta, $freqlow, $freqhigh, $tRDfilt, $vRDfilt;
print "\n\n$chk Samples Detected";
print "\n------------------------------------------";
print "\nFiltering and Merging Variants...\n";
print "------------------------------------------\n";
my ( @calls, @ID, @refsplit );
my ( %freq, %RD, %AD, %offset, %sequences, %refdupl, %dups );
my (
    $uID,     $varcount, $discard, $snp, $indel,
    $overlap, $rkey,     $dupkey,  $splituID
);
my $outfile  = "VariantStats.txt";
my $outfile2 = "LowQualVariants.txt";
open my $OUT,  '>', "CallStats.txt"       or die "$!";
open my $OUT2, '>', "LowQualVariants.txt" or die "$!";
print $OUT "uID\tChr\tPos\tRef\tVar\tCallFreq\ttRD\tvRD/tRD\n";
print $OUT2 "Chr\tPos\tRef\tVar\tCallFreq\ttRD\tvRD/tRD\n";

for my $file (@files) {    #For-each input file
    open my $IN, '<', $file or die "$!";
    while (<$IN>) {
        next if /^\s*#/;
        chomp $_;
        my @F       = split( "\t",      $_ );
        my @varinfo = split( m[[:,/]+], $F[9] );    #Split genotype information
        next if $varinfo[4] == 0;                   #Skip false-positives
        my $check = index( $F[4], ',' );
        if ( $check == '-1' ) {    #For each unique, mono-allelic variant
            $freq{ $F[0] }{ $F[1] }{ $F[3] }
              { $F[4] }++;         #Calculate call-frequencies
            $RD{ $F[0] }{ $F[1] }{ $F[3] }{ $F[4] } +=
              $varinfo[4];         #Cumulative total of total read-depth (tRD)
            $AD{ $F[0] }{ $F[1] }{ $F[3] }{ $F[4] } +=
              $varinfo[3];         #Cumulative total of allelic read-depth (vRD)
        }
    }
}
my ( $i, $k ) = -1;

sub overlap {                      #Subroutine to identify overlapping variants
    my ( $id, $type, $value ) = @_;
    $refsplit[ ++$i ] = [ $id, $type, $value ];
    $rkey = "$type:$value";
    $refdupl{$rkey} = [] if !exists $refdupl{$rkey};
    push @{ $refdupl{$rkey} }, $i;
    return;
}
print "Identifying overlaps...\n";
print "------------------------------------------\n";
foreach my $chrnum ( sort keys %freq ) {  #For each unique, mono-allelic variant
    foreach my $pos ( sort { $a <=> $b } keys %{ $freq{$chrnum} } ) {
        foreach my $ref ( keys %{ $freq{$chrnum}{$pos} } ) {
            foreach my $var ( keys %{ $freq{$chrnum}{$pos}{$ref} } ) {
                $varcount++;              #Total no. unique variant-count
                my $callfreq  = ( $freq{$chrnum}{$pos}{$ref}{$var} ) / $chk;
                my $readdepth = ( $RD{$chrnum}{$pos}{$ref}{$var} );
                my $vardepth  = ( $AD{$chrnum}{$pos}{$ref}{$var} ) / $readdepth;
                if (   $callfreq > $freqlow
                    && $callfreq < $freqhigh
                    && $readdepth > $tRDfilt
                    && $vardepth > $vRDfilt )
                { #Filter variants using user-specified call-frequency, tRD and vRD thresholds
                    $uID++;
                    printf( $OUT "%d\t%s\t%d\t%s\t%s\t%.3f\t%d\t%.3f\n",
                        $uID, $chrnum, $pos, $ref, $var, $callfreq, $readdepth,
                        $vardepth );
                    if (   length($ref) == length($var)
                        || length($ref) < length($var) )
                    {    #For SNPs or deletions (relative to reference)
                        overlap( $uID, $chrnum, $pos );
                    }
                    elsif ( length($ref) > length($var) )
                    {    #For insertions (relative to reference)
                        my $del = length($ref);
                        overlap( $uID, $chrnum, $pos );
                        foreach my $delsplit ( 1 .. length($ref) - 1 )
                        { #For each additional inserted base (within the reference)
                            overlap( $uID, $chrnum, $pos + $delsplit );
                        }
                    }
                }
                else {
                    $discard++;    #Total no. dicarded variants
                    printf( $OUT2 "%s\t%d\t%s\t%s\t%.3f\t%d\t%.3f\n",
                        $chrnum, $pos, $ref, $var, $callfreq, $readdepth,
                        $vardepth );
                }
            }
        }
    }
}
foreach my $entries (@refsplit) {
    my $dupkey = "$entries->[1]:$entries->[2]";
    if ( @{ $refdupl{$dupkey} } > 1 )
    {    #For any non-unique chr-pos combinations (overlaps)
        $overlap++;    #Total no. overlapping variants
        $dups{ @$entries[0] } = {};    #Store uID of all overlapping variants
    }
}
my $seqio =
  Bio::SeqIO->new( -file => $fasta );    #Read and store .FASTA chromosomes
while ( my $seqobj = $seqio->next_seq ) {
    my $id  = $seqobj->display_id;
    my $seq = $seqobj->seq;
    $sequences{$id} = $seq;
}
close $OUT;
close $OUT2;
print "Constructing variant reference...\n";
print "------------------------------------------\n";
open my $IN2, '<', "CallStats.txt" or die "$!";
<$IN2> for ( 1 .. 1 );                   #Skip headline
open my $OUT3, '>', "VariantTable.txt"         or die "$!";
open my $OUT4, '>', "VariantRef.fa"            or die "$!";
open my $OUT5, '>', "VariantRefChromSizes.txt" or die "$!";
print $OUT3 "uID\tchrom\tpos_c\tpos_k\tseq_c\tseq_sk\ttype_c\ttype_k\n";

while (<$IN2>) {
    chomp $_;
    my @F2 = split( "\t", $_ );          #Split each tab-delimited field
    next if exists( $dups{ $F2[0] } );   #Skip overlapping variants
    if ( defined $offset{ $F2[1] } )
    {    #Offset counters for each chromosome (INDEL-dependent position shifts)
    }
    else {
        $offset{ $F2[1] } = 0;
    }
    $splituID++;
    if ( length( $F2[3] ) == length( $F2[4] ) ) {    #For SNPs
        substr( $sequences{ $F2[1] }, ( $F2[2] - 1 + $offset{ $F2[1] } ), 1 ) =
          $F2[4];                                    #Ref->Var SNP substitution
        $snp++;                                      #Total no. SNPs
        printf( $OUT3 "%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n",
            $splituID, $F2[1], $F2[2], $F2[2] + $offset{ $F2[1] },
            $F2[3], $F2[4], "s", "s"
        );
    }
    elsif ( length( $F2[3] ) < length( $F2[4] ) )
    {    #For deletions (relative to reference)
        substr( $sequences{ $F2[1] }, ( $F2[2] - 1 + $offset{ $F2[1] } ), 1 ) =
          $F2[4];    #Insertion of additional variant bases
        $indel++;    #Total no. INDELs
        printf( $OUT3 "%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n",
            $splituID, $F2[1], $F2[2], $F2[2] + $offset{ $F2[1] },
            $F2[3], substr( $F2[4], 0, 1 ),
            "d", "i"
        );
        foreach my $inssplit ( 1 .. length( $F2[4] ) - 1 )
        {            #Base-by-base split of insertion
            printf( $OUT3 "%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n",
                $splituID, $F2[1], "-", $F2[2] + $offset{ $F2[1] } + $inssplit,
                "-", substr( $F2[4], $inssplit, 1 ),
                "d", "i"
            );
        }
        $offset{ $F2[1] } +=
          length( $F2[4] ) - length( $F2[3] );    #Calculate position offset
    }
    elsif ( length( $F2[3] ) > length( $F2[4] ) )
    {    #For insertions (within the reference)
        my $del = length( $F2[3] );
        substr( $sequences{ $F2[1] }, ( $F2[2] - 1 + $offset{ $F2[1] } ), $del )
          = $F2[4];    #Deletion of inserted bases
        $indel++;      #Total no. INDELs
        printf( $OUT3 "%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n",
            $splituID, $F2[1], $F2[2],
            $F2[2] + $offset{ $F2[1] },
            substr( $F2[3], 0, 1 ),
            $F2[4], "i", "d"
        );
        foreach my $delsplit ( 1 .. length( $F2[3] ) - 1 ) {
            printf( $OUT3 "%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n",
                $splituID, $F2[1], $F2[2] + $delsplit,
                "-", substr( $F2[3], $delsplit, 1 ),
                "-", "i", "d"
            );
        }
        $offset{ $F2[1] } -= $del - 1;    #Calculate position offset
    }
}
for my $chr ( sort keys %sequences ) {    #Construct variant .FASTA file
    print $OUT4 ">$chr\n$sequences{$chr}\n";
    print $OUT5 "$chr\t", length( $sequences{$chr} ), "\n";
}
my $run_time = time() - $^T;
print "Total Variants: $varcount\n";
print "Failed: $discard\n";
print "Overlapping: $overlap\n";
print "Passed: ", $varcount - $discard, " (SNPs: $snp, INDELs: $indel)\n";
print "------------------------------------------\n";
print "Run Completed\n";
print "Processing Runtime: $run_time Seconds\n";
print "------------------------------------------\n\n\n";
