#!/gsc/bin/perl

use warnings;
use strict;

use Getopt::Long;
use FileHandle;

=head1 VERSION

    Version 1.01

=head1 SYNOPSIS

    fpfilter is an accessory script to VarScan 2 (http://varscan.sourceforge.net) that filters false positives using information from bam-readcounts (https://github.com/genome/bam-readcount)

=head1 USAGE    fpfilter.pl [varScan file] [bam-readcounts output] OPTIONS

=head3 COMMANDS:

=cut


#        OPTIONS:
#        --output-basename  The basename for output files. Two will be created: basename.pass and basename.fail



our $VERSION = '1.01';

my $usage = qq{USAGE: fpfilter.pl [varScan file] [bam-readcounts output] OPTIONS
        OPTIONS:
        --output-basename   The basename for output files. Two will be created: basename.pass and basename.fail
};

die $usage if(!$ARGV[1]);


## Define filtering parameters ##

my $min_read_pos = 0.10;
my $max_read_pos = 1 - $min_read_pos;
my $min_var_freq = 0.05;
my $min_var_count = 4;

my $min_strandedness = 0.01;
my $max_strandedness = 1 - $min_strandedness;

my $max_mm_qualsum_diff = 50;
my $max_mapqual_diff = 30;
my $max_readlen_diff = 25;
my $min_var_dist_3 = 0.20;
my $max_var_mm_qualsum = 100;


## Parse arguments ##

my $output_basename = "outfile";
my $verbose = 0;
my $ret = parse_arguments();


execute();

################################################################################

=head2	parse_arguments

=cut
################################################################################

sub parse_arguments
{
 #   my $output_alignments, my $output_snps, my $output_indels, my $fasta_file;#, my $quality_file;
#    my $min_identity, my $default_qual_score, my $min_qual_score, my $primer_trim, my $verbose;

    my $result = GetOptions (
                                "output-basename=s"   => \$output_basename,
                                "verbose=s"   => \$verbose,
    );    

}




################################################################################

=head2	execute [ALIGNMENTS] [OPTIONS]

    Run the filter 


=cut


sub execute
{
    my %stats = ();
    $stats{'num_variants'} = $stats{'num_fail_pos'} = $stats{'num_fail_strand'} = $stats{'num_fail_varcount'} = $stats{'num_fail_varfreq'} = $stats{'num_fail_mmqs'} = $stats{'num_fail_var_mmqs'} = $stats{'num_fail_mapqual'} = $stats{'num_fail_readlen'} = $stats{'num_fail_dist3'} = $stats{'num_pass_filter'} = 0;
    
    ## Load the read counts ##
    
    my %readcounts_by_position = ();

    my $rc_input = new FileHandle ($ARGV[1]);
    my $lineCounter = 0;

    while (<$rc_input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;
            my ($chrom, $position) = split(/\t/, $line);
            $readcounts_by_position{"$chrom\t$position"} = $line;
    }
    
    close($rc_input);


    ## Open the output files ##
    
    open(PASS, ">$output_basename.pass") or die "Can't open output file: $!\n";
    open(FAIL, ">$output_basename.fail") or die "Can't open output file: $!\n";

    ## Parse the input file ##


    my $input = new FileHandle ($ARGV[0]);
    $lineCounter = 0;

    while (<$input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;
            my ($chrom, $position, $ref, $var) = split(/\t/, $line);
            $ref = uc($ref);
            $var = uc($var);

            if(!($var =~ /[ACGT]/)) {
                $var = iupac_to_base($ref, $var);
            }
            
            $stats{'num_variants'}++;
            
            if($readcounts_by_position{"$chrom\t$position"})
            {
                my $readcounts = $readcounts_by_position{"$chrom\t$position"};
                my $ref_result = read_counts_by_allele($readcounts, $ref);
                my $var_result = read_counts_by_allele($readcounts, $var);
                
                if($ref_result && $var_result)
                {
                        ## Parse out the bam-readcounts details for each allele. The fields should be: ##
                        #num_reads : avg_mapqual : avg_basequal : avg_semq : reads_plus : reads_minus : avg_clip_read_pos : avg_mmqs : reads_q2 : avg_dist_to_q2 : avgRLclipped : avg_eff_3'_dist
                        my ($ref_count, $ref_map_qual, $ref_base_qual, $ref_semq, $ref_plus, $ref_minus, $ref_pos, $ref_subs, $ref_mmqs, $ref_q2_reads, $ref_q2_dist, $ref_avg_rl, $ref_dist_3) = split(/\t/, $ref_result);
                        my ($var_count, $var_map_qual, $var_base_qual, $var_semq, $var_plus, $var_minus, $var_pos, $var_subs, $var_mmqs, $var_q2_reads, $var_q2_dist, $var_avg_rl, $var_dist_3) = split(/\t/, $var_result);

                        my $ref_strandedness = my $var_strandedness = 0.50;
                        $ref_dist_3 = 0.5 if(!$ref_dist_3);

                        ## Use conservative defaults if we can't get mismatch quality sums ##
                        $ref_mmqs = 50 if(!$ref_mmqs);
                        $var_mmqs = 0 if(!$var_mmqs);
                        my $mismatch_qualsum_diff = $var_mmqs - $ref_mmqs;

                        ## Determine map qual diff ##

                        my $mapqual_diff = $ref_map_qual - $var_map_qual;


                        ## Determine difference in average supporting read length ##

                        my $readlen_diff = $ref_avg_rl - $var_avg_rl;


                        ## Determine ref strandedness ##

                        if(($ref_plus + $ref_minus) > 0) {
                            $ref_strandedness = $ref_plus / ($ref_plus + $ref_minus);
                            $ref_strandedness = sprintf("%.2f", $ref_strandedness);
                        }

                        ## Determine var strandedness ##

                        if(($var_plus + $var_minus) > 0) {
                            $var_strandedness = $var_plus / ($var_plus + $var_minus);
                            $var_strandedness = sprintf("%.2f", $var_strandedness);
                        }
                        
                        if($var_count && ($var_plus + $var_minus))
                        {
                            ## We must obtain variant read counts to proceed ##

                            my $var_freq = $var_count / ($ref_count + $var_count);

                            ## FAILURE 1: READ POSITION ##
                            if(($var_pos < $min_read_pos)) { # || $var_pos > $max_read_pos)) {
                                print FAIL "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tReadPos<$min_read_pos\n";
                                print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tReadPos<$min_read_pos\n"if ($verbose);
                                $stats{'num_fail_pos'}++;
                            }

                            ## FAILURE 2: Variant is strand-specific but reference is NOT strand-specific ##
                            elsif(($var_strandedness < $min_strandedness || $var_strandedness > $max_strandedness) && ($ref_strandedness >= $min_strandedness && $ref_strandedness <= $max_strandedness)) {
                                ## Print failure to output file if desired ##
                                print FAIL "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tStrandedness: Ref=$ref_strandedness Var=$var_strandedness\n";
                                print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tStrandedness: Ref=$ref_strandedness Var=$var_strandedness\n"if ($verbose);
                                $stats{'num_fail_strand'}++;
                            }

                            ## FAILURE : Variant allele count does not meet minimum ##
                            elsif($var_count < $min_var_count) {
                                print FAIL "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarCount:$var_count\n";
                                print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarCount:$var_count\n" if ($verbose);
                                $stats{'num_fail_varcount'}++;
                            }

                            ## FAILURE : Variant allele frequency does not meet minimum ##
                            elsif($var_freq < $min_var_freq) {
                                print FAIL "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarFreq:$var_freq\n";
                                print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarFreq:$var_freq\n" if ($verbose);
                                $stats{'num_fail_varfreq'}++;
                            }

                            ## FAILURE 3: Paralog filter for sites where variant allele mismatch-quality-sum is significantly higher than reference allele mmqs
                            elsif($mismatch_qualsum_diff> $max_mm_qualsum_diff) {
                                ## Print failure to output file if desired ##
                                print FAIL "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tMismatchQualsum:$var_mmqs-$ref_mmqs=$mismatch_qualsum_diff\n";
                                print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tMismatchQualsum:$var_mmqs-$ref_mmqs=$mismatch_qualsum_diff" if ($verbose);
                                $stats{'num_fail_mmqs'}++;
                            }

                            ## FAILURE 4: Mapping quality difference exceeds allowable maximum ##
                            elsif($mapqual_diff > $max_mapqual_diff) {
                                print FAIL "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tMapQual:$ref_map_qual-$var_map_qual=$mapqual_diff\n";
                                print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tMapQual:$ref_map_qual-$var_map_qual=$mapqual_diff" if ($verbose);
                                $stats{'num_fail_mapqual'}++;
                            }

                            ## FAILURE 5: Read length difference exceeds allowable maximum ##
                            elsif($readlen_diff > $max_readlen_diff) {
                                print FAIL "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tReadLen:$ref_avg_rl-$var_avg_rl=$readlen_diff\n";
                                print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tReadLen:$ref_avg_rl-$var_avg_rl=$readlen_diff" if ($verbose);
                                $stats{'num_fail_readlen'}++;
                            }

                            ## FAILURE 5: Read length difference exceeds allowable maximum ##
                            elsif($var_dist_3 < $min_var_dist_3) {
                                print FAIL "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarDist3:$var_dist_3\n";
                                print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarDist3:$var_dist_3\n" if ($verbose);
                                $stats{'num_fail_dist3'}++;
                            }

                            elsif($max_var_mm_qualsum && $var_mmqs > $max_var_mm_qualsum) {
                                print FAIL "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarMMQS: $var_mmqs > $max_var_mm_qualsum\n";
                                print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarMMQS: $var_mmqs > $max_var_mm_qualsum\n" if ($verbose);
                                $stats{'num_fail_var_mmqs'}++;
                            }

                            ## SUCCESS: Pass Filter ##
                            else {
                                $stats{'num_pass_filter'}++;
                                ## Print output, and append strandedness information ##
                                print PASS "$line\n";
                                print "$line\tPASS\n" if($verbose);
                            }

                        }
                }
                else
                {
                    $stats{'num_no_readcounts'}++;
                    print FAIL "$line\tno_readcounts\n";
                }
            }
            else
            {
                $stats{'num_no_readcounts'}++;
                print FAIL "$line\tno_readcounts\n";
            }
   

    }
    
    close($input);
    
    close(PASS);
    close(FAIL);

    ## Print filtering stats ##

    print $stats{'num_variants'} . " variants\n";
    print $stats{'num_no_readcounts'} . " failed to get readcounts for variant allele\n";
    print $stats{'num_fail_pos'} . " had read position < $min_read_pos\n";
    print $stats{'num_fail_strand'} . " had strandedness < $min_strandedness\n";
    print $stats{'num_fail_varcount'} . " had var_count < $min_var_count\n";
    print $stats{'num_fail_varfreq'} . " had var_freq < $min_var_freq\n";

    print $stats{'num_fail_mmqs'} . " had mismatch qualsum difference > $max_mm_qualsum_diff\n";
    print $stats{'num_fail_var_mmqs'} . " had variant MMQS > $max_var_mm_qualsum\n" if($stats{'num_fail_var_mmqs'});
    print $stats{'num_fail_mapqual'} . " had mapping quality difference > $max_mapqual_diff\n";
    print $stats{'num_fail_readlen'} . " had read length difference > $max_readlen_diff\n";
    print $stats{'num_fail_dist3'} . " had var_distance_to_3' < $min_var_dist_3\n";

    print $stats{'num_pass_filter'} . " passed the strand filter\n";

    return(0);
}


################################################################################

=head3	iupac_to_base

    Convert IUPAC ambiguity codes to variant bases


=cut


sub iupac_to_base {
    (my $allele1, my $allele2) = @_;

    return($allele2) if($allele2 eq "A" || $allele2 eq "C" || $allele2 eq "G" || $allele2 eq "T");

    if($allele2 eq "M") {
        return("C") if($allele1 eq "A");
        return("A") if($allele1 eq "C");
        return("A");    ## Default for triallelic variant
    } elsif($allele2 eq "R") {
        return("G") if($allele1 eq "A");
        return("A") if($allele1 eq "G");
        return("A");     ## Default for triallelic variant
    } elsif($allele2 eq "W") {
        return("T") if($allele1 eq "A");
        return("A") if($allele1 eq "T");
        return("A");    ## Default for triallelic variant
    } elsif($allele2 eq "S") {
        return("C") if($allele1 eq "G");
        return("G") if($allele1 eq "C");
        return("C");    ## Default for triallelic variant
    } elsif($allele2 eq "Y") {
        return("C") if($allele1 eq "T");
        return("T") if($allele1 eq "C");
        return("C");    ## Default for triallelic variant
    } elsif($allele2 eq "K") {
        return("G") if($allele1 eq "T");
        return("T") if($allele1 eq "G");
        return("G");    ## Default for triallelic variant
    }

    return($allele2);
}


################################################################################

=head3	read_counts_by_allele

    Retrieve relevant read counts for a certain allele 


=cut

sub read_counts_by_allele {
    (my $line, my $allele) = @_;

    my @lineContents = split(/\t/, $line);
    my $numContents = @lineContents;

    for(my $colCounter = 5; $colCounter < $numContents; $colCounter++) {
        my $this_allele = $lineContents[$colCounter];
        my @alleleContents = split(/\:/, $this_allele);
        if($alleleContents[0] eq $allele) {
            my $numAlleleContents = @alleleContents;

            return("") if($numAlleleContents < 8);

            my $return_string = "";
            my $return_sum = 0;
            for(my $printCounter = 1; $printCounter < $numAlleleContents; $printCounter++) {
                $return_sum += $alleleContents[$printCounter];
                $return_string .= "\t" if($return_string);
                $return_string .= $alleleContents[$printCounter];
            }

            return($return_string);

        }
    }

    return("");
}


=head1 AUTHOR

    Daniel C. Koboldt, << <dkoboldt at genome.wustl.edu> >>
    The Genome Institute at Washington University School of Medicine
    St. Louis, Missouri, USA

=head1 COPYRIGHT

    Copyright 2009-2012 Daniel C. Koboldt and Washington University
    All rights reserved.

=head1 LICENSE

    This program is free for non-commercial use.

=cut
