#!/usr/bin/env perl

# PODNAME: munge_vcf_for_bsa.pl
# ABSTRACT: Munge a VCF with separate samples into two bulks for BSA

## Author     : i.sealy@qmul.ac.uk
## Maintainer : i.sealy@qmul.ac.uk
## Created    : 2023-04-12

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

# Default options
my @mut_names;
my @sib_names;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Iterate over STDIN and deal with header
my ( @mut_idx, @sib_idx );
while ( my $line = <> ) {
    chomp $line;
    if ( $line =~ m/\A \#\#/xms ) {
        printf "%s\n", $line;
    }
    elsif ( $line =~ m/\A \#CHROM/xms ) {
        my (
            $chrom, $pos,    $id,   $ref,    $alt,
            $qual,  $filter, $info, $format, @samples
        ) = split /\t/xms, $line;
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tmut\tsib\n", $chrom, $pos,
          $id, $ref, $alt, $qual, $filter, $info, $format;
        my %is_mut = map { $_ => 1 } @mut_names;
        my %is_sib = map { $_ => 1 } @sib_names;
        my $idx    = 0;
        foreach my $sample (@samples) {
            if ( exists $is_mut{$sample} ) {
                push @mut_idx, $idx;
            }
            if ( exists $is_sib{$sample} ) {
                push @sib_idx, $idx;
            }
            $idx++;
        }
        last;
    }
}

# Iterate over STDIN and deal with rest of VCF
while ( my $line = <> ) {
    chomp $line;
    my (
        $chrom, $pos,    $id,   $ref,    $alt,
        $qual,  $filter, $info, $format, @samples
    ) = split /\t/xms, $line;
    next if $alt =~ m/,/xms;                         # Ignore multiallelic
    next if length $ref != 1 || length $alt != 1;    # Ignore indels
    my ( $gt_idx, $ad_idx );
    my $idx = 0;
    foreach my $field ( split /:/xms, $format ) {
        if ( $field eq 'GT' ) {
            $gt_idx = $idx;
        }
        if ( $field eq 'AD' ) {
            $ad_idx = $idx;
        }
        $idx++;
    }
    my ( %mut_allele, %sib_allele );
    my $mut_ad_0 = 0;
    my $mut_ad_1 = 0;
    my $sib_ad_0 = 0;
    my $sib_ad_1 = 0;
    foreach my $idx (@mut_idx) {
        my $sample = $samples[$idx];
        my @fields = split /:/xms, $sample;
        my $gt     = $fields[$gt_idx];
        next if $gt =~ m/\A [.]/xms;
        foreach my $allele ( split /[\/|]/xms, $gt ) {
            $mut_allele{$allele}++;
        }
        my $ad = $fields[$ad_idx];
        my ( $ad_0, $ad_1 ) = split /,/xms, $ad;
        $mut_ad_0 += $ad_0;
        $mut_ad_1 += $ad_1;
    }
    next if scalar keys %mut_allele == 0;    # No mut genotype
    foreach my $idx (@sib_idx) {
        my $sample = $samples[$idx];
        my @fields = split /:/xms, $sample;
        my $gt     = $fields[$gt_idx];
        next if $gt =~ m/\A [.]/xms;
        foreach my $allele ( split /[\/|]/xms, $gt ) {
            $sib_allele{$allele}++;
        }
        my $ad = $fields[$ad_idx];
        my ( $ad_0, $ad_1 ) = split /,/xms, $ad;
        $sib_ad_0 += $ad_0;
        $sib_ad_1 += $ad_1;
    }
    next if scalar keys %sib_allele == 0;    # No sib genotype
    my $mut_gt = '0/1';
    if ( scalar keys %mut_allele == 1 ) {
        my ($allele) = keys %mut_allele;
        $mut_gt = $allele . q{/} . $allele;
    }
    my $sib_gt = '0/1';
    if ( scalar keys %sib_allele == 1 ) {
        my ($allele) = keys %sib_allele;
        $sib_gt = $allele . q{/} . $allele;
    }
    next if $mut_gt eq '0/0' && $sib_gt eq '0/0';    # No ALT alleles
    next if $mut_gt eq '1/1' && $sib_gt eq '1/1';    # No REF alleles
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tGT:AD\t%s:%d,%d\t%s:%d,%d\n",
      $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $mut_gt, $mut_ad_0,
      $mut_ad_1, $sib_gt, $sib_ad_0, $sib_ad_1;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'mut=s@{1,}' => \@mut_names,
        'sib=s@{1,}' => \@sib_names,
        'debug'      => \$debug,
        'help'       => \$help,
        'man'        => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !@mut_names ) {
        pod2usage("--mut must be specified\n");
    }
    if ( !@sib_names ) {
        pod2usage("--sib must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

munge_vcf_for_bsa.pl

Munge a VCF with separate samples into two bulks for BSA

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a VCF file on STDIN and merges multiple samples into "mut"
and "sib" bulks ready for bulk segregant analysis with
https://github.com/xiekunwhy/bsa

Genotypes and allele depths are merged and various filters are applied
(multiallelic positions, indels, positions that aren't informative for BSA,
 etc...).

=head1 EXAMPLES

    perl \
        munge_vcf_for_bsa.pl \
        --mut M1 M2 M3 --sib S1 S2 S3 \
        < input.vcf > output.vcf

=head1 USAGE

    munge_vcf_for_bsa.pl
        [--mut sample...]
        [--sib sample...]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--mut SAMPLE...>

The "mut" samples to be combined into one sample.

=item B<--sib SAMPLE...>

The "sib" samples to be combined into one sample.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=head1 DEPENDENCIES

None

=head1 AUTHOR

=over 4

=item *

Ian Sealy <i.sealy@qmul.ac.uk>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2023 by Ian Sealy.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
