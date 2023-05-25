#!/usr/bin/env perl

# PODNAME: reformat_bsa_vcf.pl
# ABSTRACT: Reformat a BSA VCF for easier prioritisation

## Author     : i.sealy@qmul.ac.uk
## Maintainer : i.sealy@qmul.ac.uk
## Created    : 2023-05-25

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

# Default options
my ( $high, $moderate );
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Iterate over STDIN
while ( my $line = <> ) {
    next if $line =~ m/\A#/ms;
    chomp $line;
    my (
        $chrom,  $pos,  $id,     $ref, $alt, $qual,
        $filter, $info, $format, $mut, $sib
    ) = split /\t/xms, $line;
    next if $alt =~ m/,/xms;                         # Ignore multiallelic
    next if length $ref != 1 || length $alt != 1;    # Ignore indels
    my $csqs = $info;
    $csqs =~ s/.*CSQ=//xms;
    $csqs =~ s/;.*//xms;
    my @csqs      = split /,/xms, $csqs;
    my $mut_reads = $mut;
    $mut_reads =~ s/.*://xms;
    my ( $mut_ref_reads, $mut_alt_reads ) = split /,/xms, $mut_reads;
    my $sib_reads = $sib;
    $sib_reads =~ s/.*://xms;
    my ( $sib_ref_reads, $sib_alt_reads ) = split /,/xms, $sib_reads;
    next if $mut_ref_reads + $mut_alt_reads == 0;    # No mut reads
    next if $sib_ref_reads + $sib_alt_reads == 0;    # No sib reads
    next
      if $sib_alt_reads / ( $sib_ref_reads + $sib_alt_reads ) >=
      $mut_alt_reads / ( $mut_ref_reads + $mut_alt_reads )
      ;    # Ignore if more ref than alt
    my $sum = $mut_ref_reads + $mut_alt_reads + $sib_ref_reads + $sib_alt_reads;
    my $np1 =
      ( $mut_alt_reads + $sib_alt_reads ) * ( $mut_alt_reads + $mut_ref_reads )
      / $sum;
    my $np2 =
      ( $mut_alt_reads + $sib_alt_reads ) * ( $sib_alt_reads + $sib_ref_reads )
      / $sum;
    my $np3 =
      ( $mut_ref_reads + $sib_ref_reads ) * ( $mut_alt_reads + $mut_ref_reads )
      / $sum;
    my $np4 =
      ( $sib_ref_reads + $mut_ref_reads ) * ( $sib_ref_reads + $sib_alt_reads )
      / $sum;
    my $v1 =
      ( $mut_alt_reads > 0 )
      ? $mut_alt_reads * log( $mut_alt_reads / $np1 )
      : 0;
    my $v2 =
      ( $sib_alt_reads > 0 )
      ? $sib_alt_reads * log( $sib_alt_reads / $np2 )
      : 0;
    my $v3 =
      ( $mut_ref_reads > 0 )
      ? $mut_ref_reads * log( $mut_ref_reads / $np3 )
      : 0;
    my $v4 =
      ( $sib_ref_reads > 0 )
      ? $sib_ref_reads * log( $sib_ref_reads / $np4 )
      : 0;
    my $g = 2 * ( $v1 + $v2 + $v3 + $v4 );

    foreach my $csq (@csqs) {
        next
          if $high
          && $csq !~
m/\b(transcript_ablation|splice_acceptor_variant|splice_donor_variant|stop_gained|frameshift_variant|stop_lost|start_losttranscript_amplification)\b/xms
          && $csq !~ m/missense_variant.*deleterious/xms;
        next
          if $moderate
          && $csq !~
m/\b(inframe_insertion|inframe_deletion|protein_altering_variant)\b/xms
          && $csq !~ m/missense_variant.*tolerated/xms;
        printf "%s:%s\t%s/%s\t%s\t%s\t%s/%s\t%s/%s\t%.1f\n",
          $chrom, $pos, $ref, $alt, $qual, $csq, $mut_ref_reads,
          $mut_alt_reads, $sib_ref_reads, $sib_alt_reads, $g;
    }
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'high'     => \$high,
        'moderate' => \$moderate,
        'debug'    => \$debug,
        'help'     => \$help,
        'man'      => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( $high && $moderate ) {
        pod2usage("Only specify one of --high and --moderate\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

reformat_bsa_vcf.pl

Reformat a BSA VCF for easier prioritisation

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a VCF file produced by munge_vcf_for_bsa.pl on STDIN and
reformats and filters it.

=head1 EXAMPLES

    perl \
        reformat_bsa_vcf.pl --high 
        < input.vcf > output.txt

=head1 USAGE

    munge_vcf_for_bsa.pl
        [--high]
        [--moderate]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--high>

Only include variants that have a high impact consequence, as defined here:
https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html

=item B<--moderate>

Only include variants that have a moderate impact consequence, as defined here:
https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html

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
