#!/usr/bin/env perl

# PODNAME: dump_ensembl_cdna.pl
# ABSTRACT: Dump Ensembl cDNAs in FASTA format

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2014-06-04

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Bio::EnsEMBL::Registry;
use Bio::Seq;
use Bio::SeqIO;

# Default options
my $three_prime_extension = 0;
my $species               = 'Danio rerio';
my $ensembl_dbhost        = 'ensembldb.ensembl.org';
my $ensembl_dbport;
my $ensembl_dbuser = 'anonymous';
my $ensembl_dbpass;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Connect to Ensembl database
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => $ensembl_dbhost,
    -port => $ensembl_dbport,
    -user => $ensembl_dbuser,
    -pass => $ensembl_dbpass,
);

# Get genebuild version
my $genebuild_version = 'e' . Bio::EnsEMBL::ApiVersion::software_version();
warn 'Genebuild version: ', $genebuild_version, "\n" if $debug;

# Get Ensembl adaptors
my $ga = Bio::EnsEMBL::Registry->get_adaptor( $species, 'core', 'Gene' );
my $sa = Bio::EnsEMBL::Registry->get_adaptor( $species, 'core', 'Slice' );

# Ensure database connection isn't lost; Ensembl 64+ can do this more elegantly
## no critic (ProhibitMagicNumbers)
if ( Bio::EnsEMBL::ApiVersion::software_version() < 64 ) {
## use critic
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
}
else {
    Bio::EnsEMBL::Registry->set_reconnect_when_lost();
}

# Dump FASTA to STDOUT
my $fasta_out = Bio::SeqIO->new(
    -fh     => \*STDOUT,
    -format => 'fasta',
);

# Get all genes and their transcripts
my $genes = $ga->fetch_all();
foreach my $gene ( @{$genes} ) {
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript ( @{$transcripts} ) {
        my $seq = Bio::Seq->new(
            -id => ( join q{_}, $transcript->stable_id, $gene->stable_id ),
            -alphabet => 'dna',
            -seq      => $transcript->seq->seq,
        );

        # Optionally extend 3' end
        if ($three_prime_extension) {
            $seq = get_three_prime_extension_seq( $transcript, $seq,
                $three_prime_extension );
        }

        $fasta_out->write_seq($seq);
    }
    $gene->flush_Transcripts();    # Save memory
}

# Get 3' extension sequence
sub get_three_prime_extension_seq {
    my ( $transcript, $seq, $three_prime_extension_length ) = @_;

    # Get slice corresponding to 3' extension
    my ( $extension_start, $extension_end ) =
      get_three_prime_extension_coordinates( $transcript,
        $three_prime_extension_length );
    my $slice = $sa->fetch_by_region( 'toplevel', $transcript->seq_region_name,
        $extension_start, $extension_end, $transcript->strand );

    # Get overlapping genes (because extension will be reduced if overlaps)
    my $overlapping_genes = $slice->get_all_Genes();

    # Remove current gene
    @{$overlapping_genes} =
      grep { $_->stable_id ne $transcript->get_Gene->stable_id }
      @{$overlapping_genes};

    # Remove genes on opposite strand
    @{$overlapping_genes} =
      grep { $_->seq_region_strand eq $transcript->strand }
      @{$overlapping_genes};

    # Work out how much (if any) to reduce extension by
    my $min_distance = $three_prime_extension_length;
    foreach my $overlapping_gene ( @{$overlapping_genes} ) {
        my $distance;
        if ( $transcript->strand > 0 ) {

            # Forward strand
            $distance =
              $overlapping_gene->seq_region_start - $transcript->seq_region_end;
        }
        else {
            # Reverse strand
            $distance =
              $transcript->seq_region_start - $overlapping_gene->seq_region_end;
        }
        if ( $distance < $min_distance ) {
            $min_distance = $distance;
        }
    }

    # Don't extend if extension completely overlaps another gene
    return $seq if $min_distance < 1;

    # Get new slice corresponding to reduced 3' extension
    if ( $min_distance < $three_prime_extension_length ) {
        ( $extension_start, $extension_end ) =
          get_three_prime_extension_coordinates( $transcript, $min_distance );
        $slice = $sa->fetch_by_region( 'toplevel', $transcript->seq_region_name,
            $extension_start, $extension_end, $transcript->strand );
    }

    # Update sequence
    $seq->seq( $seq->seq . $slice->seq );

    return $seq;
}

# Get coordinates of 3' extension
sub get_three_prime_extension_coordinates {
    my ( $transcript, $three_prime_extension_length ) = @_;

    my ( $extension_start, $extension_end );

    if ( $transcript->strand > 0 ) {

        # Forward strand
        $extension_start = $transcript->seq_region_end + 1;
        $extension_end =
          $transcript->seq_region_end + $three_prime_extension_length;
    }
    else {
        # Reverse strand
        $extension_start =
          $transcript->seq_region_start - $three_prime_extension_length;
        $extension_end = $transcript->seq_region_start - 1;
    }

    return $extension_start, $extension_end;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'three_prime_extension=i' => \$three_prime_extension,
        'species=s'               => \$species,
        'ensembl_dbhost=s'        => \$ensembl_dbhost,
        'ensembl_dbport=i'        => \$ensembl_dbport,
        'ensembl_dbuser=s'        => \$ensembl_dbuser,
        'ensembl_dbpass=s'        => \$ensembl_dbpass,
        'debug'                   => \$debug,
        'help'                    => \$help,
        'man'                     => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

dump_ensembl_cdna.pl

Dump Ensembl cDNAs in FASTA format

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script dumps all Ensembl cDNAs in FASTA format on STDOUT. Optionally, the
3' end can be extended.

=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-75/ensembl/modules \
        dump_ensembl_cdna.pl \
        > transcriptome.fa

=head1 USAGE

    dump_ensembl_cdna.pl
        [--three_prime_extension bp]
        [--species species]
        [--ensembl_dbhost host]
        [--ensembl_dbport port]
        [--ensembl_dbuser username]
        [--ensembl_dbpass password]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--three_prime_extension INT>

Maximum number of base pairs to extend 3' end by.

=item B<--species SPECIES>

Species (defaults to Danio rerio).

=item B<--ensembl_dbhost HOST>

Ensembl MySQL database host.

=item B<--ensembl_dbport PORT>

Ensembl MySQL database port.

=item B<--ensembl_dbuser USERNAME>

Ensembl MySQL database username.

=item B<--ensembl_dbpass PASSWORD>

Ensembl MySQL database password.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=head1 DEPENDENCIES

Ensembl Perl API - http://www.ensembl.org/info/docs/api/

=head1 AUTHOR

=over 4

=item *

Ian Sealy <ian.sealy@sanger.ac.uk>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2014 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
