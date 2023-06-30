#!/usr/bin/env perl

# PODNAME: classify_genome_bases.pl
# ABSTRACT: Classify every base of a genome according to gene distance

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2014-10-19

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Readonly;
use List::Util qw( max );
use Sort::Naturally;
use Bio::EnsEMBL::Registry;
use Set::IntervalTree;

# Constants
Readonly our %RANK => (
    tss      => 5,
    internal => 4,
    promoter => 3,
    proximal => 2,
    distal   => 1,
    desert   => 0,
);
my %CLASS = reverse %RANK;

# Default options
my $species        = 'Danio rerio';
my $ensembl_dbhost = 'ensembldb.ensembl.org';
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

# Get Ensembl adaptor
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

# Get all slices
my $slices = $sa->fetch_all('toplevel');
warn scalar @{$slices}, " slices\n" if $debug;

foreach my $slice ( sort { ncmp( $a->seq_region_name, $b->seq_region_name ) }
    @{$slices} )
{
    warn 'Slice: ', $slice->name, "\n" if $debug;

    my $tree = Set::IntervalTree->new;

    # Get all genes
    my $genes = $slice->get_all_Genes( undef, 'core', 1 );    # Plus transcripts
    warn scalar @{$genes}, " genes\n" if $debug;

    foreach my $gene ( @{$genes} ) {
        warn ' Gene: ', $gene->stable_id, ' / ', $gene->biotype, ' / ',
          $gene->strand, "\n"
          if $debug;

        # Get all transcripts of gene
        my $transcripts = $gene->get_all_Transcripts();
        foreach my $transcript ( @{$transcripts} ) {
            warn '  Transcript: ', $transcript->stable_id, ' / ',
              $transcript->biotype, ' / ', $transcript->strand, "\n"
              if $debug;

            my $start  = $transcript->seq_region_start;
            my $end    = $transcript->seq_region_end;
            my $strand = $transcript->seq_region_strand;
            if ( $strand < 0 ) {
                ( $start, $end ) = ( $end, $start );
            }

            warn q{   }, $start, q{-}, $end, "\n" if $debug;

            insert(
                $tree, 'tss',
                $start - 1000 * $strand,
                $start + 1000 * $strand
            );
            if ( abs $end - $start > 1000 ) {

                # Skip internal if gene <= 1000 bp (i.e. mark gene as tss)
                insert( $tree, 'internal', $start + 1001 * $strand, $end );
            }
            insert(
                $tree, 'promoter',
                $start - 5000 * $strand,
                $start - 1001 * $strand
            );
            insert(
                $tree, 'proximal',
                $start - 10_000 * $strand,
                $start - 5001 * $strand
            );
            insert(
                $tree, 'proximal',
                $end + 1 * $strand,
                $end + 10_000 * $strand
            );
            insert(
                $tree, 'distal',
                $start - 100_000 * $strand,
                $start - 10_001 * $strand
            );
            insert(
                $tree, 'distal',
                $end + 10_001 * $strand,
                $end + 100_000 * $strand
            );
        }
    }

    # Output in BED format
    my $current_class = fetch( $tree, 1 );
    my $current_start = 1;
    foreach my $base ( 2 .. $slice->seq_region_length ) {
        my $class = fetch( $tree, $base );
        if ( $class ne $current_class ) {
            printf "%s\t%d\t%d\t%s\n", $slice->seq_region_name,
              $current_start - 1, $base - 1, $current_class;
            $current_start = $base;
            $current_class = $class;
        }
    }
    printf "%s\t%d\t%d\t%s\n", $slice->seq_region_name, $current_start - 1,
      $slice->seq_region_length, $current_class;
}

# Insert into interval tree
sub insert {
    my ( $tree, $class, $start, $end ) = @_;

    if ( $start > $end ) {
        ( $start, $end ) = ( $end, $start );
    }

    warn '    Inserting ', $class, ' from ', $start, ' to ', $end, "\n"
      if $debug;

    # Convert from Ensembl closed interval to half-open interval
    $tree->insert( $RANK{$class}, $start, $end + 1 );

    return;
}

# Fetch from interval tree and return highest rank class
sub fetch {
    my ( $tree, $base ) = @_;

    # Convert from Ensembl closed interval to half-open interval
    my $ranks = $tree->fetch( $base, $base + 1 );

    return $CLASS{ ( max @{$ranks} ) || 0 };
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'species=s'        => \$species,
        'ensembl_dbhost=s' => \$ensembl_dbhost,
        'ensembl_dbport=i' => \$ensembl_dbport,
        'ensembl_dbuser=s' => \$ensembl_dbuser,
        'ensembl_dbpass=s' => \$ensembl_dbpass,
        'debug'            => \$debug,
        'help'             => \$help,
        'man'              => \$man,
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

classify_genome_bases.pl

Classify every base of a genome according to gene distance

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script classifies every base of a genome according to distance from a gene.
The classification is output in BED format. The classifications (in order of
precedence) are:

=over 8

=item I<tss>

From -1000 bp to +1000 bp from a transcription start site (TSS).

=item I<internal>

From +1001 bp from a TSS to the end of the gene.

=item I<promoter>

From -5000 bp from a TSS to -1001 bp.

=item I<proximal>

From -10,000 bp from a TSS to -5001 bp or from +1 bp from the end of the gene to
+10,000 bp from the end of the gene.

=item I<distal>

From -100,000 bp from a TSS to -10,001 bp or from +10,001 bp beyond the end of
the gene to +100,000 bp.

=item I<desert>

More than 100,000 bp from the TSS or end of the gene.

=back

=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-77/ensembl/modules \
        classify_genome_bases.pl \
        > ens77classification.bed

=head1 USAGE

    classify_genome_bases.pl
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
Set::IntervalTree

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
