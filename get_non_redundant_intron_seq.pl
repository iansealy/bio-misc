#!/usr/bin/env perl

# PODNAME: get_non_redundant_intron_seq.pl
# ABSTRACT: Get non-redundant intron sequences for list of Ensembl genes

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2014-03-07

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Bio::EnsEMBL::Registry;

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

# Iterate over each gene
while ( my $ensembl_id = <> ) {
    chomp $ensembl_id;
    my $gene = $ga->fetch_by_stable_id($ensembl_id);

    if ( !$gene ) {
        printf {*STDERR} "Can't get gene for %s\n", $ensembl_id;
        next;
    }

    # Get coordinates for all introns of all transcripts of a gene
    my @intron_coords;
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript ( @{$transcripts} ) {
        my $introns = $transcript->get_all_Introns();
        foreach my $intron ( @{$introns} ) {
            push @intron_coords,
              [ $intron->seq_region_start, $intron->seq_region_end ];
        }
    }

    # Check for genes without introns
    if ( !@intron_coords ) {
        printf {*STDERR} "No introns for %s\n", $ensembl_id;
        next;
    }

    # Collapse to non-overlapping set of intron coordinates
    @intron_coords =
      sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @intron_coords;
    my @nr_intron_coords;
    my $current_min_start = $intron_coords[0]->[0];
    my $current_max_end   = $intron_coords[0]->[1];
    foreach my $coords (@intron_coords) {
        my ( $start, $end ) = @{$coords};
        if ( $start > $current_max_end && $start - $current_max_end != 1 ) {

            # Introns don't overlap and aren't adjacent
            push @nr_intron_coords, [ $current_min_start, $current_max_end ];
            $current_min_start = $start;
            $current_max_end   = $end;
        }
        elsif ( $end > $current_max_end ) {

            # Extend intron
            $current_max_end = $end;
        }
    }
    push @nr_intron_coords, [ $current_min_start, $current_max_end ];

    # Reverse introns if on reverse strand
    if ( $gene->seq_region_strand == -1 ) {
        @nr_intron_coords =
          reverse sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
          @nr_intron_coords;
    }

    # Get sequence for each intron and concatenate
    my $concat_seq;
    foreach my $coords (@nr_intron_coords) {
        my ( $start, $end ) = @{$coords};
        my $slice =
          $sa->fetch_by_region( 'toplevel', $gene->seq_region_name, $start,
            $end, $gene->seq_region_strand );
        $concat_seq .= $slice->seq;
    }

    # Output sequence
    printf ">%s\n%s\n", $ensembl_id, $concat_seq;
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

get_non_redundant_intron_seq.pl

Get non-redundant intron sequences for list of Ensembl genes

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a list of Ensembl gene IDs on STDIN (one per line) and outputs
non-redundant intron sequence for each gene in FASTA format on STDOUT.

=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-75/ensembl/modules \
        get_non_redundant_intron_seq.pl \
        < ensembl_ids.txt > introns.fa

=head1 USAGE

    get_non_redundant_intron_seq.pl
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
