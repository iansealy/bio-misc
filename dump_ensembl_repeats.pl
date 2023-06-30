#!/usr/bin/env perl

# PODNAME: dump_ensembl_repeats.pl
# ABSTRACT: Dump Ensembl repeats in BED format

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2014-10-27

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Sort::Naturally;
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

    # Get all repeats
    my $repeats = $slice->get_all_RepeatFeatures();
    foreach my $repeat ( @{$repeats} ) {
        my $name = join q{,}, $repeat->repeat_consensus->repeat_type,
          $repeat->repeat_consensus->repeat_class,
          $repeat->repeat_consensus->name;
        my $strand = q{.};
        if ( $repeat->seq_region_strand == -1 ) {
            $strand = q{-};
        }
        elsif ( $repeat->seq_region_strand == 1 ) {
            $strand = q{+};
        }
        printf "%s\t%d\t%d\t%s\t0\t%s\n", $repeat->seq_region_name,
          $repeat->seq_region_start, $repeat->seq_region_end, $name, $strand;
    }
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

dump_ensembl_repeats.pl

Dump Ensembl repeats in BED format

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script dumps all Ensembl repeats in BED format on STDOUT.

=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-75/ensembl/modules \
        dump_ensembl_repeats.pl \
        > repeats.bed

=head1 USAGE

    dump_ensembl_repeats.pl
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
