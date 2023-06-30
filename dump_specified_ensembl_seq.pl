#!/usr/bin/env perl

# PODNAME: dump_specified_ensembl_seq.pl
# ABSTRACT: Dump sequence from Ensembl in FASTA format for specified regions

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2014-11-24

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
my $species        = 'Danio rerio';
my $ensembl_dbhost = 'ensembldb.ensembl.org';
my $ensembl_dbport;
my $ensembl_dbuser = 'anonymous';
my $ensembl_dbpass;
my @regions;
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

# Dump FASTA to STDOUT
my $fasta_out = Bio::SeqIO->new(
    -fh     => \*STDOUT,
    -format => 'fasta',
);

# Dump sequence for each specified region
foreach my $region (@regions) {
    if ( $region !~ m/\A [^:]+ : \d+ : \d+ : -?1 \z/xms ) {
        croak "$region is not a valid region (chr:start:stop:strand)";
    }
    my ( $chr, $start, $end, $strand ) = split /:/xms, $region;
    my $slice = $sa->fetch_by_region( 'toplevel', $chr, $start, $end, $strand );
    my $seq = Bio::Seq->new(
        -id       => $region,
        -alphabet => 'dna',
        -seq      => $slice->seq,
    );
    $fasta_out->write_seq($seq);
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
        'region=s{1,}'     => \@regions,
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

dump_specified_ensembl_seq.pl

Dump sequence from Ensembl in FASTA format for specified regions

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script dumps sequence from Ensembl in FASTA format on STDOUT for specific
regions.

=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-75/ensembl/modules \
        dump_specified_ensembl_seq.pl \
        --region 1:100000-200000:1 5:200000-300000:-1 \
        > regions.fa

=head1 USAGE

    dump_specified_ensembl_seq.pl
        [--species species]
        [--ensembl_dbhost host]
        [--ensembl_dbport port]
        [--ensembl_dbuser username]
        [--ensembl_dbpass password]
        [--region region...]
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

=item B<--region REGION...>

The region (or regions) for which sequence is required. In the form
"chr:start:stop:strand", where strand is 1 or -1.

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
