#!/usr/bin/env perl

# PODNAME: add_ensembl_gene_annotation.pl
# ABSTRACT: Add annotation for Ensembl genes

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2016-01-06

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Bio::EnsEMBL::Registry;

# Default options
my $prefix         = 'ENSDARG';
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

# Ensure database connection isn't lost; Ensembl 64+ can do this more elegantly
## no critic (ProhibitMagicNumbers)
if ( Bio::EnsEMBL::ApiVersion::software_version() < 64 ) {
## use critic
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
}
else {
    Bio::EnsEMBL::Registry->set_reconnect_when_lost();
}

# Iterate over STDIN and add annotation
while ( my $line = <> ) {
    chomp $line;
    my ( @chrs, @starts, @ends, @strands, @biotypes, @names, @descriptions );
    while ( $line =~ m/($prefix\d+)/xmsg ) {
        my $gene = $ga->fetch_by_stable_id($1);
        if ($gene) {
            push @chrs,     $gene->seq_region_name;
            push @starts,   $gene->seq_region_start;
            push @ends,     $gene->seq_region_end;
            push @strands,  $gene->seq_region_strand;
            push @biotypes, $gene->biotype;
            push @names, ( $gene->external_name || $gene->stable_id );
            push @descriptions, ( $gene->description || q{} );
        }
    }
    if ( !scalar @biotypes ) {
        @chrs         = (q{});
        @starts       = (q{});
        @ends         = (q{});
        @strands      = (q{});
        @biotypes     = (q{});
        @names        = (q{});
        @descriptions = (q{});
    }
    my $chr         = join q{,}, @chrs;
    my $start       = join q{,}, @starts;
    my $end         = join q{,}, @ends;
    my $strand      = join q{,}, @strands;
    my $biotype     = join q{,}, @biotypes;
    my $name        = join q{,}, @names;
    my $description = join q{,}, @descriptions;
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $line, $chr, $start, $end,
      $strand, $biotype, $name, $description;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'prefix=s'         => \$prefix,
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

add_ensembl_gene_annotation.pl

Add annotation for Ensembl genes

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes input containing Ensembl gene stable IDs and adds some
annotation.

=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-75/ensembl/modules \
        add_ensembl_gene_annotation.pl < input.txt

=head1 USAGE

    add_ensembl_gene_annotation.pl
        [--prefix prefix]
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

=item B<--prefix PREFIX>

Gene ID prefix (defaults to ENSDARG).

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

This software is Copyright (c) 2017 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
