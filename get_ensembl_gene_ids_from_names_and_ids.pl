#!/usr/bin/env perl

# PODNAME: get_ensembl_gene_ids_from_names_and_ids.pl
# ABSTRACT: Convert list of names and IDs to Ensembl gene IDs

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

# Get Ensembl adaptor
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

# Cache IDs and names for all genes
my %id2id;
my $genes = $ga->fetch_all_by_biotype('protein_coding');
foreach my $gene ( @{$genes} ) {

    # Cache stable ID
    $id2id{ $gene->stable_id }{ $gene->stable_id } = 1;

    # Cache external name
    if ( $gene->external_name ) {
        $id2id{ $gene->external_name }{ $gene->stable_id } = 1;
    }

    # Cache ZFIN IDs
    my $entries = $gene->get_all_DBEntries('ZFIN_ID');
    foreach my $entry ( @{$entries} ) {
        $id2id{ $entry->primary_id }{ $gene->stable_id } = 1;
        $id2id{ $entry->display_id }{ $gene->stable_id } = 1;
    }
}

# Convert names and IDs on STDIN
my %ensembl_id;
while ( my $line = <> ) {
    chomp $line;
    my @ids = split /\s* , \s*/xms, $line;
    foreach my $id (@ids) {
        $id =~ s/\A \s+//xms;
        $id =~ s/\s+ \z//xms;
        if ( exists $id2id{$id} ) {
            foreach my $ensembl_id ( keys %{ $id2id{$id} } ) {
                $ensembl_id{$ensembl_id} = 1;
            }
        }
        else {
            printf {*STDERR} "Can't find Ensembl ID for %s\n", $id;
        }
    }
}

# Output Ensembl IDs
foreach my $ensembl_id ( sort keys %ensembl_id ) {
    printf "%s\n", $ensembl_id;
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

get_ensembl_gene_ids_from_names_and_ids.pl

Convert list of names and IDs to Ensembl gene IDs

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a list of gene names and/or ZFIN IDs (one per line and/or
separated with commas) on STDIN and converts them to a list of Ensembl gene IDs
on STDOUT.

=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-75/ensembl/modules \
        get_ensembl_gene_ids_from_names_and_ids.pl \
        < genes.txt > ensembl_ids.txt

=head1 USAGE

    get_ensembl_gene_ids_from_names_and_ids.pl
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
