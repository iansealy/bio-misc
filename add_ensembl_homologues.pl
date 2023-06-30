#!/usr/bin/env perl

# PODNAME: add_ensembl_homologues.pl
# ABSTRACT: Add Ensembl homologues to table of data

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2015-02-16

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use List::Util qw( sum );
use Bio::EnsEMBL::Registry;

# Default options
my $species = 'Danio rerio';
my @homologous_species;
my $gene_field     = 1;
my $has_header     = 0;
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
my $gma =
  Bio::EnsEMBL::Registry->get_adaptor( 'multi', 'compara', 'GeneMember' );
my $ha = Bio::EnsEMBL::Registry->get_adaptor( 'multi', 'compara', 'Homology' );

# Ensure database connection isn't lost; Ensembl 64+ can do this more elegantly
## no critic (ProhibitMagicNumbers)
if ( Bio::EnsEMBL::ApiVersion::software_version() < 64 ) {
## use critic
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
}
else {
    Bio::EnsEMBL::Registry->set_reconnect_when_lost();
}

# Print header line if present
if ($has_header) {
    my $line = <>;
    chomp $line;
    printf "%s\t%s\n", $line, join "\t",
      map { $_ . ' Ensembl ID' } @homologous_species;
}

# Iterate over STDIN and add homologues
while ( my $line = <> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;
    my $gene_id = $fields[ $gene_field - 1 ];
    my @all_homologous_genes;
    my $member = $gma->fetch_by_source_stable_id( 'ENSEMBLGENE', $gene_id );
    for my $homologous_species (@homologous_species) {
        if ( !$member ) {
            push @all_homologous_genes, q{};
        }
        else {
            my @homologous_genes;
            my $homologies =
              $ha->fetch_all_by_Member_paired_species( $member,
                $homologous_species );
            foreach my $homology ( @{$homologies} ) {
                foreach my $homology_member ( @{ $homology->gene_list() } ) {
                    next if $homology_member->stable_id eq $gene_id;
                    push @homologous_genes, $homology_member->stable_id;
                }
            }
            push @all_homologous_genes, join q{,}, sort @homologous_genes;
        }
    }
    printf "%s\t%s\n", $line, join "\t", @all_homologous_genes;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'species=s'                 => \$species,
        'homologous_species=s@{1,}' => \@homologous_species,
        'gene_field=i'              => \$gene_field,
        'header'                    => \$has_header,
        'ensembl_dbhost=s'          => \$ensembl_dbhost,
        'ensembl_dbport=i'          => \$ensembl_dbport,
        'ensembl_dbuser=s'          => \$ensembl_dbuser,
        'ensembl_dbpass=s'          => \$ensembl_dbpass,
        'debug'                     => \$debug,
        'help'                      => \$help,
        'man'                       => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !@homologous_species ) {
        pod2usage("--homologous_species must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

add_ensembl_homologues.pl

Add Ensembl homologues to table of data

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a tab-delimited file on STDIN and adds new final columns
containing the Ensembl IDs of genes in the requested species homologous to the
gene in the specified column.

=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-75/ensembl/modules \
        add_ensembl_homologues.pl \
        --species "Danio rerio" \
        --homologous_species "Homo sapiens" "Mus musculus" \
        --gene_field 1 \
        --header \
        < input.tsv > output.tsv

=head1 USAGE

    add_ensembl_homologues.pl
        [--species species]
        [--homologous_species species...]
        [--gene_field int]
        [--header]
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

The species whose Ensembl gene IDs are in the input file.

=item B<--homologous_species SPECIES...>

The species whose Ensembl gene IDs should be added to the output file.

=item B<--gene_field INT>

The field containing the Ensembl gene IDs in the input file.

=item B<--header>

Input has a header line.

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

None

=head1 AUTHOR

=over 4

=item *

Ian Sealy <ian.sealy@sanger.ac.uk>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2015 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
