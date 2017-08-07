#!/usr/bin/env perl

# PODNAME: convert_two_col_to_upset.pl
# ABSTRACT: Convert two columns of data into UpSet format

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2017-08-07

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

# Default options
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

my %upset;
my %is_set;

while ( my $line = <> ) {
    chomp $line;
    ## no critic (ProhibitAmbiguousNames)
    my ( $set, $id ) = split /\s+/xms, $line;
    ## use critic
    $upset{$id}{$set} = 1;
    $is_set{$set} = 1;
}
my @sets = sort keys %is_set;

printf "%s\r\n", join q{,}, 'Gene', @sets;
foreach my $id ( sort keys %upset ) {
    printf "%s\r\n", join q{,}, $id, map { $upset{$id}{$_} || 0 } @sets;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'debug' => \$debug,
        'help'  => \$help,
        'man'   => \$man,
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

convert_two_col_to_upset.pl

Convert two columns of data into UpSet format

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes two column input (sets and IDs) on STDIN and converts to the
format required by UpSet.

=head1 EXAMPLES

    perl convert_two_col_to_upset.pl < cols.tsv > upset.csv

=head1 USAGE

    convert_two_col_to_upset.pl
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

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

This software is Copyright (c) 2017 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
