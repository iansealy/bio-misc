#!/usr/bin/env perl

# PODNAME: run_deseq2_rnaseq.pl
# ABSTRACT: Run DESeq2 on RNA-Seq counts

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2017-05-16

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use English qw( -no_match_vars );
use POSIX qw( WIFEXITED);
use File::Spec;
use File::Path qw( make_path );

# Default options
my $counts_file  = 'deseq2/counts.txt';
my $samples_file = 'deseq2/samples.txt';
my $output_dir   = 'deseq2';
my @comparisons;
my $remove_other_conditions;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Get samples
my %condition_for;
my %is_condition;
my %group_for;
my %is_group;
my @all_samples;
open my $samples_fh, '<', $samples_file;    ## no critic (RequireBriefOpen)
my $header = <$samples_fh>;
while ( my $line = <$samples_fh> ) {
    chomp $line;
    my ( $sample, $condition, $group ) = split /\t/xms, $line;
    $condition_for{$sample}   = $condition;
    $is_condition{$condition} = 1;
    if ($group) {
        $group_for{$sample} = $group;
        $is_group{$group}   = 1;
    }
    push @all_samples, $sample;
}
close $samples_fh;

# Remove common prefix from conditions
my $prefix = get_common_prefix( keys %is_condition );
%is_condition = ();
foreach my $sample ( keys %condition_for ) {
    $condition_for{$sample} =~ s/\A $prefix //xms;
    $is_condition{ $condition_for{$sample} } = 1;
}

my @all_conditions = sort keys %is_condition;
my @all_groups     = sort keys %is_group;

# Get counts
my @genes;
my %counts_for;
open my $counts_fh, '<', $counts_file;    ## no critic (RequireBriefOpen)
$header = <$counts_fh>;
chomp $header;
my ( undef, @count_samples ) = split /\t/xms, $header;
while ( my $line = <$counts_fh> ) {
    chomp $line;
    my ( $gene, @counts ) = split /\t/xms, $line;
    push @genes, $gene;
    foreach my $i ( 0 .. ( scalar @count_samples ) - 1 ) {
        push @{ $counts_for{ $count_samples[$i] } }, $counts[$i];
    }
}

# Assume some comparisons
if ( !@comparisons ) {
    if ( scalar @all_conditions == 1 ) {
        confess "Only one condition (@all_conditions)";
    }

    my ($wt)  = grep { m/(\b|_)wt \z/xms } @all_conditions;
    my ($het) = grep { m/(\b|_)het \z/xms } @all_conditions;
    my ($hom) = grep { m/(\b|_)hom \z/xms } @all_conditions;
    my ($sib) = grep { m/(\b|_)sib \z/xms } @all_conditions;
    my ($mut) = grep { m/(\b|_)mut \z/xms } @all_conditions;

    if ( scalar @all_conditions == 2 ) {
        @comparisons =
            $wt  && $het ? ("$het:$wt")
          : $wt  && $hom ? ("$hom:$wt")
          : $het && $hom ? ("$hom:$het")
          : $sib && $mut ? ("$mut:$sib")
          :                ();
    }
    ## no critic (ProhibitMagicNumbers)
    elsif ( scalar @all_conditions == 3 && $wt && $het && $hom ) {
        ## use critic
        push @comparisons, "$het:$wt";
        push @comparisons, "$hom:$wt";
        push @comparisons, "$hom:$het";
        push @comparisons, "$hom:$het,$wt";
        push @comparisons, "$hom,$het:$wt";
    }
}

foreach my $comparison (@comparisons) {
    my ( $all_exp, $all_con ) = split /:/xms, $comparison;
    confess "Experimental condition missing from $comparison" if !$all_exp;
    confess "Control condition missing from $comparison"      if !$all_con;
    my ( $exp, $exp_name ) = split /=/xms, $all_exp;
    my ( $con, $con_name ) = split /=/xms, $all_con;
    if ( !$exp_name ) {
        $exp_name = $exp;
        $exp_name =~ s/,/_/xmsg;
    }
    if ( !$con_name ) {
        $con_name = $con;
        $con_name =~ s/,/_/xmsg;
    }
    my @exp = split /,/xms, $exp;
    my @con = split /,/xms, $con;
    my %rename;
    foreach my $condition (@exp) {
        confess "Unknown condition ($condition) in $comparison"
          if !$is_condition{$condition};
        $rename{$condition} = $exp_name;
    }
    foreach my $condition (@con) {
        confess "Unknown condition ($condition) in $comparison"
          if !$is_condition{$condition};
        $rename{$condition} = $con_name;
    }
    my $dir = File::Spec->catdir( $output_dir, $exp_name . '_vs_' . $con_name );
    next if -e ( $dir . '.done' );
    make_path($dir);

    # Write new samples file
    my $new_samples_file = File::Spec->catfile( $dir, 'samples.txt' );
    open $samples_fh, '>', $new_samples_file;    ## no critic (RequireBriefOpen)
    printf {$samples_fh} "\tcondition%s\n", ( @all_groups ? "\tgroup" : q{} );
    foreach my $sample (@all_samples) {
        my $condition = $condition_for{$sample};
        next if $remove_other_conditions && !exists $rename{$condition};

        printf {$samples_fh} "%s\t%s%s\n", $sample,
          $rename{ $condition_for{$sample} },
          ( @all_groups ? "\t" . $group_for{$sample} : q{} );
    }
    close $samples_fh;

    # Write new counts file
    my $new_counts_file = File::Spec->catfile( $dir, 'counts.txt' );
    open $counts_fh, '>', $new_counts_file;    ## no critic (RequireBriefOpen)
    printf {$counts_fh} "\t%s\n", ( join "\t", @all_samples );
    foreach my $i ( 0 .. ( scalar @genes ) - 1 ) {
        my @counts;
        foreach my $sample (@all_samples) {
            my $condition = $condition_for{$sample};
            next if $remove_other_conditions && !exists $rename{$condition};
            push @counts, $counts_for{$sample}->[$i];
        }
        printf {$counts_fh} "%s\t%s\n", $genes[$i], ( join "\t", @counts );
    }
    close $counts_fh;

    # Write R script
    my $design = (@all_groups) ? 'group + condition' : 'condition';
    ## no critic (RequireBriefOpen)
    open my $r_fh, '>', File::Spec->catfile( $dir, 'deseq2.R' );
    ## use critic
    print {$r_fh} <<"EOF";    ## no critic (RequireCheckedSyscalls)
suppressWarnings(library(tcltk))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gplots))
countData <- read.table( "$new_counts_file", header=TRUE, row.names=1 )
samples <- read.table( "$new_samples_file", header=TRUE, row.names=1 )
dds <- DESeqDataSetFromMatrix(countData, samples, design = ~ $design)
dds <- DESeq(dds)
write.table(sizeFactors(dds), file="$dir/size-factors.txt", col.names=FALSE, quote=FALSE, sep="\\t")
write.table(counts(dds, normalized=TRUE), file="$dir/normalised-counts.txt", col.names=FALSE, quote=FALSE, sep="\\t")
res <- results(dds, contrast=c("condition", "$exp_name", "$con_name"))
out <- data.frame(pvalue=res\$pvalue, padj=res\$padj, log2fc=res\$log2FoldChange, row.names=rownames(res))
write.table(out, file="$dir/output.txt", col.names=FALSE, row.names=TRUE, quote=FALSE, sep="\\t")
pdf("$dir/qc.pdf")
plotMA(dds)
rld <- rlogTransformation(dds, blind=TRUE)
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(assay(rld)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
heatmap.2(as.matrix(dist(t(assay(rld)))), trace="none", col = rev(hmcol), margin=c(13, 13))
print(plotPCA(rld, intgroup=c("condition")))
plotDispEsts(dds)
dev.off()
quit()
EOF
    close $r_fh;

    # Run R script under LSF
    printf "Running %s\n", "$dir/deseq2.R";
    my $cmd = <<"EOF";
bsub -o $dir/deseq2.o -e $dir/deseq2.e -R'select[mem>4000] rusage[mem=4000]' -M4000 "Rscript $dir/deseq2.R"
EOF
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";
}

sub get_common_prefix {
    my @strings = @_;

    my $common_prefix_idx = 0;
  CHAR: while (1) {
        my $char = substr $strings[0], $common_prefix_idx, 1;
        foreach my $string (@strings) {
            last CHAR if substr( $string, $common_prefix_idx, 1 ) ne $char;
        }
        $common_prefix_idx++;
    }

    return substr $strings[0], 0, $common_prefix_idx;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'counts_file=s'           => \$counts_file,
        'samples_file=s'          => \$samples_file,
        'output_dir'              => \$output_dir,
        'comparisons=s@{,}'       => \@comparisons,
        'remove_other_conditions' => \$remove_other_conditions,
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

run_deseq2_rnaseq.pl

Run DESeq2 on RNA-seq counts

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes an RNA-Seq counts file and samples file and runs DESeq2 using
specified or all conditions.

=head1 EXAMPLES

    perl \
        run_deseq2_rnaseq.pl \
        --counts_file counts.txt --samples_file samples.txt

    perl \
        run_deseq2_rnaseq.pl \
        --counts_file counts.txt --samples_file samples.txt \
        --output_dir deseq2 \
        --comparisons hom:het hom:wt het:wt hom=mut:het,wt=sib \
        --remove_other_conditions

=head1 USAGE

    convert_to_biolayout.pl
        [--counts_file file]
        [--samples_file file]
        [--output_dir dir]
        [--comparisons comparison...]
        [--remove_other_conditions]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--counts_file FILE>

RNA-Seq counts file (e.g. counts.txt).

=item B<--samples_file FILE>

DESeq2 samples file (e.g. samples.txt). The order of samples in the samples file
determines the order of the columns in the output.

=item B<--output_dir DIR>

Directory in which to create output directories.

=item B<--comparisons COMPARISONS>

Condition comparisons. Each comparison is a pair of exerpimental and controls
conditions (in that order) separated by a colon (e.g. hom:wt). If multiple
conditions are to be combined then separate them with a comma (e.g. hom:wt,het).
To rename a condition, append an equals sign (e.g. hom=mut:het,wt=sib).

=item B<--remove_other_conditions>

Remove other conditions from the counts file prior to running DESeq2.

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