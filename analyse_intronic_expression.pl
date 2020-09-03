#!/usr/bin/env perl

# PODNAME: analyse_intronic_expression.pl
# ABSTRACT: Analyse intronic expression

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2018-02-28

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Sort::Naturally;
use Bio::EnsEMBL::Registry;
use Bio::DB::HTS;
use Readonly;

# Constants
Readonly our $MAPQ_THRESHOLD => 10;

# Default options
my $bam_file;
my $repeats_file;
my $perfect_matches;
my $slice_regexp;
my $species        = 'Mus musculus';
my $ensembl_dbhost = 'ensembldb.ensembl.org';
my $ensembl_dbport;
my $ensembl_dbuser = 'anonymous';
my $ensembl_dbpass;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Connnect to Ensembl database
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

my $hts = Bio::DB::HTS->new( -bam => $bam_file );

printf "%s\n", join "\t", 'Gene ID', 'Gene Biotype', 'Gene Overlap Count',
  'Transcript ID', 'Transcript Biotype', 'Transcript Overlap Count', 'Chr',
  'Strand', 'Intron Start', 'Intron End', 'Intron Enclosed Count',
  'Intron Enclosed FPKM', 'Prev Exon ID', 'Prev Exon Start', 'Prev Exon End',
  'Prev Exon Overlap Count', 'Prev Exon Overlap FPKM', 'Next Exon ID',
  'Next Exon Start',        'Next Exon End', 'Next Exon Overlap Count',
  'Next Exon Overlap FPKM', 'Repeat Names',  'Repeat Coords',
  'Repeat Overlap Counts',  'Repeat Overlap FPKMs',
  'Repeat Overlap Total Count',    'Repeat Overlap Total FPKM',
  'Intron Segment Coords',         'Intron Segment Enclosed Counts',
  'Intron Segment Enclosed FPKMs', 'Intron Segment Enclosed Total Count',
  'Intron Segment Enclosed Total FPKM';
my $slices = $sa->fetch_all('toplevel');
warn scalar @{$slices}, " slices\n" if $debug;
my $total_count = 0;
foreach my $slice ( sort { ncmp( $a->seq_region_name, $b->seq_region_name ) }
    @{$slices} )
{
    warn 'Slice: ', $slice->name, "\n" if $debug;
    my $genes = $slice->get_all_Genes;
    foreach my $gene ( @{$genes} ) {
        my $gene_count = count_overlapping( $hts, $gene->seq_region_name,
            $gene->seq_region_start, $gene->seq_region_end, $gene->strand );
        $total_count += $gene_count;
    }
}
warn 'Total count: ', $total_count, "\n" if $debug;
foreach my $slice ( sort { ncmp( $a->seq_region_name, $b->seq_region_name ) }
    @{$slices} )
{
    next
      if defined $slice_regexp
      && $slice->seq_region_name !~ m/$slice_regexp/xms;
    warn 'Slice: ', $slice->name, "\n" if $debug;
    my $repeat_cache = load_repeats( $repeats_file, $slice->seq_region_name );
    my $genes = $slice->get_all_Genes;
    foreach my $gene ( @{$genes} ) {
        ## no critic (RequireCarping)
        warn sprintf "Gene: %s %s:%d-%d:%d\n", $gene->stable_id,
          $gene->seq_region_name, $gene->seq_region_start,
          $gene->seq_region_end,  $gene->seq_region_strand
          if $debug;
        ## use critic
        my $gene_count = count_overlapping( $hts, $gene->seq_region_name,
            $gene->seq_region_start, $gene->seq_region_end, $gene->strand );
        printf "# %s\n", join "\t", $gene->stable_id, $gene->biotype,
          $gene_count;
        my $transcripts = $gene->get_all_Transcripts();
        foreach my $transcript (
            sort {
                     $a->seq_region_start <=> $b->seq_region_start
                  || $a->seq_region_end <=> $b->seq_region_end
            } @{$transcripts}
          )
        {
            ## no critic (RequireCarping)
            warn sprintf "Transcript: %s %s:%d-%d:%d\n",
              $transcript->stable_id,        $transcript->seq_region_name,
              $transcript->seq_region_start, $transcript->seq_region_end,
              $transcript->seq_region_strand
              if $debug;
            ## use critic
            my $transcript_count = count_overlapping(
                $hts,
                $transcript->seq_region_name,
                $transcript->seq_region_start,
                $transcript->seq_region_end,
                $transcript->strand
            );
            printf "## %s\n", join "\t", $gene->stable_id, $gene->biotype,
              $gene_count, $transcript->stable_id, $transcript->biotype,
              $transcript_count;
            my $introns = $transcript->get_all_Introns();
            if ( !@{$introns} ) {
                printf "### %s\tno-introns\n", join "\t", $gene->stable_id,
                  $gene->biotype, $gene_count, $transcript->stable_id,
                  $transcript->biotype, $transcript_count;
                next;
            }
            foreach my $intron ( @{$introns} ) {
                ## no critic (RequireCarping)
                warn sprintf "Intron: %s:%d-%d:%d\n", $intron->seq_region_name,
                  $intron->seq_region_start, $intron->seq_region_end,
                  $intron->seq_region_strand
                  if $debug;
                ## use critic
                my ( $intron_count, $intron_fpkm ) = count_enclosed(
                    $hts,                       $slice->seq_region_name,
                    $intron->seq_region_start,  $intron->seq_region_end,
                    $intron->seq_region_strand, $total_count
                );
                my ( $prev_exon_count, $prev_exon_fpkm ) = count_overlapping(
                    $hts,
                    $slice->seq_region_name,
                    $intron->prev_Exon->seq_region_start,
                    $intron->prev_Exon->seq_region_end,
                    $intron->seq_region_strand,
                    $total_count
                );
                my ( $next_exon_count, $next_exon_fpkm ) = count_overlapping(
                    $hts,
                    $slice->seq_region_name,
                    $intron->next_Exon->seq_region_start,
                    $intron->next_Exon->seq_region_end,
                    $intron->seq_region_strand,
                    $total_count
                );
                my $repeats = get_overlapping_repeats(
                    $repeat_cache->{ $intron->seq_region_strand },
                    $intron->seq_region_start, $intron->seq_region_end );
                foreach my $repeat (
                    sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
                    @{$repeats} )
                {
                    ## no critic (RequireCarping)
                    warn sprintf "Uncollapsed Repeat: %s:%d-%d:%d\n",
                      $intron->seq_region_name, $repeat->[0], $repeat->[1],
                      $intron->seq_region_strand
                      if $debug;
                    ## use critic
                }
                $repeats =
                  collapse_repeats( $repeats, $intron->seq_region_start,
                    $intron->seq_region_end );
                my $repeat_total_count = 0;
                my $repeat_total_fpkm  = 0;
                my @repeat_counts;
                my @repeat_fpkms;
                my @repeat_names;
                my @repeat_coords;
                my $intron_segment_total_count = 0;
                my $intron_segment_total_fpkm  = 0;
                my @intron_segment_counts;
                my @intron_segment_fpkms;
                my @intron_segment_coords;
                my $intron_segment_start = $intron->seq_region_start;
                my $intron_segment_end;
                my $intron_segment_count;
                my $intron_segment_fpkm;
                my $prev_repeat_end;

                foreach my $repeat (
                    sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
                    @{$repeats} )
                {
                    ## no critic (RequireCarping)
                    warn sprintf "Collapsed Repeat: %s:%d-%d:%d\n",
                      $intron->seq_region_name, $repeat->[0], $repeat->[1],
                      $intron->seq_region_strand
                      if $debug;
                    ## use critic
                    my ( $repeat_count, $repeat_fpkm ) =
                      count_overlapping( $hts, $slice->seq_region_name,
                        $repeat->[0], $repeat->[1], $intron->seq_region_strand,
                        $total_count );
                    push @repeat_counts, $repeat_count;
                    push @repeat_fpkms,  $repeat_fpkm;
                    push @repeat_names,  $repeat->[2];
                    push @repeat_coords, $repeat->[0] . q{-} . $repeat->[1];
                    $intron_segment_end = $repeat->[0] - 1;

                    ## no critic (ProhibitDeepNests)
                    if ( $repeat->[0] > $intron->seq_region_start ) {
                        ## use critic
                        ( $intron_segment_count, $intron_segment_fpkm ) =
                          count_enclosed(
                            $hts,
                            $slice->seq_region_name,
                            $intron_segment_start,
                            $intron_segment_end,
                            $intron->seq_region_strand,
                            $total_count
                          );
                        push @intron_segment_counts, $intron_segment_count;
                        push @intron_segment_fpkms,  $intron_segment_fpkm;
                        push @intron_segment_coords,
                          $intron_segment_start . q{-} . $intron_segment_end;
                    }
                    $intron_segment_start = $repeat->[1] + 1;
                    $prev_repeat_end      = $repeat->[1];
                }
                $intron_segment_end = $intron->seq_region_end;
                if ( !defined $prev_repeat_end
                    || $prev_repeat_end < $intron->seq_region_end )
                {
                    ( $intron_segment_count, $intron_segment_fpkm ) =
                      count_enclosed(
                        $hts,                       $slice->seq_region_name,
                        $intron_segment_start,      $intron_segment_end,
                        $intron->seq_region_strand, $total_count
                      );
                    push @intron_segment_counts, $intron_segment_count;
                    push @intron_segment_fpkms,  $intron_segment_fpkm;
                    push @intron_segment_coords,
                      $intron_segment_start . q{-} . $intron_segment_end;
                }

                # Get totals
                if ( @{$repeats} ) {
                    my @repeat_pairs = map { [ $_->[0], $_->[1] ] } @{$repeats};
                    ( $repeat_total_count, $repeat_total_fpkm ) =
                      count_overlapping_multiple( $hts,
                        $slice->seq_region_name, \@repeat_pairs,
                        $intron->seq_region_strand, $total_count );
                }
                if (@intron_segment_coords) {
                    my @intron_segment_pairs =
                      map { [ split /-/xms ] } @intron_segment_coords;
                    ( $intron_segment_total_count, $intron_segment_total_fpkm )
                      = count_enclosed_multiple( $hts, $slice->seq_region_name,
                        \@intron_segment_pairs, $intron->seq_region_strand,
                        $total_count );
                }

                $intron_segment_total_count =
                  @intron_segment_counts ? $intron_segment_total_count : q{-};
                $intron_segment_total_fpkm =
                  @intron_segment_counts ? $intron_segment_total_fpkm : q{-};
                my $intron_segment_counts =
                  @intron_segment_counts
                  ? join q{,}, @intron_segment_counts
                  : q{-};
                my $intron_segment_fpkms =
                  @intron_segment_fpkms
                  ? join q{,}, @intron_segment_fpkms
                  : q{-};
                my $intron_segment_coords =
                  @intron_segment_coords
                  ? join q{,}, @intron_segment_coords
                  : q{-};
                $repeat_total_count =
                  @repeat_counts ? $repeat_total_count : q{-};
                $repeat_total_fpkm = @repeat_counts ? $repeat_total_fpkm : q{-};
                my $repeat_counts =
                  @repeat_counts
                  ? join q{,}, @repeat_counts
                  : q{-};
                my $repeat_fpkms =
                  @repeat_fpkms
                  ? join q{,}, @repeat_fpkms
                  : q{-};
                my $repeat_names =
                  @repeat_names
                  ? join q{,}, @repeat_names
                  : q{-};
                my $repeat_coords =
                  @repeat_coords
                  ? join q{,}, @repeat_coords
                  : q{-};
                printf "%s\n", join "\t", $gene->stable_id, $gene->biotype,
                  $gene_count, $transcript->stable_id, $transcript->biotype,
                  $transcript_count, $slice->seq_region_name,
                  $intron->seq_region_strand, $intron->seq_region_start,
                  $intron->seq_region_end, $intron_count, $intron_fpkm,
                  $intron->prev_Exon->stable_id,
                  $intron->prev_Exon->seq_region_start,
                  $intron->prev_Exon->seq_region_end, $prev_exon_count,
                  $prev_exon_fpkm, $intron->next_Exon->stable_id,
                  $intron->next_Exon->seq_region_start,
                  $intron->next_Exon->seq_region_end, $next_exon_count,
                  $next_exon_fpkm, $repeat_names, $repeat_coords,
                  $repeat_counts,  $repeat_fpkms, $repeat_total_count,
                  $repeat_total_fpkm,          $intron_segment_coords,
                  $intron_segment_counts,      $intron_segment_fpkms,
                  $intron_segment_total_count, $intron_segment_total_fpkm;
            }
        }
        $gene->flush_Transcripts();    # Save memory
    }
}

# Load repeats from RepeatMasker file
sub load_repeats {
    my ( $file, $chr ) = @_;

    my %cache;

    open my $fh, '<', $file;           ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    $header = <$fh>;
    $header = <$fh>;
    while ( my $line = <$fh> ) {
        chomp $line;
        $line =~ s/\A \s+//xms;        # Remove leading space
        my @fields = split /\s+/xms, $line;
        ## no critic (ProhibitMagicNumbers)
        next if $fields[4] ne $chr;
        my ( $start, $end, $strand, $name, $class_family ) =
          @fields[ 5, 6, 8, 9, 10 ];
        ## use critic
        $strand =~ s/[+]/1/xms;
        $strand =~ s/C/-1/xms;
        push @{ $cache{$strand} },
          [ $start, $end, $name . q{:} . $class_family ];
    }
    close $fh;

    return \%cache;
}

sub get_overlapping_repeats {
    my ( $cache, $start, $end ) = @_;

    my @repeats = grep { $_->[0] <= $end && $_->[1] >= $start } @{$cache};

    return \@repeats;
}

# Collapse overlapping and adjacent repeats
sub collapse_repeats {
    my ( $repeats, $region_start, $region_end ) = @_;

    return $repeats if !@{$repeats};

    my @collapsed_repeats;

    @{$repeats} =
      sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$repeats};
    my $first_repeat = shift @{$repeats};
    my $current_repeat_start =
      $first_repeat->[0] < $region_start ? $region_start : $first_repeat->[0];
    my $current_repeat_end =
      $first_repeat->[1] > $region_end ? $region_end : $first_repeat->[1];
    my %current_repeat_name = ( $first_repeat->[2] => 1 );
    my $name;

    foreach my $repeat ( @{$repeats} ) {
        if ( $repeat->[0] < $region_start ) {
            $repeat->[0] = $region_start;
        }
        if ( $repeat->[1] > $region_end ) {
            $repeat->[1] = $region_end;
        }
        if ( $repeat->[0] > $current_repeat_end + 1 ) {
            $name = join q{;}, sort keys %current_repeat_name;
            push @collapsed_repeats,
              [ $current_repeat_start, $current_repeat_end, $name, ];
            $current_repeat_start = $repeat->[0];
            $current_repeat_end   = $repeat->[1];
            %current_repeat_name  = ( $repeat->[2] => 1 );
        }
    }
    $name = join q{;}, sort keys %current_repeat_name;
    push @collapsed_repeats,
      [ $current_repeat_start, $current_repeat_end, $name, ];

    return \@collapsed_repeats;
}

# Count read pairs where one read overlaps a stranded region
sub count_overlapping {    ## no critic (ProhibitManyArgs)
    ## no critic (ProhibitReusedNames)
    my ( $hts, $region_chr, $region_start, $region_end, $region_strand,
        $total_count, $length )
      = @_;
    ## use critic

    confess sprintf 'Region %s:%d-%d:%d malformed', $region_chr, $region_start,
      $region_end, $region_strand
      if $region_start > $region_end;

    my $alignments = $hts->features(
        -seq_id   => $region_chr,
        -start    => $region_start,
        -end      => $region_end,
        -iterator => 1,
    );
    my %ignore;
    my $count = 0;
    while ( my $alignment = $alignments->next_seq ) {
        next if $alignment->get_tag_values('FLAGS') =~ m/\bDUPLICATE\b/xms;
        next if $alignment->qual < $MAPQ_THRESHOLD;
        next if $alignment->strand != $region_strand;
        next if $perfect_matches && @{ $alignment->cigar_array } > 1;
        next
          if exists $ignore{ $alignment->qname };  # Ignore mate of counted read
        if ( !$alignment->munmapped ) {
            $ignore{ $alignment->qname } = 1;
        }
        $count++;
    }

    return $count if !$total_count;

    # Also return FPKM
    if ( !$length ) {
        $length = $region_end - $region_start + 1;
    }
    return ( $count, $count / ( $total_count / 1e6 ) / $length * 1000 );
}

# Count read pairs where one read is entirely enclosed in a stranded region
sub count_enclosed {    ## no critic (ProhibitManyArgs)
    ## no critic (ProhibitReusedNames)
    my ( $hts, $region_chr, $region_start, $region_end, $region_strand,
        $total_count, $length )
      = @_;
    ## use critic

    confess sprintf 'Region %s:%d-%d:%d malformed', $region_chr, $region_start,
      $region_end, $region_strand
      if $region_start > $region_end;

    my $alignments = $hts->features(
        -seq_id   => $region_chr,
        -start    => $region_start,
        -end      => $region_end,
        -iterator => 1,
    );
    my %ignore;
    my $count = 0;
    while ( my $alignment = $alignments->next_seq ) {
        next if $alignment->start < $region_start;
        next if $alignment->end > $region_end;
        next if $alignment->get_tag_values('FLAGS') =~ m/\bDUPLICATE\b/xms;
        next if $alignment->qual < $MAPQ_THRESHOLD;
        next if $alignment->strand != $region_strand;
        next if $perfect_matches && @{ $alignment->cigar_array } > 1;
        next
          if exists $ignore{ $alignment->qname };  # Ignore mate of counted read

        if ( !$alignment->munmapped ) {
            $ignore{ $alignment->qname } = 1;
        }
        $count++;
    }

    return $count if !$total_count;

    # Also return FPKM
    if ( !$length ) {
        $length = $region_end - $region_start + 1;
    }
    return ( $count, $count / ( $total_count / 1e6 ) / $length * 1000 );
}

# Count read pairs where one read overlaps any of a number of stranded regions
sub count_overlapping_multiple {    ## no critic (ProhibitManyArgs)
    ## no critic (ProhibitReusedNames)
    my ( $hts, $region_chr, $region_pairs, $region_strand, $total_count,
        $length )
      = @_;
    ## use critic

    my %ignore;
    my $count = 0;
    foreach my $region_pair ( @{$region_pairs} ) {
        my ( $region_start, $region_end ) = @{$region_pair};
        confess sprintf 'Region %s:%d-%d:%d malformed', $region_chr,
          $region_start, $region_end, $region_strand
          if $region_start > $region_end;

        my $alignments = $hts->features(
            -seq_id   => $region_chr,
            -start    => $region_start,
            -end      => $region_end,
            -iterator => 1,
        );
        while ( my $alignment = $alignments->next_seq ) {
            next if $alignment->get_tag_values('FLAGS') =~ m/\bDUPLICATE\b/xms;
            next if $alignment->qual < $MAPQ_THRESHOLD;
            next if $alignment->strand != $region_strand;
            next if $perfect_matches && @{ $alignment->cigar_array } > 1;
            next
              if exists $ignore{ $alignment->qname };
            if ( !$alignment->munmapped ) {
                $ignore{ $alignment->qname } = 1;
            }
            $count++;
        }
    }

    return $count if !$total_count;

    # Also return FPKM
    if ( !$length ) {
        $length = 0;
        foreach my $region_pair ( @{$region_pairs} ) {
            my ( $region_start, $region_end ) = @{$region_pair};
            $length += $region_end - $region_start + 1;
        }
    }
    return ( $count, $count / ( $total_count / 1e6 ) / $length * 1000 );
}

# Count read pairs where one read is entirely enclosed in one of a number of
# stranded regions
sub count_enclosed_multiple {    ## no critic (ProhibitManyArgs)
    ## no critic (ProhibitReusedNames)
    my ( $hts, $region_chr, $region_pairs, $region_strand, $total_count,
        $length )
      = @_;
    ## use critic

    my %ignore;
    my $count = 0;
    foreach my $region_pair ( @{$region_pairs} ) {
        my ( $region_start, $region_end ) = @{$region_pair};
        confess sprintf 'Region %s:%d-%d:%d malformed', $region_chr,
          $region_start, $region_end, $region_strand
          if $region_start > $region_end;

        my $alignments = $hts->features(
            -seq_id   => $region_chr,
            -start    => $region_start,
            -end      => $region_end,
            -iterator => 1,
        );
        while ( my $alignment = $alignments->next_seq ) {
            next if $alignment->start < $region_start;
            next if $alignment->end > $region_end;
            next if $alignment->get_tag_values('FLAGS') =~ m/\bDUPLICATE\b/xms;
            next if $alignment->qual < $MAPQ_THRESHOLD;
            next if $alignment->strand != $region_strand;
            next if $perfect_matches && @{ $alignment->cigar_array } > 1;
            next
              if exists $ignore{ $alignment->qname };

            if ( !$alignment->munmapped ) {
                $ignore{ $alignment->qname } = 1;
            }
            $count++;
        }
    }

    return $count if !$total_count;

    # Also return FPKM
    if ( !$length ) {
        $length = 0;
        foreach my $region_pair ( @{$region_pairs} ) {
            my ( $region_start, $region_end ) = @{$region_pair};
            $length += $region_end - $region_start + 1;
        }
    }
    return ( $count, $count / ( $total_count / 1e6 ) / $length * 1000 );
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'bam_file=s'       => \$bam_file,
        'repeats_file=s'   => \$repeats_file,
        'perfect_matches'  => \$perfect_matches,
        'slice_regexp=s'   => \$slice_regexp,
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

    if ( !$bam_file ) {
        pod2usage("--bam_file must be specified\n");
    }
    if ( !$repeats_file ) {
        pod2usage("--repeats_file must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

analyse_intronic_expression.pl

Analyse intronic expression

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script analyses intronic expression.

=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-88/ensembl/modules \
        analyse_intronic_expression.pl \
            --bam_file sample.bam --repeats_file repeats.txt

=head1 USAGE

    analyse_intronic_expression.pl
        [--bam_file file]
        [--repeats_file file]
        [--perfect_matches]
        [--slice_regexp regexp]
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

=item B<--bam_file FILE>

BAM file to analyse.

=item B<--repeats_file FILE>

RepeatMasker repeats file.

=item B<--perfect_matches>

Only count perfectly matching reads (in terms of alignment).

=item B<--slice_regexp REGEXP>

Regular expression for limiting slices.

=item B<--species SPECIES>

Species (defaults to Mus musculus).

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

This software is Copyright (c) 2018 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
