#!/usr/bin/env perl

# $Id: AnalyseGmapSry.pl 55 2013-05-15 11:41:39Z s187512 $

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";

use List::Util;

use Verbose;
use Verbose::ProgressBar;

use Gmap::Record;
use Gmap::Parser;

our $VERSION = '0.01';

=head1 NAME

AnalyseGmapSry - Parse Gmap summary file and create stats.

=cut

=head1 CHANGELOG

=head2 0.01

=over

=item [Init]

=back

=cut

=head1 TODO

=over

=back

=cut





##------------------------------------------------------------------------##

=head1 SYNOPSIS

  perl AnalyseGmapSry.pl --in <SRYFILE> <OPTIONS>

=cut

=head1 OPTIONS

=cut

=over

=cut

my %opt;

=item [--in=<FASTA/FASTQ>]

Input FASTA/FASTQ file. Default STDIN.

=cut

$opt{'in=s'} = \(my $opt_in = undef);

=item [--out=<STRING>]

Output file. Default off. Specify '-' for STDOUT.

=cut

$opt{'out=s'} = \(my $opt_out = undef);

=item [--ref-cov=<-/FILE>]

Calculate and output per base coverage of reference sequences. Only 
 sequences with at least one hit are reported. Default off. Specify 
 '-' for STDOUT, or a FILE to write to.

=cut

$opt{'ref-cov=s'} = \(my $opt_ref_cov = undef);

=item [--allow-edge-mappings] [OFF]

When dealing with fragmented or circular references, reads might be 
 partially mapped to ends of reference sequences. This introduces
 false errors since overhangs of these mappings will be treated as
 dropped and might also look like chimeras. To prevent these artefacts
 mappings over reference ends are by default ignored.

=cut

$opt{'allow-edge-mappings!'} = \(my $opt_edge_mappings = undef);

=item --verbose=<INT>

Toggle verbose level, default 2, which outputs statistics and progress.
Set 1 for statistics only, 0 for no verbose output and 3 for debug level
messages.

=cut 

$opt{'verbose=i'} = \(my $opt_verbose = 3);

=item [--help]

Display this help

=cut

$opt{'help|?'} = \(my $opt_help);

=back

=cut

GetOptions(%opt) or pod2usage(1);
pod2usage(1) if $opt_help;

if($opt_verbose > 2){
	use Data::Dumper;
	$Data::Dumper::Sortkeys = 1;
} 

##------------------------------------------------------------------------##

=head1 MAIN

=cut

=head2 globals

=cut

my $V = Verbose->new(
	level => 1,
	report_level => $opt_verbose,
	line_width => 80,
	
);

my ($MA, $MM, $INS, $DEL, $Dropped,$Unmapped, $Record_count, $Chimera_count) = (0,0,0,0,0,0,0,0);
my %REF_COV;

##------------------------------------------------------------------------##


=head2 check input file

=cut

$V->verbose("Reading ".( $opt_in ? $opt_in : "<&STDIN" ));

my $gp = Gmap::Parser->new(file => $opt_in)->check_format() 
	|| $V->exit("Data does not look like gmap summary"); # defaults to &STDIN


=head2 loop

=cut

while(my $r = $gp->next_record){
	$Record_count++;
	unless(@{$r->{paths}}){
		$Unmapped++;
	}else{
		my $p0 = $r->{paths}[0];
		
		unless($opt_edge_mappings){
			# check if mapping is to close to an end
			next if $p0->{qry_len} > $p0->{ref_hit_to};
			next if $p0->{qry_len} > $p0->{ref_len} - $p0->{ref_hit_from};
		}
		
		$MA+=$p0->{aln_ma};
		$MM+=$p0->{aln_mm};
		$INS+=$p0->{aln_ins};
		$DEL+=$p0->{aln_del};

		if($opt_ref_cov){
			# init ref cov
			unless(exists $REF_COV{$p0->{ref_id}}){
				$REF_COV{$p0->{ref_id}} = [(0)x $p0->{ref_len}];
			}
			map{$_++}@{$REF_COV{$p0->{ref_id}}}[$p0->{ref_hit_from}..$p0->{ref_hit_to}];
		}
		if($r->{chimera}){
			my $p1 = $r->{paths}[1];
			$Chimera_count++;
			$MA+=$p1->{aln_ma};
			$MM+=$p1->{aln_mm};
			$INS+=$p1->{aln_ins};
			$DEL+=$p1->{aln_del};

			$Dropped+=($p0->{qry_len}-$p0->{qry_hit_len}-$p1->{qry_hit_len});
			
			if($opt_ref_cov){
				# init ref cov
				unless(exists $REF_COV{$p1->{ref_id}}){
					$REF_COV{$p1->{ref_id}} = [(0)x $p1->{ref_len}];
				}
				map{$_++}@{$REF_COV{$p1->{ref_id}}}[$p1->{ref_hit_from}..$p1->{ref_hit_to}];
			}
		}else{
			$Dropped+=($p0->{qry_len}-$p0->{qry_hit_len});
		}
	}
}

my $TOT = $MM+$MA+$INS+$DEL+$Dropped;
my $pat="%10s: %10s %8s %8.2f%%\n";
printf "%10s: %10s %8s %8s\n", '', qw(bp reads resp_%);
printf $pat, 'Matches', $MA, '', $MA/$TOT*100;
printf $pat, 'Mismatches', $MM, '', $MM/$TOT*100;
printf $pat, 'Inserts', $INS, '', $INS/$TOT*100;
printf $pat, 'Deletions', $DEL, '', $DEL/$TOT*100;
printf $pat, 'Dropped', $Dropped, '', $Dropped/$TOT*100;

printf $pat, 'Unmapped', '', $Unmapped, $Unmapped/$Record_count*100;
printf $pat, 'Chimeras', '', $Chimera_count, $Chimera_count/$Record_count*100;


if($opt_ref_cov){
	my $rcfh;
	unless($opt_ref_cov eq '-'){
		open($rcfh, '>', $opt_ref_cov) or die "$!: $opt_ref_cov";
	}
	select $rcfh if $rcfh;
	print "count\tcov\tref\n";
	
	foreach (sort keys %REF_COV){
		my %cov;
		map{$cov{$_}++}@{$REF_COV{$_}};
		my $max = List::Util::max(keys %cov);
		for(my $i=0; $i<=$max;$i++){
			print exists $cov{$i} ? $cov{$i} : 0, "\t$i\t$_\n";
		}
	}
	
	close $rcfh if $rcfh;
	select STDOUT;
}


=pod
>m111006_202713_42141_c100202382555500000315044810141104_s1_p0/8043/4690_6133 SUBSTR:26,1265 SUBSTR:11,1239
    Coverage: 48.9 (query length: 1239 bp)
    Coverage: 51.1 (query length: 1239 bp)
>m111006_202713_42141_c100202382555500000315044810141104_s1_p0/8075/6853_7408 SUBSTR:64,505 SUBSTR:0,478
    Coverage: 99.0 (query length: 478 bp)
    Coverage: 95.4 (query length: 478 bp)
>m111006_202713_42141_c100202382555500000315044810141104_s1_p0/8043/1822_3120 SUBSTR:21,1263 SUBSTR:24,1228
    Coverage: 48.5 (query length: 1228 bp)
    Coverage: 51.5 (query length: 1228 bp)

=cut

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut
