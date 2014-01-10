#!/usr/bin/env perl

# $Id: AnalyseGmapSry.pl 55 2013-05-15 11:41:39Z s187512 $

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;
use Log::Log4perl qw(:easy :no_extra_logdie_message);

use FindBin qw($RealBin);
use List::Util;



use lib "$RealBin/../lib/";
use Gmap::Record;
use Gmap::Parser;

our $VERSION = '0.01';

Log::Log4perl->init(\<<'CFG');
	log4perl.logger.main				= DEBUG, Screen
	log4perl.appender.Screen			= Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr		= 0
	log4perl.appender.Screen.layout		= PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [blat.pl] %m%n

CFG

my $L = Log::Log4perl::get_logger();
$L->level($INFO);

=head1 NAME

AnalyseGmapSry - Parse Gmap summary file and create stats.

=cut

=head1 CHANGELOG

=head2 0.02

=over

=item [Change] Log::Log4perl replaces Verbose module.

=item [Change] Condensed Getopt reading, values to hash only.

=back

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

=over

=item [--in=<FASTA/FASTQ>]

Input FASTA/FASTQ file. Default STDIN.

=item [--out=<STRING>]

Output file. Default off. Specify '-' for STDOUT.

=item [--ref-cov=<-/FILE>]

Calculate and output per base coverage of reference sequences. Only 
 sequences with at least one hit are reported. Default off. Specify 
 '-' for STDOUT, or a FILE to write to.

=item [--allow-edge-mappings] [OFF]

When dealing with fragmented or circular references, reads might be 
 partially mapped to ends of reference sequences. This introduces
 false errors since overhangs of these mappings will be treated as
 dropped and might also look like chimeras. To prevent these artefacts
 mappings over reference ends are by default ignored.

=item --quiet | --debug

Decrease/increase verbose level.

=item [--help]

Display this help

=back

=cut

my %opt = (
);

Getopt::Long::Configure("no_ignore_case");
GetOptions(\%opt, qw(
	in=s
	out=s
	ref_cov|ref-cov=s
	allow_edge_mappings|allow-edge-mappings!
	quiet
	debug
	help|h
)) or pod2usage(1);

pod2usage(1) if $opt{help};


##------------------------------------------------------------------------##

=head1 MAIN

=cut

=head2 globals

=cut

my %S = (
	p0 => [], # 0 == chimeric
);
my ($Unmapped, $Record_count, $Chimera_count, $Edge_mapped) = (0,0,0);
my %REF_COV;

##------------------------------------------------------------------------##


=head2 check input file

=cut

$L->info("Reading ".( $opt{in} ? $opt{in} : "<&STDIN" ));

my $gp = Gmap::Parser->new(file => $opt{in})->check_format() 
	|| $L->exit("Data does not look like gmap summary"); # defaults to &STDIN


=head2 loop

=cut

while(my $r = $gp->next_record){
	$Record_count++;
	my $pc = @{$r->{paths}};
	unless($pc){
		$Unmapped++;
	}else{
		if ($r->{chimera}){
			$pc = 0;
		};
		
		# analyse primary path: path "0""
		my $p0 = $r->{paths}[0];
		
		## exclude
		# edge mappings
		unless($opt{allow_edge_mappings}){
			# check if mapping is to close to an end
			if(
				$p0->{qry_len} > $p0->{ref_hit_to}
				|| 
				$p0->{qry_len} > $p0->{ref_len} - $p0->{ref_hit_from}
			){
				$Edge_mapped++;
				next; 
			}
		}
		# N's (ref/qry)
		
		$S{p0}[$pc]{ma} +=$p0->{aln_ma};
		$S{p0}[$pc]{mm} +=$p0->{aln_mm};
		$S{p0}[$pc]{in} +=$p0->{aln_ins};
		$S{p0}[$pc]{de} +=$p0->{aln_del};

#		if($opt{ref_cov}){
#			# init ref cov
#			unless(exists $REF_COV{$p0->{ref_id}}){
#				$REF_COV{$p0->{ref_id}} = [(0)x $p0->{ref_len}];
#			}
#			map{$_++}@{$REF_COV{$p0->{ref_id}}}[$p0->{ref_hit_from}..$p0->{ref_hit_to}];
#		}
		if($r->{chimera}){
			my $p1 = $r->{paths}[1];
			$Chimera_count++;
			$S{p0}[$pc]{ma} +=$p1->{aln_ma};
			$S{p0}[$pc]{mm} +=$p1->{aln_mm};
			$S{p0}[$pc]{in} +=$p1->{aln_ins};
			$S{p0}[$pc]{de} +=$p1->{aln_del};

			$S{p0}[$pc]{dr} += ($p0->{qry_len}-$p0->{qry_hit_len}-$p1->{qry_hit_len});
			
#			if($opt{ref_cov}){
#				# init ref cov
#				unless(exists $REF_COV{$p1->{ref_id}}){
#					$REF_COV{$p1->{ref_id}} = [(0)x $p1->{ref_len}];
#				}
#				map{$_++}@{$REF_COV{$p1->{ref_id}}}[$p1->{ref_hit_from}..$p1->{ref_hit_to}];
#			}
		}else{
			$S{p0}[$pc]{dr} += ($p0->{qry_len}-$p0->{qry_hit_len});
		}
	}
}

#print Dumper($S{p0});

foreach my $i (1..@{$S{p0}}-1, 0){
	my $p = $i || "Chimera";
	print Dumper($p, $S{p0}[$i]);
}

my $pat ="%10s: %10s %8s %8.2f%%\n";

#my $TOT = $MM+$MA+$INS+$DEL+$Dropped;
#printf "%10s: %10s %8s %8s\n", '', qw(bp reads resp_%);
#printf $pat, 'Matches', $MA, '', $MA/$TOT*100;
#printf $pat, 'Mismatches', $MM, '', $MM/$TOT*100;
#printf $pat, 'Inserts', $INS, '', $INS/$TOT*100;
#printf $pat, 'Deletions', $DEL, '', $DEL/$TOT*100;
#printf $pat, 'Dropped', $Dropped, '', $Dropped/$TOT*100;

printf $pat, 'Unmapped', '', $Unmapped, $Unmapped/$Record_count*100;
printf $pat, 'Chimeras', '', $Chimera_count, $Chimera_count/$Record_count*100;
unless($opt{allow_edge_mappings}){
	printf $pat, 'Edgemapped', '', $Edge_mapped, $Edge_mapped/$Record_count*100;
}


if($opt{ref_cov}){
	my $rcfh;
	unless($opt{ref_cov} eq '-'){
		open($rcfh, '>', $opt{ref_cov}) or die "$!: $opt{ref_cov}";
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
