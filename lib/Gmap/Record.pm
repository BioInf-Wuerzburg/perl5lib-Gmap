package Gmap::Record;

use warnings;
use strict;

# $Id$

# preference libs in same folder over @INC
use lib '../';

use Verbose;

use overload
	'.' => \&cat,
	'""' => \&string;


our $VERSION = '0.01';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;

##------------------------------------------------------------------------##

=head1 NAME 

Gmap::Record.pm

=head1 DESCRIPTION

Class for handling records of Gmap summary files.

=head1 SYNOPSIS

=cut

=head1 CHANGELOG

=head2 0.01

=over

=item [Init] Initial module. Provides Constructor and generic accessor 
 methods.

=back

=cut


=head1 TODO

=over

=item Synopsis

=back

=cut



##------------------------------------------------------------------------##

=head1 Class Attributes

=cut

=head2 $V

Verbose messages are handled using the Verbose.pm module. To 
 customize verbose message behaviour, overwrite the attribute with
 another Verbose object created with the Verbose module.

=cut

our $V = Verbose->new();

##------------------------------------------------------------------------##

=head1 Class METHODS

=cut

##------------------------------------------------------------------------##


=head1 Constructor METHOD

=head2 new

Create a new FASTA seq object. Either provide a single STRING 
 containing one FASTA record or a key => value construct, for which
 either C<seq_head> or C<id> is required.

  my $seq = ">seq1 blub\nATGC\n";
  Fasta::Seq->new($seq);
  # or
  Fasta::Seq->new(
  	id => 'seq1',
  	desc => 'blub',
  	seq => "ATGC",
  );
  # or
  Fasta::Seq->new(
  	seq_head => '>seq1 blub'
  	seq => "ATGC"
  );


=cut

sub new{
	my $proto = shift;
	my $self;
	my $class;
	# object method -> clone + overwrite
	if($class = ref $proto){
		$self = bless ({%$proto, @_}, $class)
	}else{ # init emtpy obj
		$class = $proto;
		$self = {
				id => '',
				desc => '',
				byte_offset => undef,
			};
	}
	
	if(@_){
		if(@_%2){ # input is string to split
			my %self;
			my ($paths, $alns);
			(@self{'id','desc'}, $paths, $alns) = shift =~ m/
				(?:>?(\S*))					# id, >? for records
				(?:[^\S\n]([^\n]+))?\n		# desc, optional
				(Paths.+)\n(Alignments.+)	# Paths and Alignments
			/xs;
								# . includes \n
			my $paths_desc;
			($paths_desc, $paths) = split(/\n/, $paths, 2);
			my @paths = split(/\n\n/, $paths);
			if(@paths){
				$paths[-1] =~ s/\n+$//;
				
				if($paths_desc =~ /chimera/){
					$self->{chimera}++;
					($self->{chimera_breakpoint}) = $paths_desc =~ /breakpoint at (\d+)/;
				}else{
					$self->{chimera} = 0;
					$self->{chimera_breakpoint} = undef;
				}
				
				for(my $i=0; $i<@paths; $i++){
					my ($path_desc, $body) = split(/\n/, $paths[$i], 2);
					#$self->{paths}[$i]{desc} = $path_desc;
					@{$self->{paths}[$i]}{qw(
						qry_hit_from
						qry_hit_to
						qry_hit_len
						ref_hit_from
						ref_hit_to
						ref_hit_strand
						ref_hit_len
					)} = $path_desc =~ m/(\d+)\.\.(\d+)\D+(\d+).+?([\d,]+)\.\.([\d,]+)\D+?(-?)(\d+)/; 
					
					$self->{paths}[$i]{ref_hit_from} =~ tr/,//d;
					$self->{paths}[$i]{ref_hit_to} =~ tr/,//d;
					if($self->{paths}[$i]{ref_hit_strand} eq '-'){
						$self->{paths}[$i]{ref_hit_rev} = 1;
					}else{
						$self->{paths}[$i]{ref_hit_rev} = 0;
						$self->{paths}[$i]{ref_hit_strand} = '+';
					}
					
					
					#$self->{paths}[$i]{body} = $body;
					my @body = split(/\n/, $body);
					chomp(@body);
					my %body = map{
						s/^\s+//; 
						my($k, $v) = split(/: /, $_, 2); 
						$k=~tr/ /_/; $k=lc($k); 
						($k, $v)
					}@body;
					
					# coverage
					($body{qry_cov}, $body{qry_len}) = $body{coverage} =~ m/(^[.\d]+)\D+(\d+)/;
					delete $body{coverage};
					# trimmed coverage
					($body{qry_trimmed_cov}, $body{qry_trimmed_len}, $body{qry_trimmed_from}, $body{qry_trimmed_to}) = $body{trimmed_coverage} =~ m/(^[.\d]+)\D+(\d+)\D+(\d+)\.\.(\d+)/;
					delete $body{trimmed_coverage};
					# percent_identity
					($body{aln_idy}, $body{aln_ma}, $body{aln_mm}, $body{aln_ins}, $body{aln_del}) = $body{percent_identity} =~ m/(^[.\d]+)\D+(\d+)\D+(\d+)\D+(\d+)\D+(\d+)/;
					delete $body{percent_identity};
					
					$self->{paths}[$i] = {%{$self->{paths}[$i]}, %body};
				}	
			}else{
				$self->{paths} = [];
			}
			
			# TODO
			#$self->{alignments} = $alns;
			
			$self = {
				%$self,
				%self,
				@_	# overwrite defaults
			};
		}else{
			$self = {
				%$self,
				@_	# overwrite defaults
			};
			
		}
		# make sure, id has not leading '>'
		$self->{id} =~ s/^>// if $self->{id};
		$self->{desc} = '' unless defined $self->{desc};
	}
	
	return bless $self, $class;
}







##------------------------------------------------------------------------##

=head1 Object METHODS

=cut


=head2 string([LINEWIDTH])

Get entire sequence as FASTA string. Provide optional line width.

=cut

sub string{
	my ($self, $lw) = @_;
}


##------------------------------------------------------------------------##

=head1 Accessor METHODS

=head2 seq_head

Get seq_head (>id desc).

=cut

sub seq_head{
	my ($self, $seq_head) = @_;
	return $self->{seq_desc} ? $self->{id}.' '.$self->{desc} : $self->{id};
}

=head2 byte_offset

Get/Set the byte_offset.

=cut

sub byte_offset{
	my ($self, $byte_offset) = @_;
	$self->{byte_offset} = $byte_offset if defined($byte_offset);
	return $self->{byte_offset};
}


=head2 id

Get/set the seqs id. Also updates seq_head.

=cut

sub id{
	my ($self, $id) = @_;
	if(defined $id){
		$self->{id} = $id;
	}
	return $self->{id};
}

=head2 desc

Get/set the seqs description, also updates seq_head.

=cut

sub desc{
	my ($self, $desc) = @_;
	if(defined $desc){
		$self->{desc} = $desc;
	}
	return $self->{desc};
}


=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



