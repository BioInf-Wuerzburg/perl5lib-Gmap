package Gmap::Parser;

# $Id: Parser.pm 55 2013-05-15 11:41:39Z s187512 $

use warnings;
use strict;

use IO::File;

# preference libs in same folder over @INC
use lib '../';

use Gmap::Record 0.01;

use List::Util;

our $VERSION = '0.01';
our ($REVISION) = '$Revision: 55 $' =~ /(\d+)/;
our ($MODIFIED) = '$Date: 2013-05-15 13:41:39 +0200 (Wed, 15 May 2013) $' =~ /Date: (\S+\s\S+)/;



##------------------------------------------------------------------------##

=head1 NAME 

Gmap::Parser.pm

=head1 DESCRIPTION

Parser module to read Gmap summary files.

=head1 SYNOPSIS

  my $gp = Gmap::Parser->new(file => $opt_in)->check_format() 
    || $V->exit("Data does not look like gmap summary"); # defaults to &STDIN
  
  while(my $r = $gp->next_record){
    # do something with the record
  }

=cut

=head1 CHANGELOG

=head2 0.01

=over

=item [Feature] << $fp->check_format >> now reads and B<unreads> the first
 char form input to determine format. Unreading makes it safe to use on 
 STDIN without using up stuff from the stream.

=item [Init] Initial Parser module. Provides Constructor and generic accessor
 methods

=back

=cut

=head1 TODO

=over

=back

=cut




##------------------------------------------------------------------------##

=head1 Constructor METHOD

=head2 new

Initialize a parser object. Takes parameters in key => value format. 
 Parameters are:

  fh => undef,
  file => undef,
  mode => '<',   # read, 
                 # '+>': read+write (clobber file first)
                 # '+<': read+write (append)
                 # '>' : write (clobber file first)
                 # '>>': write (append)
  # NOTE: writing modes are considered experimental

=cut

sub new{
	my $class = shift;
	
	# defaults=
	my $self = {
		fh => undef,
		file => undef,
		mode => '<',
		_is_fh => undef, # 0 => FILE, 1 => PIPE, 2 => SCALAR
		@_	# overwrite defaults
	};

	if($self->{file} && $self->{fh}){
		die sprintf("%s: %s",(caller 0)[3],"Can only take either file or fh!");
	}

	bless $self, $class;

	# open file in read/write mode
	if($self->{file}){
		my $fh;
		open ( $fh , $self->{mode}, $self->{file}) or die sprintf("%s: %s, %s",(caller 0)[3],$self->{file}, $!);
		$self->{fh} = $fh;
		if(ref $self->{file} eq 'SCALAR'){
			$self->{_is_fh} = 2;
		}elsif(-f $self->{file}){
			$self->{_is_fh} = 0;
		}else{
			die sprintf("%s: %s",(caller 0)[3],"file neither plain file nor SCALAR reference");
		}
	}else{
		if($self->{fh}){
			$self->fh($self->{fh});
		}else{
			open(my $fh, "<&STDIN") or die $!;
			$self->{fh} = $fh;
			$self->{_is_fh} = 1;
		}
	}
	
	return $self;
}

sub DESTROY{
	# just to be sure :D
	my $self = shift;
	close $self->fh;
}



##------------------------------------------------------------------------##

=head1 Object METHODS

=cut

=head2 next_record

Loop through file and return next 'Gmap::Record' object.

=cut


sub next_record{
	my ($self) = @_;
	
	my $fh = $self->{fh};
	local $/ = "\n>";
	
	# return fasta seq object
	my $byte_offset = tell($fh);
	my $r = <$fh>;
	unless(defined $r){
		seek($fh,0,0); # reset to file start
		return 
	};
	chomp($r);
	return Gmap::Record->new($r, byte_offset => $byte_offset);
}




=head2 check_format

Takes a peek at the first entry in the file and checks wether the format of 
 the input looks like FASTQ (leading @). Returns the Parser object on 
 success, undef on failure. Does not modify the input stream, therefore
 can be used on STDIN safely.

NOTE: It only works at the start of the input. This means for pipes, use it
 before you start reading, on files use it either in the beginning or seek
 to the start to perform the check.

=cut

sub check_format{
	my ($self) = @_;
	my $fh = $self->fh;
	die sprintf("%s: %s",(caller 0)[3],"Format checking only works at the start of the file") 
		if $fh->tell;
	my $c =$fh->getc(); # read first char
	$fh->ungetc(ord($c)); # unread first char
	# read first char
	return $c eq '>' ? $self : undef;
}


=head2 seek

Set the filehandle to the specified byte offset. Takes two
optional arguments "POSITION" (0), "WHENCE" (0), see perl "seek" for more.
Returns 'true' on success.

NOTE: this operation does only work on real files, not on STDIN.

=cut

sub seek{
	my ($self, $offset, $whence) = (@_, 0, 0);
	return seek($self->fh, $offset, $whence);
}


=head2 is_fh

Determine the type of the filehandle. Without parameter, returns 0 for 
 handle to FILE, 1 for PIPE, and 2 for a handle to a SCALAR.
Alternatively you can provide the name of the type or the corresponding
 INT as single parameter. In these cases, the methods returns 1, if the
 type is matched, 0 otherwise.

  $fp->is_fh();  # 0,1 or 2
  $fp->is_fh('FILE') # 1 if filehandle is a handle to a FILE, 0 otherwise
  $fp->is_fh(0) # 1 if filehandle is a handle to a FILE, 0 otherwise

=cut

sub is_fh{
	my ($self, $type) = @_;
	
	my %type = (
		'FILE' => 0,
		'PIPE' => 1,
		'SCALAR' => 2,
		0 => 0,
		1 => 1,
		2 => 2,
	);
	
	if(defined $type){
		die sprintf("%s: %s",(caller 0)[3],"unknown type $type") 
			unless exists $type{$type};
		
		return $type{$type} == $self->{_is_fh} ? 1 : 0;
		
	}

	return $self->{_is_fh};
}


=head2 append_record

NOT YET WORKING, Gmap::Record is missing string method

Append an sequence to the file, provided as object or string. Returns the
 byte offset position in the file.

NOTE: In case a string is provided, make sure it contains trailing newline 
 since no further test is performed.

=cut

sub append_record{
	my ($self, $seq) = @_;
	my $pos = tell($self->{fh});
	print {$self->{fh}} "$seq";
	return $pos;
}


=head2 append_tell

NOT YET WORKING, Gmap::Record is missing string method

Return the byte offset of the current append filehandle position

=cut

sub append_tell{
	my ($self) = @_;
	return tell($self->{fh});
}



##------------------------------------------------------------------------##

=head1 Accessor METHODS

=cut

=head2 fh

Get/Set the file handle. Determines the source read from - FILE, PIPE, STRING.

=cut

sub fh{
	my ($self, $fh) = @_;
	
	if($fh){
		if(-f $fh){
			$self->{_is_fh} = 0;
		}elsif(-p $fh or  -t $fh){
			$self->{_is_fh} = 1;
		}else{
			$self->{_is_fh} = 2;
		}
		$self->{fh} = $fh 
	}
	
	return $self->{fh};
}




##------------------------------------------------------------------------##

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



