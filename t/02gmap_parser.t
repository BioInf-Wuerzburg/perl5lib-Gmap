#!/usr/bin/env perl

# $Id$

use strict;
use warnings;

use Test::More;
use Data::Dumper;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";

BEGIN { use_ok('Gmap::Parser'); }

# globals
my $Class = 'Gmap::Parser';
my ($gp, $gp_string);


#--------------------------------------------------------------------------#

=head2 sample data

=cut

my $Dat_file = "01gmap_record.dat"; # data
my $Dmp_file = "01gmap_record.dmp"; # data structure dumped

my $Dat = do { local $/; local @ARGV = $Dat_file; <> }; # slurp data to string
my @Dat = split(/(?<=\n)(?=>)/, $Dat);

my @Dmp = do "$Dmp_file"; # read and eval the dumped structure


#--------------------------------------------------------------------------#

=head2 new

=cut

# NOTE: cannot use is_deeply as it crashes on overloaded objects...
# string

subtest 'new and DESTROY' => sub{
	$gp = new_ok($Class);
	
	#  default attributes
	is($gp->{mode}, '<', 'Default attribute "mode"');
	is($gp->{file}, undef, 'Default attribute "file"');
	ok(fileno $gp->{fh} > 2, 'Default attribute "fh"'); # copy of STDIN needs to be greater than 3
	
	# DESTROY closes fh
	my $fh = $gp->{fh};
	$gp = undef;
	is(fileno($fh), undef, 'Autoclose filehandle');
	# read from file
	$gp = new_ok($Class, [
		mode => '<', 
		file => $Dat_file,
	]);
	#  default attributes
	is($gp->{mode}, '<', 'Custom attribute "mode"');
	is($gp->{file}, $Dat_file, 'Custom attribute "file"');
};


subtest 'filehandle' => sub{
	# fh
	can_ok($gp, 'fh');
	can_ok($gp, 'is_fh');
	my $gph = $gp->{fh};
	is($gp->fh, $gph, 'Get $obj->fh()');
	is($gp->is_fh, 0, '$obj->is_fh() FILE');
	ok($gp->is_fh('FILE'), '$obj->is_fh("FILE")');
	
	is($gp->fh(\*STDIN), \*STDIN, 'Set $obj->fh()');
	is($gp->is_fh, 1, '$obj->is_fh() PIPE');
	ok($gp->is_fh('PIPE'), '$obj->is_fh("PIPE")');
	
	is($gp->fh($gph), $gph , 'Reset $obj->fh() FILE');
	is($gp->is_fh, 0, '$obj->is_fh()');

	$gp_string = Gmap::Parser->new(file => \$Dat);
	is($gp_string->is_fh(), 2, '$obj->is_fh() SCALAR');	
	ok($gp_string->is_fh("SCALAR"),'$obj->is_fh("SCALAR")');	
	
};

# next_record
subtest 'next_record' => sub{
	can_ok($gp, 'next_record');
	for(my $i=0; $i<@Dat; $i++){
		is($gp->next_record->{id}, $Dmp[$i]->{id}, 'Get $obj->next_record()');	
	}
	is($gp->next_record, undef, 'Get $obj->next_record() eof');
	is($gp->next_record->{id}, $Dmp[0]->{id}, 'Get $obj->next_record() restart after eof');
};


# check_format
subtest 'check_format' => sub{
	can_ok($gp, 'check_format');
	is(Gmap::Parser->new(file => \("Not a gmap record"))->check_format, 
		undef, '$obj->check_format on not Gmap is undef');
	isa_ok(Gmap::Parser->new(file => $Dat_file)->check_format, 
		$Class, '$obj->check_format on Gmap');
};

done_testing();
