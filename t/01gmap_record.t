#!/usr/bin/env perl

# $Id$

use strict;
use warnings;

use Test::More;
use Test::Deep;
use Data::Dumper;

use File::Basename;
use FindBin qw($RealBin);
use lib "$RealBin/../lib/";

#--------------------------------------------------------------------------#

=head2 sample data

=cut

#--------------------------------------------------------------------------#

=head2 load module

=cut

BEGIN { use_ok('Gmap::Record'); }

my $Class = 'Gmap::Record';

#--------------------------------------------------------------------------#

=head2 sample data

=cut

(my $Dat_file = $FindBin::RealScript) =~ s/t$/dat/; # data
(my $Dmp_file = $FindBin::RealScript) =~ s/t$/dmp/; # data structure dumped

my $Dat = do { local $/; local @ARGV = $Dat_file; <> }; # slurp data to string
my @Dat = split(/(?<=\n)(?=>)/, $Dat);
my @Dmp = do "$Dmp_file"; # read and eval the dumped structure

#--------------------------------------------------------------------------#

=head2 new

=cut

# NOTE: cannot use is_deeply as it crashes on overloaded objects...
# string
subtest 'new' => sub {
	for(my $i=0; $i< @Dat; $i++){
		my $set = "set".($i+1);
		# from string
		my $rec = new_ok($Class, [$Dat[$i]], "$set from string");
		my $rec_cp = {%$rec}; #  copy obj to hash to make is_deeply work with overloaded obj..
		is_deeply( $rec_cp, $Dmp[$i], "$set from string attributes" );
		# from hash
		my $rec_hash = new_ok($Class, [%{$Dmp[$i]}], "$set from hash");
		$rec_cp = {%$rec_hash}; #  copy obj to hash to make is_deeply work with overloaded obj..
		is_deeply( $rec_cp, $Dmp[$i], "$set from hash attributes" );
		# from clone
		my $rec_clone = $rec->new();
		isa_ok($rec_clone, $Class, "$set cloned");
		$rec_cp = {%$rec_clone}; #  copy obj to hash to make is_deeply work with overloaded obj..
		is_deeply( $rec_cp, $Dmp[$i], "$set cloned attributes " );
	}
};


done_testing();

__END__


