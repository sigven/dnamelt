#!/usr/bin/perl -w

use strict;
use DNAMelting::Meltingprofiles;

no warnings;

sub get_temperature_profile{
	 my $seq = shift;	 
	 
	 my $helicity = 0.5;
	 my $algorithm = "Approx";
	 my $nat_cons = 0.0156;
	 my $thermodyn = "Blossey and Carlon 2003";
	 
	 my @temps = &DNAMelting::Meltingprofiles::calc_tm_prof($helicity,$algorithm,$thermodyn,$seq,$nat_cons);
	#print length($seq) . "=" . scalar(@temps) . '\n';
	#print
	if($temps[0] eq ''){
             shift @temps;
          } 

	 return @temps;
}

print join(' ', get_temperature_profile($ARGV[0]));


