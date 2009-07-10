#!/usr/bin/perl -w

use Getopt::Std;
use BioX::SeqUtils::RandomSequence;

getopts('l:a:c:t:g:');  
check_opts();

my $randomizer = BioX::SeqUtils::RandomSequence->new({ l => $opt_l, 
                                                       a => $opt_a,
                                                       c => $opt_c,
                                                       g => $opt_g,
                                                       t => $opt_t });
print $randomizer->rand_pro(), "\n";

exit;

sub usage { print "   USAGE: ./randomize.pp -l<length> -a<rate> -c<rate> -g<rate> -t<rate>\n\n"; exit; }

sub defaults { $opt_a = $opt_c = $opt_g = $opt_t = 1; }

sub check_opts {
	if ( !$opt_a or !$opt_c or !$opt_g or !$opt_t ) { defaults(); }
	if ( !$opt_l ) { $opt_l = $opt_a + $opt_c + $opt_g + $opt_t; }
}

