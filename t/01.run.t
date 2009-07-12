use Test::More tests => 5;

BEGIN {
use_ok( 'BioX::SeqUtils::RandomSequence' );
}

my $randomizer = BioX::SeqUtils::RandomSequence->new();
my $test = $randomizer->rand_nuc(), "\n";
ok( $test, "random nucleotide");

my $test = $randomizer->rand_pro(), "\n";
ok( $test, "random protein");

my $test = $randomizer->rand_pro_set(), "\n";
ok( $test, "random protein set (scalar)");

my ($test, $dummy) = $randomizer->rand_pro_set(), "\n";
ok( $test, "random protein set (list)");

diag( "Testing BioX::SeqUtils::RandomSequence $BioX::SeqUtils::RandomSequence::VERSION" );
