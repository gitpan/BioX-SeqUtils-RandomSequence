use Test::More tests => 4;

BEGIN {
use_ok( 'BioX::SeqUtils::RandomSequence' );
}

my $randomizer = BioX::SeqUtils::RandomSequence->new();
my $test = $randomizer->rand_nuc(), "\n";
ok( $test, "random nucleotide");

my $test = $randomizer->rand_pro(), "\n";
ok( $test, "random protein");

my $test = $randomizer->rand_pro_set(), "\n";
ok( $test, "random protein set");

diag( "Testing BioX::SeqUtils::RandomSequence $BioX::SeqUtils::RandomSequence::VERSION" );
