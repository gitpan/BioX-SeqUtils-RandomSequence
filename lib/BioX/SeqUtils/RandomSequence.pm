package BioX::SeqUtils::RandomSequence;
use Class::Std;
use Class::Std::Utils;
use Bio::Tools::CodonTable;

use warnings;
use strict;
use Carp;

use version; our $VERSION = qv('0.9.3');

{
        my %type_of      :ATTR( :get<type>     :set<type>     :default<'n'>    :init_arg<y> );
        my %length_of    :ATTR( :get<length>   :set<length>   :default<'4'>    :init_arg<l> );
        my %table_of     :ATTR( :get<table>    :set<table>    :default<'1'>    :init_arg<s> );
        my %a_freq_of    :ATTR( :get<a_freq>   :set<a_freq>   :default<'1'>    :init_arg<a> );
        my %c_freq_of    :ATTR( :get<c_freq>   :set<c_freq>   :default<'1'>    :init_arg<c> );
        my %g_freq_of    :ATTR( :get<g_freq>   :set<g_freq>   :default<'1'>    :init_arg<g> );
        my %t_freq_of    :ATTR( :get<t_freq>   :set<t_freq>   :default<'1'>    :init_arg<t> );
        my %tmpl_of      :ATTR( :get<tmpl>     :set<tmpl>     :default<''>     );
                
        sub START {
                my ($self, $ident, $arg_ref) = @_;
		$self->_check_type();
		$self->_reset_tmpl();
                return;
        }

	sub rand_seq {
		my ( $self, $arg_ref ) = @_;

		# Set parameters redefined by this method
		$self->_args_to_attributes( $arg_ref );

		my $type = $self->get_type();
		if    ( $type =~ m/^n/ ) { return $self->rand_nuc(); }
		elsif ( $type =~ m/^d/ ) { return $self->rand_nuc({ l => 2}); }
		elsif ( $type =~ m/^p/ ) { return $self->rand_pro(); }
		elsif ( $type =~ m/^s/ ) { return $self->rand_pro_set(); }
                return;
        }

	sub rand_nuc {
		my ( $self, $arg_ref ) = @_;

		# Set parameters redefined by this method
		$self->_args_to_attributes( $arg_ref );

		# Create random nucleotide sequence of specified length
		my $tmpl_length = $self->get_length();
		my $nucleotide  = $self->randomize_tmpl();
		while ( length($nucleotide) < $tmpl_length ) { $nucleotide .= $self->randomize_tmpl(); }
		$nucleotide     =~ s/^([ACGT]{$tmpl_length}).*$/$1/;
                return $nucleotide;
        }

	sub rand_pro {
		my ( $self, $arg_ref ) = @_;

		# Set parameters redefined by this method
		$self->_args_to_attributes( $arg_ref );

		my $opt_l       = $arg_ref->{l} ? $arg_ref->{l} * 3 : $self->get_length() * 3;
		my $seq         = $self->rand_nuc({ l => $opt_l });
		my $codon_table = Bio::Tools::CodonTable->new( -id => $self->get_table() );
		if ( $codon_table->name() eq '' ) { print "  Error: Codon Table " . $self->get_table() . " not defined.\n"; exit; }
		my $protein     = $codon_table->translate( $seq );
		
                return $protein;
        }

	sub rand_pro_set {
		my ( $self, $arg_ref ) = @_;

		# Set parameters redefined by this method
		$self->_args_to_attributes( $arg_ref );

		my $opt_l       = $arg_ref->{l} ? $arg_ref->{l} * 3 : $self->get_length() * 3;
		my $seq         = $self->rand_nuc({ l => $opt_l + 1 });
		my $seq1        = $seq; $seq1  =~ s/.$//;  # Remove the last base
		my $seq2        = $seq; $seq2  =~ s/^.//;  # Remove the first base
		
		my $codon_table = Bio::Tools::CodonTable->new( -id => $self->get_table() );
		if ( $codon_table->name() eq '' ) { print "  Error: Codon Table " . $self->get_table() . " not defined.\n"; exit; }
		my $protein1    = $codon_table->translate( $seq1 );
		my $protein2    = $codon_table->translate( $seq2 );
		
                return wantarray() ? ( $protein1, $protein2 ) : [ $protein1, $protein2 ];
        }

	sub randomize_tmpl {
		my ( $self ) = @_;
		my $tmpl     = $self->get_tmpl();
		my @tmpl     = @$tmpl;
		for ( my $i = @tmpl; $i >= 0; --$i ) {
			my $j = int rand ($i + 1);
			next if $i == $j;
			@tmpl[$i,$j] = @tmpl[$j,$i];
		}
		no warnings 'all';
		return join("", @tmpl);
	}

	sub _args_to_attributes {
		my ( $self, $arg_ref ) = @_;

		if ( defined $arg_ref->{y} ) { $self->set_type( $arg_ref->{y} ); }
		if ( defined $arg_ref->{l} ) { $self->set_length( $arg_ref->{l} ); }

		my $freq_changed       = 0;

		if ( defined $arg_ref->{a} ) { $self->set_a_freq( $arg_ref->{a} ); $freq_changed++; }
		if ( defined $arg_ref->{c} ) { $self->set_c_freq( $arg_ref->{c} ); $freq_changed++; }
		if ( defined $arg_ref->{g} ) { $self->set_g_freq( $arg_ref->{g} ); $freq_changed++; }
		if ( defined $arg_ref->{t} ) { $self->set_t_freq( $arg_ref->{t} ); $freq_changed++; }
		
		# All frequencies must be set together or they are all reset to 1
		if ( $freq_changed && $freq_changed < 4 ) { $self->set_a_freq( 1 ); $self->set_c_freq( 1 ); $self->set_g_freq( 1 ); $self->set_t_freq( 1 ); }

		$self->_check_type();
		$self->_reset_tmpl() if $freq_changed;
		return;
	}

	sub _check_type {
		my ( $self, $arg_ref ) = @_;
		my $type = $self->get_type(); 
		if ( $type =~ m/^[^ndps]/ )   { print STDERR "  Error: Type (y) must be n, d, p, or s.\n"; exit; }
		return;
	}

	sub _reset_tmpl {
		my ( $self, $arg_ref ) = @_;
		my @tmpl = split( //, 'A' x $self->get_a_freq() . 
		                      'C' x $self->get_c_freq() . 
				      'G' x $self->get_g_freq() . 
				      'T' x $self->get_t_freq() );
		$self->set_tmpl( \@tmpl );
		return;
	}
}

1; # Magic true value required at end of module
__END__

=head1 NAME

BioX::SeqUtils::RandomSequence - Creates a random nuc or prot sequence with given nuc frequencies

=head1 VERSION

This document describes BioX::SeqUtils::RandomSequence version 0.9.3

=head1 SYNOPSIS

The package includes scripts for random nucleotide, dinucleotide, protein, and protein set. The length and frequency parameters should always be integers.

To create a nucleotide:

    ./random-nucleotide.pp                               # Defaults: length 60, all frequencies .25
    ./random-nucleotide.pp -l2200 -a23 -c27 -g27 -t23    # Enrich GC content with length 2200

To create a dinucleotide:

    ./random-dinucleotide.pp                             # Defaults: length 2, all frequencies .25
    ./random-dinucleotide.pp -a225 -c275 -g275 -t225     # Enrich GC content ~ more 

To create a protein:

    ./random-protein.pp                                  # Defaults: length 60, all frequencies .25
    ./random-protein.pp -l2200 -a23 -c27 -g27 -t23       # Enrich underlying GC content, aa length 2200

To create a protein set (with common DNA shifted by one base):

    ./random-protein-set.pp                              # Defaults: length 60, all frequencies .25
    ./random-protein-set.pp -l2200 -a23 -c27 -g27 -t23   # Enrich underlying GC content 

Additionally, a "master script" uses a tYpe parameter for any:

    ./random-sequence.pp -yn -l100                       # Type n nucleotide
    ./random-sequence.pp -yd                             # Type d dinucleotide
    ./random-sequence.pp -yp -l100                       # Type p protein
    ./random-sequence.pp -ys -l100                       # Type s protein set

This module uses Bio::Tools::CodonTable for translations, and the parameter s can be used to change from the default (1) Standard:

    ./random-protein.pp -l2200 -s2                       # Non-standard codon table

In script, each sequence type can be accessed using the "y" (tYpe) parameter with rand_seq(). The default is "nucleotide". The type may be set in new() or any of the rand_X() methods. All four frequencies are set to "1" by default ( so that the probablity of each A, C, G, T is 0.25 ).

    use BioX::SeqUtils::RandomSequence;

    my $randomizer = BioX::SeqUtils::RandomSequence->new({ l => $length, 
                                                           s => 1,
                                                           y => "nucleotide",
                                                           a => $a_frequency,
                                                           c => $c_frequency,
                                                           g => $g_frequency,
                                                           t => $t_frequency });
    print $randomizer->rand_seq(), "\n";

You can use the same randomizer object to create all types of sequences, by passing the changing parameters with each call.

    my $nuc_short     = $randomizer->rand_seq({ y => 'n', l => 21 });
    my $nuc_long      = $randomizer->rand_seq({ l => 2200 });          # Still nucleotide
    my $nuc_richer    = $randomizer->rand_seq({ a => 225, 
                                                c => 275, 
						g => 275, 
						t => 225 });           # Still length 2200
    my $protein_now   = $randomizer->rand_seq({ y => 'p' });           # Still richer GC
    my $dinuc_for_fun = $randomizer->rand_seq({ y => 'd',
                                                a => 1 });             # Missing bases resets all freq to 1
    my $protein_new   = $randomizer->rand_seq({ y => 'p',
                                                s => 3 });             # Use codon table 'Yeast Mitochondrial'

Type "protein" creates a protein of the given length l by creating a random nucleotide sequence with the given nucleotide frequencies of length l * 3, which is translated into a protein. The default length is 60. 
    
    my $randomizer = BioX::SeqUtils::RandomSequence->new();
    print $randomizer->rand_seq({ y = "protein" }), "\n";
    
Type "set" creates a test protein set each with the given length l by creating a random nucleotide sequence with the given nucleotide frequencies of length l * 3 + 1, removing the first base for sequence 1 and removing the last base for sequence 2, then translating them into proteins. 

    print join( " ", @{ $randomizer->rand_seq({ y = "set" }) }, "\n";
  
The indvidual methods may be preferred:

    my $nucleotide    = $randomizer->rand_nuc();
    my $dinucleotide  = $randomizer->rand_nuc({ l => 2 });
    my $protein       = $randomizer->rand_pro();

The rand_pro_set() method uses wantarray(), and will either return a list or list reference (scalar) depending on the context:

    my ($pro1, $pro2) = $randomizer->rand_pro_set();
    my $protein_set   = $randomizer->rand_pro_set();

=head1 DESCRIPTION

Create random nucleotide and protein sequences.

=head1 CONFIGURATION AND ENVIRONMENT

None.

=head1 DEPENDENCIES

Class::Std;
Class::Std::Utils;
Bio::Tools::CodonTable;

=head1 INCOMPATIBILITIES

None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to
C<bug-biox-sequtils-randomsequence@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.

=head1 AUTHOR

Roger A Hall  C<< <rogerhall@cpan.org> >>

=head1 LICENSE AND COPYRIGHT

Copyleft (c) 2009, Roger A Hall C<< <rogerhall@cpan.org> >>. All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

