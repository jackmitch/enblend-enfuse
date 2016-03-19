# name:         Quote.pm
# synopsis:     Put strings in quotes
# author:       Dr. Christoph L. Spiel
# perl version: 5.20.2


# This file is part of Enblend.
# Licence details can be found in the file COPYING.


package Quote;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(quote quote2 gnu_style);


sub quote {
    my ($quote, @strings) = @_;

    map {qq($quote$_$quote)} @strings;
}


sub quote2 {
    my ($opening_quote, $closing_quote, @strings) = @_;

    map {qq($opening_quote$_$closing_quote)} @strings;
}


sub gnu_style {
    quote2(q(`), q('), @_);
}


1;
