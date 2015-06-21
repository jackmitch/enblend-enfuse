# name:         TexAux.pm
# synopsis:     Helper functions for TeX/LaTeX processing
# author:       Dr. Christoph L. Spiel
# perl version: 5.20.2


# This file is part of Enblend.
# Licence details can be found in the file COPYING.


package TexAux;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(escape lift);


sub escape {
    my $string = shift;

    $string =~ s#([\$&%_^{}])#\\$1#g;

    return $string;
}


sub lift {
    my $string = shift;

    if ($string =~ m#^https?:#) {
        $string = "\\url{$string}";
    }

    return $string;
}


1;
