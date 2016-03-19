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
our @EXPORT = qw(escape lift default_name);


sub escape {
    my $string = shift;

    $string =~ s#([\$&%_^{}])#\\$1#g;

    return $string;
}


sub lift {
    my $string = shift;

    if ($string =~ m#^https?:#) {
        $string = "\\url{$string}";
    } else {
        $string =~ s#:#:\\feasiblebreak #g;
        $string =~ s#(?<!^)/(?!$)#\\slash #g;
    }

    return $string;
}


my %defaults =
  (EMPTY_MACRO_VALUE => '\\definedbypreprocessor',
   FIND_MACRO_NAME => '\\hashfind',
   INSERT_MACRO_NAME => '\\hashinsert',
   KEY_PREFIX => 'val:',
   UNDEF_MACRO_VALUE => '\\undefinedbypreprocessor');


sub default_name {
    my $id = shift;

    return $defaults{$id};
}


1;
