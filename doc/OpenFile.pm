# name:         OpenFile.pm
# synopsis:     open file or device in the OO way
# author:       Dr. Christoph L. Spiel
# perl version: 5.20.2


# This file is part of Enblend.
# Licence details can be found in the file COPYING.


package OpenFile;

use strict;
use warnings;

use Carp;
use English;
use IO::File;
use IO::Handle;

use Quote ();


sub open_file {
    my ($file, $mode) = @_;

    my $handle;

    $mode = 'r' unless $mode;

    if ($file eq '-') {
        $handle = new IO::Handle;
        if ($mode eq 'r') {
            $handle->fdopen(fileno(STDIN), 'r') or
              croak("cannot open standard input: $OS_ERROR");
        } elsif ($mode eq 'w') {
            $handle->fdopen(fileno(STDOUT), 'w') or
              croak("cannot open standard output: $OS_ERROR");
        } else {
            croak("internal error, unknown file mode ", Quote::gnu_style($mode));
        }
    } else {
        $handle = new IO::File($file, $mode) or
          croak("cannot open ", Quote::gnu_style($file), ": $OS_ERROR");
    }

    return $handle;
}


1;
