# SYNOPSIS
#
#   AX_CXX_AS_CONST
#
# DESCRIPTION
#
#   If the C++ compiler supports std::as_const(), define HAVE_AS_CONST.
#
# LICENSE
#
#   Copyright (c) 2015 Christoph L. Spiel <cspiel@users.sourceforge.net>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 1

AU_ALIAS([AC_CXX_AS_CONST], [AX_CXX_AS_CONST])
AC_DEFUN([AX_CXX_AS_CONST],
[dnl
  AC_CACHE_CHECK([whether the compiler supports as_const()], ax_cv_cxx_as_const,
  [dnl
    AC_LANG_PUSH([C++])
    AC_COMPILE_IFELSE([dnl
      AC_LANG_PROGRAM(dnl
        [
#include <string>
#include <type_traits>
#include <utility>
        ],
        [std::string s("foo");
const std::string& const_s = std::as_const(s)])],
      ax_cv_cxx_as_const=yes,
      ax_cv_cxx_as_const=no)
    AC_LANG_POP([C++])
  ])
  AS_IF([test "X$ax_cv_cxx_as_const" = Xyes],
    [AC_DEFINE(HAVE_AS_CONST,,[define if the compiler supports as_const()])])
])
