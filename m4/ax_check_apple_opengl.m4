# -*- mode: autoconf -*-
#
# AX_CHECK_APPLE_OPENGL
#
# Check for Apple's OpenGL framework.
#
# version: 1.0
# author: Dr. Christoph L. Spiel <cspiel@users.sourceforge.net>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# As a special exception, the you may copy, distribute and modify the
# configure scripts that are the output of Autoconf when processing
# the Macro.  You need not follow the terms of the GNU General Public
# License when using or distributing such scripts.
#
AC_DEFUN([AX_CHECK_APPLE_OPENGL],
[

AC_DEFINE([HAVE_APPLE_OPENGL_FRAMEWORK], 1,
          [Use the Apple OpenGL framework.])

AC_PATH_X
ACX_PTHREAD

# GL
dnl GL_CFLAGS=
GL_LIBS="${GL_LIBS} -framework OpenGL -framework AGL"

# GLU
AC_CACHE_CHECK([for OpenGL Utility library],
               [ax_cv_check_glu_libglu],
               [ax_cv_check_glu_libglu="no"
                ax_save_CPPFLAGS="${CPPFLAGS}"
                CPPFLAGS="${GL_CFLAGS} ${CPPFLAGS}"
                ax_save_LIBS="${LIBS}"
                LIBS=""
                ax_check_libs="-lglu32 -lGLU"
                for ax_lib in ${ax_check_libs}; do
                    if test X$ax_compiler_ms = Xyes; then
                        ax_try_lib=`echo $ax_lib | sed -e 's/^-l//' -e 's/$/.lib/'`
                    else
                        ax_try_lib="${ax_lib}"
                    fi
                    LIBS="${ax_try_lib} ${GL_LIBS} ${ax_save_LIBS}"
                    dnl libGLU typically links with libstdc++ on POSIX
                    dnl platforms.  However, setting the language to
                    dnl C++ means that test program source is named
                    dnl "conftest.cc"; and Microsoft cl does not know
                    dnl what to do with such a file.
                    AC_LANG_PUSH([C++])
                    if test X$ax_compiler_ms = Xyes; then
                        AC_LANG_PUSH([C])
                    fi
                    AC_LINK_IFELSE(
                    [AC_LANG_PROGRAM(
[[
#if HAVE_WINDOWS_H && defined(_WIN32)
#include<windows.h>
#endif
#include <GL/glu.h>
]],
[[
gluBeginCurve(0)
]])],
                    [ax_cv_check_glu_libglu="${ax_try_lib}"; break])
                    if test X$ax_compiler_ms = Xyes; then
                        AC_LANG_POP([C])
                    fi
                    AC_LANG_POP([C++])
                done
                LIBS=${ax_save_LIBS}
                CPPFLAGS=${ax_save_CPPFLAGS}])

if test "X${ax_cv_check_glu_libglu}" = "Xno"; then
    no_glu="yes"
    GLU_CFLAGS=""
    GLU_LIBS=""
else
    GLU_LIBS="${ax_cv_check_glu_libglu} ${GL_LIBS}"
fi

# GLUT
GLUT_CFLAGS="${GLU_CFLAGS}"
GLUT_LIBS="-framework GLUT -lobjc ${GL_LIBS}"

# finish
AC_SUBST([GL_CFLAGS])
AC_SUBST([GL_LIBS])
AC_SUBST([GLU_CFLAGS])
AC_SUBST([GLU_LIBS])
AC_SUBST([GLUT_CFLAGS])
AC_SUBST([GLUT_LIBS])

])dnl
