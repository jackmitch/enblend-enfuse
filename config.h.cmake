/*
 * \file config.h
 * This file is part of enblend.
 * Licence details can be found in the file COPYING.
 *
 * This is the compilation configuration file for enblend.
 * You might want to change some of the defaults if something goes wrong
 * during the compilation.
 */

#ifndef _CONFIG_H
#define _CONFIG_H

/* Default verbosity level of Enblend and Enfuse. A value of zero reduces the
   message output to only warnings and errors. Values of one and more make
  Enblend and Enfuse report progress in detail. */
#define DEFAULT_VERBOSITY 1

${CMAKE_HEADER_EXISTS}
${CMAKE_FUNCTION_EXISTS}

/* Define to 1 if the `closedir' function returns void instead of `int'. */
#cmakedefine CLOSEDIR_VOID 1

/* Define if you have POSIX threads libraries and header files. */
#cmakedefine HAVE_PTHREAD 1

/* Define to 1 if the system has the type `ptrdiff_t'. */
#cmakedefine HAVE_PTRDIFF_T 1

/* define to 1 if you have dlfcn.h header file and dl lib. */
#cmakedefine HAVE_DL 1

/* define to 1 if you have opencl.h header file in "CL" dir. */
#cmakedefine HAVE_CL_CL_HPP 1

/* define to 1 if you have opencl.h header file in "OpenCL" dir. */
#cmakedefine HAVE_OPENCL_CL_HPP 1

/* Define to 1 if the system has the type `_Bool'. */
/* #undef HAVE__BOOL */

/* define the correct restrict keyword, Empty or one of
  __restrict __restrict__ _Restrict restrict */
#define RESTRICT ${RESTRICT}

/* Name of package */
#define PACKAGE "enblend-enfuse"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://bugs.launchpad.net/enblend"

/* Define to the full name of this package. */
#define PACKAGE_NAME "enblend-enfuse"

/* Define to the full name and version of this package. */
#cmakedefine PACKAGE_STRING "${PACKAGE_STRING}"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "enblend-enfuse"

/* Define to the version of this package. */
#define PACKAGE_VERSION "${PACKAGE_VERSION_STRING}"

#define PACKAGE_URL "${PACKAGE_URL}"

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

/* Define as the return type of signal handlers (`int' or `void'). */
#cmakedefine RETSIGTYPE ${RETSIGTYPE}

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS 1

/* Define to 1 if strerror_r returns char *. */
#cmakedefine STRERROR_R_CHAR_P 1

/* workaround for older boost versions <1.55 */
#cmakedefine HAVE_BOOST_FALLTHROUGH 1
#ifndef HAVE_BOOST_FALLTHROUGH
#define BOOST_FALLTHROUGH ((void) 0)
#endif

/* Version number of package */
#define VERSION "${ENBLEND_VERSION_ONLY}"

/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
#cmakedefine WORDS_BIGENDIAN 1

/* Define to `long int' if <sys/types.h> does not define. */
#cmakedefine HAVE_OFF_T 1
#if ! defined(HAVE_OFF_T)
  #define off_t long int
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
#cmakedefine HAVE_SIZE_T 1
#if ! defined(HAVE_SIZE_T)
  #define size_t unsigned int
#endif
#endif

/* Define to the implicit search path for OpenCL kernels. */
#define DEFAULT_OPENCL_PATH "${DEFAULT_OPENCL_PATH}"

/* Prefer separate OpenCL kernels or use build-in strings. */
#cmakedefine PREFER_SEPARATE_OPENCL_SOURCE 1
