/*
** Copyright (C) 2001 Erik de Castro Lopo <erikd AT mega-nerd DOT com>
**
** Permission to use, copy, modify, distribute, and sell this file for any
** purpose is hereby granted without fee, provided that the above copyright
** and this permission notice appear in all copies.  No representations are
** made about the suitability of this software for any purpose.  It is
** provided "as is" without express or implied warranty.
*/

/* Version 1.1 */


/*============================================================================
** On Intel Pentium processors (especially PIII and probably P4), converting
** from float to int is very slow. To meet the C specs, the code produced by
** most C compilers targeting Pentium needs to change the FPU rounding mode
** before the float to int conversion is performed.
**
** Changing the FPU rounding mode causes the FPU pipeline to be flushed. It
** is this flushing of the pipeline which is so slow.
**
** Fortunately the ISO C99 specifications define the functions lrint, lrintf,
** llrint and llrintf which fix this problem as a side effect.
**
** On Unix-like systems, the configure process should have detected the
** presence of these functions. If they weren't found we have to replace them
** here with a standard C cast.
*/

/*
** The C99 prototypes for lrint and lrintf are as follows:
**
** long int lrintf (float x) ;
** long int lrint  (double x) ;
*/

//#include "config.h"

/* The presence of the required functions are detected during the configure
** process and the values HAVE_LRINT and HAVE_LRINTF are set accordingly in
** the config.h file.
*/

#ifndef __FLOAT_CAST_H__
#define __FLOAT_CAST_H__

#if (HAVE_LRINT && HAVE_LRINTF)

/* These defines enable functionality introduced with the 1999 ISO C
** standard. They must be defined before the inclusion of math.h to
** engage them. If optimisation is enabled, these functions will be
** inlined. With optimisation switched off, you have to link in the
** maths library using -lm.
*/

#define _ISOC9X_SOURCE  1
#define _ISOC99_SOURCE  1

#define __USE_ISOC9X    1
#define __USE_ISOC99    1

#include        <math.h>

#elif (defined (_MSC_VER))

/* Win32 doesn't seem to have these functions.
** Therefore implement inline versions of these functions here.
*/

#if defined(_M_X64)

#include <intrin.h>

__forceinline long long int llrint (double flt) {
    return (long long int) _mm_cvtsd_si64x(_mm_loadu_pd(&flt));
}

__forceinline long int lrint (double flt) {
    return (long int) _mm_cvtsd_si64x(_mm_loadu_pd(&flt));
}

#else

__forceinline long int lrint (double flt) {
    long int intgr;
    __asm {
        fld flt
        fistp intgr
    }
    return intgr;
}

__forceinline long int lrintf(float flt) {
    long int intgr;
    __asm {
        fld flt
        fistp intgr
    }
    return intgr;
}

__forceinline long int lrintl(long double flt) {
    long int intgr;
    __asm {
        fld flt
        fistp intgr
    }
    return intgr;
}

__forceinline long long int llrint(double flt) {
    long long int intgr;
    __asm {
        fld flt
        fistp intgr
    }
    return intgr;
}

__forceinline long long int llrintf(float flt) {
    long long int intgr;
    __asm {
        fld flt
        fistp intgr
    }
    return intgr;
}

__forceinline long long int llrintl(long double flt) {
    long long int intgr;
    __asm {
        fld flt
        fistp intgr
    }
    return intgr;
}

#endif /* _M_X64 */

#endif

#endif /* __FLOAT_CAST_H__ */
