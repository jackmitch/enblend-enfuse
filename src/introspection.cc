/*
 * Copyright (C) 2015, 2016 Christoph L. Spiel
 *
 * This file is part of Enblend.
 *
 * Enblend is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Enblend is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Enblend; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include <iostream>
#include <string>
#include <regex>

#include <gsl/gsl_version.h>             // GSL_VERSION
#include <lcms2.h>                       // LCMS_VERSION
#include <vigra/imageinfo.hxx>           // VIGRA_VERSION, impexListExtensions(), impexListFormats()

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "filespec.h"
#include "global.h"
#include "signature.h"
#include "dynamic_loader.h"    // HAVE_DYNAMICLOADER_IMPL
#include "openmp_def.h"        // OPENMP
#include "opencl.h"            // OPENCL

#include "introspection.h"

#ifdef _MSC_VER
#include <delayimp.h>
#endif

extern const std::string command;
extern Signature sig;
extern int Verbose;


namespace introspection
{
    void
    printVersion(int argc, char** argv)
    {
        std::cout << command << " " << VERSION << "\n\n";

        if (Verbose >= VERBOSE_VERSION_REPORTING)
        {
#ifdef ENFUSE_SOURCE
            {
                std::cout << "Extra feature: dynamic linking support: " <<
#ifdef HAVE_DYNAMICLOADER_IMPL
                    "yes" <<
#else
                    "no" <<
#endif
                    "\n";
            }
#endif // ENFUSE_SOURCE

// #ifdef CACHE_IMAGES
//             std::cout << "Extra feature: image cache: yes\n";
//             {
// #ifdef WIN32
//                 char lpPathBuffer[MAX_PATH];
//                 const DWORD dwRetVal = GetTempPath(MAX_PATH, lpPathBuffer);
//                 if (dwRetVal <= MAX_PATH && dwRetVal != 0)
//                 {
//                     std::cout << "  - cache file located in \"" << lpPathBuffer << "\"\n";
//                 }
// #else
//                 const char* tmpdir = getenv("TMPDIR");
//                 std::cout << "  - environment variable TMPDIR ";
//                 if (tmpdir == nullptr)
//                 {
//                     std::cout << "not set, cache file in default directory \"/tmp\"\n";
//                 } else
//                 {
//                     std::cout << "set, cache file located in \"" << tmpdir << "\"\n";
//                 }
// #endif
//             }
// #else
            std::cout << "Extra feature: image cache: no\n";
// #endif

#ifdef OPENMP
            const bool have_dynamic = have_openmp_dynamic();
            std::cout <<
                "Extra feature: OpenMP: yes\n" <<
                "  - version " << OPENMP_YEAR << '-' << OPENMP_MONTH << "\n" <<
                "  - " << (have_dynamic ? "" : "no ") <<
                "support for dynamic adjustment of the number of threads;\n" <<
                "    dynamic adjustment " <<
                (have_dynamic && omp_get_dynamic() ? "enabled" : "disabled") << " by default\n" <<
                "  - using " <<
                omp_get_num_procs() << " processor" << (omp_get_num_procs() >= 2 ? "s" : "") << " and up to " <<
                omp_get_max_threads() << " thread" << (omp_get_max_threads() >= 2 ? "s" : "") << "\n" <<
                "  - allocating thread-local dynamic memory with " << OMP_MALLOC_FUNCTIONS << "\n";
#else
            std::cout << "Extra feature: OpenMP: no\n";
#endif

#ifdef OPENCL
            std::cout << "Extra feature: OpenCL: yes\n";
#ifdef _MSC_VER
            // catching errors with delay loading of opencl.dll
            // remember to reset __pfnDliFailureHook2 before
            // otherwise the failure hooks exits the program before __except is reached
            __try
            {
                ocl::print_opencl_information();
            }
            __except(GetExceptionCode() == VcppException(ERROR_SEVERITY_ERROR, ERROR_MOD_NOT_FOUND) ||
                GetExceptionCode() == VcppException(ERROR_SEVERITY_ERROR, ERROR_PROC_NOT_FOUND))
            {
                std::cout << "                       but not available on this system\n";
            }
#else
            ocl::print_opencl_information();
#endif
#else
            std::cout << "Extra feature: OpenCL: no\n";
#endif
            std::cout << "\n";
        }

        std::cout <<
            "Copyright (C) 2004-2009 Andrew Mihal.\n" <<
            "Copyright (C) 2009-2016 Christoph Spiel.\n" <<
            "\n" <<
            "License GPLv2+: GNU GPL version 2 or later <http://www.gnu.org/licenses/gpl.html>\n" <<
            "This is free software: you are free to change and redistribute it.\n" <<
            "There is NO WARRANTY, to the extent permitted by law.\n" <<
            "\n" <<
            "Written by Andrew Mihal, Christoph Spiel and others." <<
            std::endl;

        exit(0);
    }


    // With the exception of TIFF, VIFF, PNG, and PNM all file types
    // store only 1 byte (gray scale and mapped RGB) or 3 bytes (RGB)
    // per pixel.
    //
    // PNG can store UInt8 and UInt16 values, and supports 1 and 3
    // channel images. One additional alpha channel is also supported.
    //
    // PNM can store 1 and 3 channel images with UInt8, UInt16 and
    // UInt32 values in each channel.
    //
    // TIFF and VIFF are additionally able to store short and long
    // integers (2 or 4 bytes) and real values (32 bit float and 64 bit
    // double) without conversion.

    void
    printImageFormats()
    {
        std::string formats(vigra::impexListFormats());
        std::string extensions(vigra::impexListExtensions());

        const std::regex space(" ");
        const std::string replacement("\n  ");

        formats = std::regex_replace(formats, space, replacement);
        extensions = std::regex_replace(extensions, space, replacement);

        std::cout <<
            "Following image formats are supported by " << command << ":\n" <<
            "  " << formats << "\n" <<
            "It automatically recognizes the image file extensions:\n" <<
            "  " << extensions << "\n" <<
            "and -- where supported by the image format -- accepts per-channel depths:\n" <<
#ifndef DEBUG_8BIT_ONLY
            "  " << "8 bits unsigned integral\n" <<
#endif
            "  " << "16 bits unsigned or signed integral\n" <<
            "  " << "32 bits unsigned or signed integral\n" <<
            "  " << "32 bits floating-point\n" <<
            "  " << "64 bits floating-point\n" <<
            std::endl;

        exit(0);
    }


    void
    printSignature()
    {
        std::cout.flush();
        std::wcout << sig.message() << L"\n\n";
        std::wcout.flush();

        exit(0);
    }


    void
    printGlobbingAlgos()
    {
        const enblend::algorithm_list algos = enblend::known_globbing_algorithms();

        std::cout << "Following globbing algorithms are supported:\n";
        for (auto i = algos.begin(); i != algos.end(); ++i)
        {
            std::cout <<
                "  " << i->first << "\n" <<
                "    " << i->second << "\n";
        }
        std::cout << std::endl;

        exit(0);
    }


    void
    printSoftwareComponents()
    {
        std::cout << "Compiler\n  " <<
            // IMPLEMENTATION NOTE: Order matters as both CLang and ICC
            // identify themselves as GCC, too.
#if defined(__clang__)
            "clang++ " << __clang_major__ << '.' << __clang_minor__ << '.' << __clang_patchlevel__ <<
#elif defined(__ICC) || defined(__INTEL_COMPILER)
            "icpc " << __GNUC__ << '.' << __GNUC_MINOR__ << '.' << __GNUC_PATCHLEVEL__ <<
#elif defined(__GNUG__)
            "g++ " << __GNUC__ << '.' << __GNUC_MINOR__ << '.' << __GNUC_PATCHLEVEL__ <<
#elif defined(__IBMC__) || defined(__IBMCPP__)
            "xlc++ " << __IBMCPP__ <<
#elif defined(_MSC_VER)
            "MS Visual Studio " << _MSC_FULL_VER <<
#elif defined(__PGI)
            "pgcpp " << __PGIC__ << '.' << __PGIC_MINOR__ << '.' << __PGIC_PATCHLEVEL__ <<
#else
            "unknown" <<
#endif
            "\n" <<
#ifdef OPENMP
            "  implementing OpenMP standard of " << _OPENMP / 100 << '-' << _OPENMP % 100 << "\n" <<
#endif
            "\n";

#ifdef OPENCL
        std::cout << "OpenCL APIs\n" <<
#ifdef CL_VERSION_1_0
            "  1.0\n" <<
#endif
#ifdef CL_VERSION_1_1
            "  1.1\n" <<
#endif
#ifdef CL_VERSION_1_2
            "  1.2\n" <<
#endif
#ifdef CL_VERSION_1_3
            "  1.3\n" <<
#endif
#ifdef CL_VERSION_2_0
            "  2.0\n" <<
#endif
#ifdef CL_VERSION_2_1
            "  2.1\n" <<
#endif
            "\n";
#endif // OPENCL

        std::cout << "Libraries\n" <<
            "  GSL:        " << GSL_VERSION << "\n" <<
            //"  JPEG:       " << "\n" <<
            "  Little CMS: " << LCMS_VERSION / 1000 << '.' << (LCMS_VERSION / 10) % 100 << '.' << LCMS_VERSION % 10 << "\n" <<
            //"  PNG:        " << "\n" <<
            //"  OpenEXR:    " << "\n" <<
            //"  TIFF:       " << "\n" <<
            "  Vigra:      " << VIGRA_VERSION << "\n" <<
            std::endl;

        exit(0);
    }
} // end namespace introspection
