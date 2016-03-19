/*
 * Copyright (C) 2016 T. Modes
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

// handling of error message for delayed loading of opencl.dll

#ifndef DELAYHELPER_H_INCLUDED
#define DELAYHELPER_H_INCLUDED

#pragma comment(lib, "delayimp")

#include <iostream>
#define _STLP_VERBOSE_AUTO_LINK
#define _USE_MATH_DEFINES
#define NOMINMAX
#define VC_EXTRALEAN
#include <windows.h>
#undef DIFFERENCE
#include <delayimp.h>
static char* openclLib = "opencl.dll";

FARPROC WINAPI delayHookFailureFunc(unsigned dliNotify, PDelayLoadInfo pdli)
{
    switch (dliNotify)
    {
        case dliFailLoadLib:
            // LoadLibrary failed.
            if (pdli)
            {
                std::cerr << command << ": error: Could not load \"" << pdli->szDll << "\".\n";
                if (_strnicmp(pdli->szDll, openclLib, strlen(openclLib)) == 0)
                {
                    std::cerr << command << ": info: Remove option \"--gpu\" and try again.\n";
                }
            }
            else
            {
                std::cerr << command << ": error: Could not load an unknown dll.\n";
            }
            exit(1);
            break;

        case dliFailGetProc:
            // GetProcAddress failed.
            if (pdli)
            {
                if (pdli->dlp.fImportByName)
                {
                    std::cerr << command << ": error: Function \"" << pdli->dlp.szProcName << "\" was not found in \"" << pdli->szDll << "\".\n";
                }
                else
                {
                    std::cerr << command << ": error: Function ordinal " << pdli->dlp.dwOrdinal << " was not found in \"" << pdli->szDll << "\".\n";
                }
            }
            else
            {
                std::cerr << command << ": error: Unknown error in GetProcAddress for unknown dll.\n";
            }
            exit(1);
            break;
    }
    return (FARPROC)NULL;
}

PfnDliHook __pfnDliFailureHook2;

#endif // DELAYHELPER_H_INCLUDED