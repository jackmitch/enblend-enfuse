/*
 * Copyright (C) 2009-2014 Dr. Christoph L. Spiel
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


#include <cstring>
#include <string>
#include <sstream>
#ifdef _WIN32
#include <windows.h>
#endif  // _WIN32

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "error_message.h"


namespace enblend
{
    static std::string
    noErrorMessageAvailable(int anErrorNumber)
    {
        std::ostringstream oss;
        oss << "No detailed error message available; decimal error #" << anErrorNumber;
        return oss.str();
    }


    std::string
    errorMessage(int anErrorNumber)
    {
#if !defined(HAVE_STRERROR) && !defined(HAVE_STRERROR_R)
        return noErrorMessageAvailable(anErrorNumber);
#endif

#if !defined(HAVE_STRERROR_R)
        std::string message(strerror(anErrorNumber));
#else
        const size_t buffer_size = 65536U;
        std::string message(buffer_size, '\0');
#if defined(STRERROR_R_CHAR_P)
        char* const message_begin = strerror_r(anErrorNumber, &message[0], buffer_size);
        if (message_begin != &message[0])
        {
            message = std::string(message_begin);
        }
        else
        {
            message.resize(message.find('\0'));
        }
#else
        const int return_code = strerror_r(anErrorNumber, &message[0], buffer_size);
        if (return_code != 0)
        {
            std::ostringstream oss;
            oss <<
                "Conversion of decimal error #" << anErrorNumber <<
                " to an error message failed with decimal return code " << return_code;
            return oss.str();
        }
#endif // STRERROR_R_CHAR_P
#endif // HAVE_STRERROR_R

#if defined(_WIN32)
        // strerror translates only errors in C runtime library
        // now translate errors in Windows API
        if (message.compare("Unknown error") == 0)
        {
            message.clear();
        }
        if (anErrorNumber && message.empty())
        {
            // first check for errors in Windows API
            LPSTR messageBuffer = nullptr;
            DWORD size = FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                NULL, anErrorNumber, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&messageBuffer, 0, NULL);
            if (size)
            {
                // convert codepage to get correct display of all characters, e.g. umlaute
                LPSTR OEMmessageBuffer = (LPSTR)LocalAlloc(LPTR, (size + 1)*sizeof(char));
                if (OEMmessageBuffer)
                {
                    if (CharToOemBuff(messageBuffer, OEMmessageBuffer, size))
                    {
                        message = std::string(OEMmessageBuffer, size);
                    }
                }
                LocalFree(OEMmessageBuffer);
            }
            LocalFree(messageBuffer);
        }
#endif // _WIN32    
        
        return message.empty() ? noErrorMessageAvailable(anErrorNumber) : message;
    }
} // namespace enblend
