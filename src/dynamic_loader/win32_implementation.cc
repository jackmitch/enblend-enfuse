/*
 * Copyright (C) 2015 Thomas Modes
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


#include <string>

#include <boost/lexical_cast.hpp>

#include "global.h"  // for enblend::trim

#include "win32_implementation.h"


WinDynamicLoaderImplementation::WinDynamicLoaderImplementation(const std::string& a_library_name) :
    super(a_library_name), handle_(nullptr)
{}


void
WinDynamicLoaderImplementation::open()
{
    if (handle_)
    {
        throw super::error("already open");
    }
    handle_ = LoadLibrary(library_name().c_str());
    if (!handle_)
    {
        throw super::error(GetLastErrorString());
    }
}


void
WinDynamicLoaderImplementation::close()
{
    if (!handle_)
    {
        throw super::error("not open");
    }
    if (!FreeLibrary(handle_))
    {
        throw super::error(GetLastErrorString());
    }
}


void*
WinDynamicLoaderImplementation::resolve(const std::string& symbol_name) const
{
    if (!handle_)
    {
        throw super::error("not open");

    }
    void* symbol = (void*)GetProcAddress(handle_, symbol_name.c_str());
    if (symbol == nullptr)
    {
        throw super::error(GetLastErrorString());
    }

    return symbol;
}


const std::string
WinDynamicLoaderImplementation::GetLastErrorString()
{
    LPTSTR lpMsgBuf;
    DWORD lastError = GetLastError();
    std::string errorMsg;

    if (FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                      nullptr, lastError, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPTSTR)&lpMsgBuf, 0, nullptr) > 0)
    {
        errorMsg = lpMsgBuf;

        // remove leading and trailing white spaces
        enblend::trim(errorMsg);
        LocalFree(lpMsgBuf);
    }
    else
    {
        errorMsg = "Unknown error";
    }

    // add numeric error code
    errorMsg.append(" (Code: ");
    errorMsg.append(boost::lexical_cast<std::string>(lastError));
    errorMsg.append(")");

    return errorMsg;
}
