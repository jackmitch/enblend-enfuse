/*
 * Copyright (C) 2009-2015 Dr. Christoph L. Spiel
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
#include <sstream>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "error_message.h"

namespace enblend
{
namespace detail
{
    std::string
    noErrorMessageAvailable(int anErrorNumber)
    {
        std::ostringstream oss;
        oss << "No detailed error message available; decimal error #" << anErrorNumber;
        return oss.str();
    }
}


std::string
errorMessage(int anErrorNumber)
{
#if !defined(HAVE_STRERROR) && !defined(HAVE_STRERROR_R)
    return detail::noErrorMessageAvailable(anErrorNumber);
#endif

#if !defined(HAVE_STRERROR_R)
    std::string message(strerror(anErrorNumber));
#else
    const size_t buffer_size = 65536;
    std::string message(buffer_size, '\0');
#if defined(STRERROR_R_CHAR_P)
    char* messageptr = strerror_r(anErrorNumber, &message[0], buffer_size);
    if (messageptr != &message[0])
    {
        message = std::string(messageptr);
    }
    else
    {
        message.resize(message.find('\0'));
    }
#else
    int return_code = strerror_r(anErrorNumber, &message[0], buffer_size);
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

    if (message.length() == 0)
    {
        return detail::noErrorMessageAvailable(anErrorNumber);
    }
    else
    {
        return message;
    }
}

}

