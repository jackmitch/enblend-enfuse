#include <cstring>
#include <sstream>

#include <boost/scoped_ptr.hpp>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "error_message.h"

namespace enblend
{

std::string
errorMessage(int anErrorNumber)
{
#if HAVE_STRERROR_R || HAVE_STRERROR
#if HAVE_STRERROR_R
    const size_t size = 4096;
    boost::scoped_ptr<char> message_buffer(new char[size]);
    const char* message = strerror_r(anErrorNumber, message_buffer.get(), size);
#elif HAVE_STRERROR
    const char* message = strerror(anErrorNumber);
#endif

    if (strlen(message) == 0)
    {
        std::ostringstream oss;
        oss << "No detailed error message available, error #" << anErrorNumber;
        return oss.str();
    }
    else
    {
        return std::string(message);
    }
#else
    std::ostringstream oss;
    oss << "No detailed error message available, error #" << anErrorNumber;
    return oss.str();
#endif
}

}
