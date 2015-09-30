#include "win32_implementation.h"


#ifdef WIN32


#include <string>

#include <boost/lexical_cast.hpp>

#include "global.h"  // for enblend::trim


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


#endif // WIN32
