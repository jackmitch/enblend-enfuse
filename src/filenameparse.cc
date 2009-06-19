/*
 * Copyright (C) 2009 Dr. Christoph L. Spiel
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


// Life is tough and then you die.  -- Jack Dempsey


#include <string>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "filenameparse.h"


#ifdef HAVE_WINDOWS_H
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif


#ifdef HAVE_BOOST_FILESYSTEM
#include "boost/filesystem.hpp"

typedef boost::filesystem::basic_path<std::string, boost::filesystem::path_traits> basic_path;
#endif


namespace enblend {

std::string
extractDirname(const std::string& aFilename)
{
#ifdef HAVE_BOOST_FILESYSTEM
    const basic_path path(aFilename);
    const std::string directory(path.branch_path().string());
    return directory.empty() ? "." : directory;
#else
    const std::string::size_type separator = aFilename.rfind(PATH_SEPARATOR);
    return (separator == std::string::npos) ? "." : aFilename.substr(0, separator);
#endif
}


std::string
extractBasename(const std::string& aFilename)
{
#ifdef HAVE_BOOST_FILESYSTEM
    const basic_path path(aFilename);
    return path.leaf();
#else
    const std::string::size_type separator = aFilename.rfind(PATH_SEPARATOR);
    return
        (separator == std::string::npos) ?
        aFilename :
        aFilename.substr(separator + 1, aFilename.length() - separator - 1);
#endif
}


std::string
extractFilename(const std::string& aFilename)
{
#ifdef HAVE_BOOST_FILESYSTEM
    const basic_path path(aFilename);
    return basename(path);
#else
    const std::string::size_type separator = aFilename.rfind(PATH_SEPARATOR);
    const std::string::size_type dot = aFilename.rfind(".");
    if (separator == std::string::npos)
    {
        return (dot == std::string::npos) ? aFilename : aFilename.substr(0, dot);
    }
    else
    {
        return
            (dot == std::string::npos) ?
            aFilename.substr(separator + 1, aFilename.length() - separator - 1) :
            aFilename.substr(separator + 1, dot - separator - 1);
    }
#endif
}


std::string
extractExtension(const std::string& aFilename)
{
#ifdef HAVE_BOOST_FILESYSTEM
    const basic_path path(aFilename);
    return extension(path);
#else
    const std::string::size_type dot = aFilename.rfind(".");
    return
        (dot == std::string::npos) ?
        "" :
        aFilename.substr(dot, aFilename.length() - dot);
#endif
}

} // namespace enblend

// Local Variables:
// mode: c++
// End:
