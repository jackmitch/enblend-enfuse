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
#ifndef __FILENAMEPARSE_H__
#define __FILENAMEPARSE_H__

#include <string>

namespace enblend {
    /** Answer the directory part of aFilename.  The function duplicates
     *  dirname(1). */
    std::string extractDirname(const std::string& aFilename);

    /** Answer the non-directory part of aFilename.  The function
     *  duplicates basename(1). */
    std::string extractBasename(const std::string& aFilename);

    /** Answer the filename part of aFilename.  This is the basename of
     *  aFilename without extension. */
    std::string extractFilename(const std::string& aFilename);

    /** Answer the extension part of aFilename including the leading
     *  dot. */
    std::string extractExtension(const std::string& aFilename);
} // namespace enblend

#endif /* __FILENAMEPARSE_H__ */

// Local Variables:
// mode: c++
// End:
