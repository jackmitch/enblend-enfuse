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
#ifndef __ERROR_MESSAGE_H__
#define __ERROR_MESSAGE_H__

#include <string>

namespace enblend
{
    /** Answer the error message associated with anErrorNumber. */
    std::string errorMessage(int anErrorNumber);
}

#endif /* __ERROR_MESSAGE_H__ */

// Local Variables:
// mode: c++
// End:
