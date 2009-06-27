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
#ifndef __GLOBAL_H__
#define __GLOBAL_H__

// Here we define macros and types that we already need in the
// definitions of global variables.


#include <sstream>


#define DEFAULT_OUTPUT_FILENAME "a.tif"


struct AlternativePercentage {
    double value;
    bool isPercentage;

    std::string str() const {
        std::ostringstream oss;
        oss << value;
        if (isPercentage) {oss << "%";}
        return oss.str();
    }
};


#define DEFAULT_TIFF_RESOLUTION 300.0f


struct TiffResolution {
    TiffResolution() : x(0.0f), y(0.0f) {}

    TiffResolution(float anXresolution, float aYresolution) :
        x(anXresolution), y(aYresolution) {}

    bool operator==(const TiffResolution& anOther) const {
        return this->x == anOther.x && this->y == anOther.y;
    }

    bool operator!=(const TiffResolution& anOther) const {
        return !operator==(anOther);
    }

    float x;
    float y;
};

#endif /* __GLOBAL_H__ */

// Local Variables:
// mode: c++
// End:
