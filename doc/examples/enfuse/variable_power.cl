// Copyright (C) 2015-2017 Christoph L. Spiel
//
// This file is part of Enblend.
//
// Enblend is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Enblend is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Enblend; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// Note: The functionality encoded in this OpenCL-file does *not*
// require "variable_power.cc", a C++-file which defines the same
// weight function using the shared-object extension mechanism of
// Enfuse.


// add missing constant if necessary
#ifndef M_LN2
#define M_LN2 0.69314718f
#endif


// default
#ifndef EXPONENT
#define EXPONENT 2.0f
#endif


float
normalized_luminance(float y)
{
    const float fwhm = 2.0f / exp(M_LN2 / EXPONENT);

    return enfuse_normalized_luminance(y) * fwhm / ENFUSE_FWHM_GAUSSIAN;
}


float
weight(float y)
{
    return max(1.0f - pow(fabs(normalized_luminance(y)), EXPONENT), 0.0f);
}
