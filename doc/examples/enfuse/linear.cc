// Copyright (C) 2014-2016 Christoph L. Spiel
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


#include <cmath>                  // std::abs()
#include <iostream>               // std::cerr

#include "exposure_weight_base.h" // macro FWHM_GAUSSIAN, class ExposureWeight


struct Linear : public ExposureWeight
{
    void initialize(double y_optimum, double width_parameter,
                    ExposureWeight::argument_const_iterator arguments_begin,
                    ExposureWeight::argument_const_iterator arguments_end) override
    {
        if (arguments_begin != arguments_end)
        {
            std::cerr << "warning: weight function \"linear\" does not take any parameters" << std::endl;
        }

        ExposureWeight::initialize(y_optimum, width_parameter * FWHM_GAUSSIAN,
                                   arguments_begin, arguments_end);
    }


    double weight(double y) override
    {
        const double z = std::abs(normalize(y));

        return z <= 1.0 ? 1.0 - z : 0.0;
    }
};


Linear linear;
