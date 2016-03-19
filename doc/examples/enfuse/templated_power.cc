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


#include <algorithm>              // std::max()
#include <cmath>                  // macro M_LN2, std::exp(), std::abs()
#include <iostream>               // std::cerr

#include "exposure_weight_base.h" // macro FWHM_GAUSSIAN, class ExposureWeight


template <int n>
inline static double
ipower(double x)
{
    return x * ipower<n - 1>(x);
}


template <>
double
ipower<0>(double)
{
    return 1.0;
}


template <int n>
struct TemplatedPower : public ExposureWeight
{
    void initialize(double y_optimum, double width,
                    ExposureWeight::argument_const_iterator arguments_begin,
                    ExposureWeight::argument_const_iterator arguments_end) override
    {
        if (arguments_begin != arguments_end)
        {
            std::cerr << "warning: weight function \"power#\" does not take any parameters" << std::endl;
        }

        const double fwhm = 2.0 / std::exp(M_LN2 / static_cast<double>(n));

        ExposureWeight::initialize(y_optimum, width * FWHM_GAUSSIAN / fwhm,
                                   arguments_begin, arguments_end);
    }


    double weight(double y) override
    {
        return std::max(1.0 - ipower<n>(std::abs(normalize(y))), 0.0);
    }
};


TemplatedPower<2> tpower2;
TemplatedPower<3> tpower3;
TemplatedPower<4> tpower4;
