// Copyright (C) 2014 Christoph L. Spiel
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
#include <cerrno>                 // errno
#include <cmath>                  // macro M_LN2, std::exp(), std::abs(), std::pow()

#include "exposure_weight_base.h" // macro FWHM_GAUSSIAN, class ExposureWeight


class VariablePower : public ExposureWeight
{
    typedef ExposureWeight super;

public:
    void initialize(double y_optimum, double width, const argument_list_t& argument_list) override
    {
        if (argument_list.empty())
        {
            exponent = 2.0;
        }
        else
        {
            char* tail;

            errno = 0;
            exponent = strtod(argument_list[0].c_str(), &tail);
            if (*tail != 0 || errno != 0)
            {
                throw super::error("non-numeric exponent");
            }
            if (exponent <= 0.0 || exponent >= 100.0)
            {
                throw super::error("exponent out of range ]0, 100]");
            }
        }

        const double fwhm = 2.0 / std::exp(M_LN2 / exponent);

        super::initialize(y_optimum, width * FWHM_GAUSSIAN / fwhm, argument_list);
    }


    double weight(double y) override
    {
        return std::max(1.0 - std::pow(std::abs(normalize(y)), exponent), 0.0);
    }

private:
    double exponent;
};


VariablePower vpower;
