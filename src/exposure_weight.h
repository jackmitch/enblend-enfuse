/*
 * Copyright (C) 2013, 2014 Dr. Christoph L. Spiel
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
#ifndef EXPOSURE_WEIGHT_INCLUDED
#define EXPOSURE_WEIGHT_INCLUDED


#include <cmath>

#include "exposure_weight_base.h"


struct Gaussian : public ExposureWeight
{
    Gaussian(double y_optimum, double width_parameter) : ExposureWeight(y_optimum, width_parameter) {}

    double weight(double y) override
    {
        const double z = normalize(y);
        return exp(-0.5 * z * z);
    }
};


struct Lorentzian : public ExposureWeight
{
    // FWHM = 2 * sqrt(2)
#define FWHM_LORENTZIAN 2.8284271247461900976033774484193961571

    Lorentzian(double y_optimum, double width_parameter) :
        ExposureWeight(y_optimum, width_parameter * FWHM_GAUSSIAN / FWHM_LORENTZIAN) {}

    double weight(double y) override
    {
        const double z = normalize(y);
        return 1.0 / (1.0 + 0.5 * z * z);
    }
};


struct HalfSinusodial : public ExposureWeight
{
    // FWHM = 2 * arccos(1/2)
#define FWHM_HALFSINUSODIAL 2.0943951023931954923084289221863352561

    HalfSinusodial(double y_optimum, double width_parameter) :
        ExposureWeight(y_optimum, width_parameter * FWHM_GAUSSIAN / FWHM_HALFSINUSODIAL) {}

    double weight(double y) override
    {
        const double z = normalize(y);
        return fabs(z) <= M_PI_2 ? cos(z) : 0.0;
    }
};


struct FullSinusodial : public ExposureWeight
{
    // FWHM = pi
#define FWHM_FULLSINUSODIAL M_PI

    FullSinusodial(double y_optimum, double width_parameter) :
        ExposureWeight(y_optimum, width_parameter * FWHM_GAUSSIAN / FWHM_FULLSINUSODIAL) {}

    double weight(double y) override
    {
        const double z = normalize(y);
        return fabs(z) <= M_PI ? 0.5 + cos(z) / 2.0 : 0.0;
    }
};


struct Bisquare : public ExposureWeight
{
    // FWHM = 2 / sqrt(sqrt(2))
#define FWHM_BISQUARE 1.6817928305074290860622509524664297901

    Bisquare(double y_optimum, double width_parameter) :
        ExposureWeight(y_optimum, width_parameter * FWHM_GAUSSIAN / FWHM_BISQUARE) {}

    double weight(double y) override
    {
        const double z = normalize(y);

        if (fabs(z) <= 1.0)
        {
            const double z2 = z * z;
            return (1.0 - z2) * (1.0 + z2);
        }
        else
        {
            return 0.0;
        }
    }
};


#endif // EXPOSURE_WEIGHT_INCLUDED


// Local Variables:
// mode: c++
// End:
