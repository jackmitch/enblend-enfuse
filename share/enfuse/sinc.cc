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
#include <array>                  // std::array<>
#include <cassert>
#include <cmath>                  // std::abs(), std::sqrt()
#include <functional>             // std::unary_function<>
#include <iostream>               // std::cerr
#include <limits>                 // std::numeric_limits<>
#include <sstream>                // std::istringstream

#include "exposure_weight_base.h" // macro FWHM_GAUSSIAN, class ExposureWeight
#include "root_finding.hh"        // root_finding::dekker()


// sinus_cardinalis(x) = if(x == 0, 1, sin(x) / x)
inline static double
sinus_cardinalis(double x)
{
    return x == 0.0 ? 1.0 : std::sin(x) / x;
}


// scaled_sinc(x) = sinus_cardinalis(Pi * x)
inline static double
scaled_sinc(double x)
{
    return sinus_cardinalis(M_PI * x);
}


// truncated_sinc(x) = if(abs(x) >= 1, 0, scaled_sinc(x))
inline static double
truncated_sinc(double x)
{
    return std::fabs(x) >= 1.0 ? 0.0 : scaled_sinc(x);
}


class Function : public std::unary_function<double, double>
{
public:
    virtual double operator()(double) const = 0;
    virtual double fwhm() const = 0;
    virtual ~Function() {}

protected:
    double find_fwhm(double x0, double x1) const
    {
        assert(x0 >= 0.0);
        const double y_half = root_finding::dekker(x0, x1, [&](double x) {return operator()(x) - 0.5;});
        assert(std::fabs(operator()(y_half) - 0.5) < std::sqrt(std::numeric_limits<double>::epsilon()));

        return 2.0 * y_half;
    }
};


class Sinc : public Function
{
public:
    double operator()(double x) const {return truncated_sinc(x);}
    double fwhm() const {return 1.2067091288032283952915076867831260818;}
};


class PowerSinc : public Function
{
public:
    PowerSinc() : argument_exponent(1.0), exponent(1.0) {}

    PowerSinc(double a, double e) : argument_exponent(a), exponent(e)
    {
        assert(argument_exponent >= 1.0);
        assert(exponent >= 0.0);
    }

    double operator()(double x) const
    {
        return std::pow(truncated_sinc(std::pow(x, argument_exponent)), exponent);
    }

    double fwhm() const {return find_fwhm(0.0, 1.0);}

private:
    const double argument_exponent;
    const double exponent;
};


////////////////////////////////////////////////////////////////////////////////


class SincWeight : public ExposureWeight
{
public:
    void initialize(double y_optimum, double width, const argument_list_t& argument_list) override
    {
        if (!argument_list.empty())
        {
            std::cerr << "warning: weight function does not take any parameters" << std::endl;
        }

        ExposureWeight::initialize(y_optimum, width * FWHM_GAUSSIAN / sinc.fwhm(), argument_list);
    }

    double weight(double y) override
    {
        return sinc(std::abs(normalize(y)));
    }

private:
    Sinc sinc;
};


class PowerSincWeight : public ExposureWeight
{
    typedef ExposureWeight super;

public:
    PowerSincWeight() : sinc(nullptr) {}
    PowerSincWeight(const PowerSincWeight&) = delete;
    PowerSincWeight& operator=(const PowerSincWeight&) = delete;
    ~PowerSincWeight() {delete sinc;}


    void initialize(double y_optimum, double width, const argument_list_t& argument_list) override
    {
        std::array<double, 2> exponents = {1.0, 1.0};

        for (size_t i = 0U; i != exponents.size(); ++i)
        {
            if (i >= argument_list.size())
            {
                break;
            }

            std::istringstream iss(argument_list.at(i));

            if (!(iss >> exponents.at(i)))
            {
                throw super::error("non-numeric exponent");
            }
        }

        sinc = new PowerSinc(exponents.at(0), exponents.at(1));

        ExposureWeight::initialize(y_optimum, width * FWHM_GAUSSIAN / sinc->fwhm(), argument_list);
    }


    double weight(double y) override
    {
        return sinc->operator()(std::abs(normalize(y)));
    }

private:
    PowerSinc* sinc;
};


////////////////////////////////////////////////////////////////////////////////


SincWeight sinc;
PowerSincWeight power_sinc;
