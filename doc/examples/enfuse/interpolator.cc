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


#include <cassert>

#include "interpolator.hh"


Interpolator::Interpolator() :
    type(nullptr), interpolator(nullptr), accelerator(nullptr)
{}


Interpolator::Interpolator(const std::vector<double>& some_xs, const std::vector<double>& some_ys,
                           interpolator_t a_type) :
    xs(some_xs),
    ys(some_ys),
    type(map_type(a_type)),
    interpolator(gsl_interp_alloc(type, xs.size())),
    accelerator(gsl_interp_accel_alloc())
{
    const size_t n {xs.size()};

    assert(n == ys.size() && "number of ordinates and abscissas differs");
    assert(n >= gsl_interp_type_min_size(type) &&
           "not enough data points for chosen type of interpolation");

    // sort pairs (x, y) according to x values.
    gsl_interp_init(interpolator, xs.data(), ys.data(), n);
}


Interpolator::Interpolator(const Interpolator& another_interpolator) :
    xs(another_interpolator.xs),
    ys(another_interpolator.ys),
    type(another_interpolator.type),
    interpolator(gsl_interp_alloc(type, xs.size())),
    accelerator(gsl_interp_accel_alloc())
{}


Interpolator&
Interpolator::operator=(const Interpolator& another_interpolator)
{
    if (this != &another_interpolator)
    {
        xs = another_interpolator.xs;
        ys = another_interpolator.ys;

        type = another_interpolator.type;

        if (interpolator)
        {
            gsl_interp_free(interpolator);
        }
        interpolator = gsl_interp_alloc(type, xs.size());

        if (accelerator)
        {
            gsl_interp_accel_free(accelerator);
        }
        accelerator = gsl_interp_accel_alloc();
    }

    return *this;
}


Interpolator::~Interpolator()
{
    gsl_interp_accel_free(accelerator);
    gsl_interp_free(interpolator);
}


void
Interpolator::set_data(const std::vector<double>& some_xs, const std::vector<double>& some_ys)
{
    assert(some_xs.size() != some_ys.size()  && "number of ordinates and abscissas differs");

    xs = some_xs;
    ys = some_ys;

    if (interpolator)
    {
        gsl_interp_free(interpolator);
    }
    interpolator = gsl_interp_alloc(type, xs.size());

    if (accelerator)
    {
        gsl_interp_accel_free(accelerator);
    }
    accelerator = gsl_interp_accel_alloc();
}


void
Interpolator::set_type(interpolator_t a_type)
{
    const gsl_interp_type* new_type {map_type(a_type)};

    if (type != new_type)
    {
        type = new_type;

        if (!xs.empty())
        {
            assert(xs.size() >= gsl_interp_type_min_size(type) &&
                   "not enough data points for chosen type of interpolation");

            if (interpolator)
            {
                gsl_interp_free(interpolator);
            }
            interpolator = gsl_interp_alloc(type, xs.size());
        }
    }
}


double
Interpolator::evaluate(double x)
{
    return gsl_interp_eval(interpolator, xs.data(), ys.data(), x, accelerator);
}


// https://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html#Interpolation-Types

const gsl_interp_type*
Interpolator::map_type(interpolator_t a_type)
{
    switch (a_type)
    {
    case CubicSpline:
        return gsl_interp_cspline;

    case Akima:
        return gsl_interp_akima;

    case Linear: // FALLTHROUGH
    default:
        return gsl_interp_linear;
    }
}
