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


#ifndef INTERPOLATOR_HH_INCLUDED
#define INTERPOLATOR_HH_INCLUDED


#include <vector>

#include <gsl/gsl_interp.h>


class Interpolator
{
public:
    typedef enum interpolator_type {Linear, CubicSpline, Akima} interpolator_t;

    Interpolator();
    Interpolator(const std::vector<double>& some_xs, const std::vector<double>& some_ys,
                 interpolator_t a_type = Linear);
    Interpolator(const Interpolator& another_interpolator);
    Interpolator& operator=(const Interpolator& another_interpolator);
    virtual ~Interpolator();

    void set_data(const std::vector<double>& some_xs, const std::vector<double>& some_ys);
    void set_type(interpolator_t a_type);

    virtual double evaluate(double x); // Non-const as method may change member `accelerator'!

private:
    static const gsl_interp_type* map_type(enum interpolator_type a_type);

    std::vector<double> xs;
    std::vector<double> ys;

    const gsl_interp_type* type;      // nothing to allocate, points to GSL data
    gsl_interp* interpolator;
    gsl_interp_accel* accelerator;
};


#endif // INTERPOLATOR_HH_INCLUDED
