/*
 * Copyright (C) 2015  Christoph L. Spiel
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

#ifndef MERSENNE_H_INCLUDED_
#define MERSENNE_H_INCLUDED_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gsl/gsl_rng.h>


class MersenneTwister
{
public:
    typedef unsigned long result_type;

    MersenneTwister();
    MersenneTwister(const MersenneTwister& another_generator);
    virtual ~MersenneTwister();

    MersenneTwister& operator=(const MersenneTwister& another_generator);

    result_type min() const {return gsl_rng_min(generator_);}
    result_type max() const {return gsl_rng_max(generator_);}

    void seed();
    void seed(result_type a_seed);

    result_type operator()() {return gsl_rng_get(generator_);}

private:
    gsl_rng* generator_;
}; // class MersenneTwister


class UniformMersenneTwister : public MersenneTwister
{
public:
    void non_deterministic_seed();

    result_type get() {return this->operator()();}
    double get_uniform() {return static_cast<double>(get()) / static_cast<double>(max());}
}; // class UniformMersenneTwister


#endif // MERSENNE_H_INCLUDED_

// Local Variables:
// mode: c++
// End:
