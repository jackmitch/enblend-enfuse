/*
 * Copyright (C) 2014 Christoph L. Spiel
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


#include <cassert>

#include <time.h>

#include "openmp_def.h"

#include "mersenne.h"


MersenneTwister::MersenneTwister() :
    generator_(gsl_rng_alloc(gsl_rng_mt19937))
{
    assert(generator_);
}


MersenneTwister::MersenneTwister(const MersenneTwister& another_generator) :
    generator_(gsl_rng_clone(another_generator.generator_))
{
    assert(generator_);
}


MersenneTwister::~MersenneTwister()
{
    gsl_rng_free(generator_);
}


MersenneTwister&
MersenneTwister::operator=(const MersenneTwister& another_generator)
{
    if (this != &another_generator)
    {
        gsl_rng_free(generator_);
        generator_ = gsl_rng_clone(another_generator.generator_);
        assert(generator_);
    }
    return *this;
}


void
MersenneTwister::seed()
{
    gsl_rng_set(generator_, gsl_rng_default_seed);
}


void
MersenneTwister::seed(result_type a_seed)
{
    gsl_rng_set(generator_, a_seed);
}


void
UniformMersenneTwister::non_deterministic_seed()
{
#ifdef DEBUG
    seed(); // apply default seed in all threads for reproducibility
#else
    unsigned seed_value = static_cast<unsigned>(1 + omp_get_thread_num());

    // ANTICIPATED CHANGE: Tap into /dev/random.
    const clock_t now = clock();
    if (now != static_cast<clock_t>(-1))
    {
        seed_value ^= static_cast<unsigned>(now);
    }

    seed(seed_value);
#endif
}
