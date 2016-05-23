/*
 * Copyright (C) 2012-2016 Dr. Christoph L. Spiel
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

#ifndef MINIMIZER_H_INCLUDED
#define MINIMIZER_H_INCLUDED


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <stdexcept>
#include <vector>

#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "optional_transitional.hpp"


class Minimizer
{
    enum {ITERATIONS_PER_DIMENSION = 100U};

public:
    Minimizer() = delete;
    explicit Minimizer(size_t a_dimension);
    Minimizer(const Minimizer& a_minimizer);
    Minimizer& operator=(const Minimizer& a_minimizer);
    virtual ~Minimizer() {}

    virtual std::string proper_name() const = 0;

    size_t dimension() const {return dimension_;}
    void set_dimension(size_t a_dimension) {dimension_ = a_dimension;}

    Minimizer* set_maximum_number_of_iterations(unsigned n);
    Minimizer* unset_maximum_number_of_iterations();
    unsigned number_of_iterations() const {return iteration_;}
    void next_iteration() {++iteration_;}

    Minimizer* set_goal(double a_goal);
    Minimizer* unset_goal();

    Minimizer* set_absolute_error(double absolute_error);
    Minimizer* unset_absolute_error();

    virtual double f_minimum() const = 0;

    virtual bool has_reached_goal() const;
    virtual bool has_reached_maximum_iteration() const;

    struct minimum_not_bracketed : public std::runtime_error
    {
        double minimum;
        double lower;
        double upper;

        explicit minimum_not_bracketed(const std::string& a_message) :
            std::runtime_error(a_message),
            minimum(0.0), lower(0.0), upper(0.0)
        {}

        minimum_not_bracketed(const std::string& a_message, double x_minimum, double x_lower, double x_upper) :
            std::runtime_error(a_message + " [(" +
                               std::to_string(x_lower) + " | " +
                               std::to_string(x_minimum) + " | " +
                               std::to_string(x_upper) + ")]"),
            minimum(x_minimum), lower(x_lower), upper(x_upper)
        {}
    };

protected:
    virtual double absolute_error() const;

private:
    size_t dimension_;
    std::optional<unsigned> maximum_iteration_;
    unsigned iteration_;
    std::optional<double> f_goal_;
    std::optional<double> absolute_error_;
};


class Minimizer1D : public Minimizer
{
public:
    Minimizer1D() = delete;
    Minimizer1D(const gsl_function& a_function, double x_minimum, double x_lower, double x_upper);
    Minimizer1D(const Minimizer1D& a_minimizer);
    Minimizer1D& operator=(const Minimizer1D& a_minimizer);
    virtual ~Minimizer1D() {gsl_min_fminimizer_free(minimizer_);}

    virtual std::string proper_name() const;

    void set_bracket(const gsl_function& a_function, double x_minimum, double x_lower, double x_upper);

    Minimizer1D* set_relative_error(double a_relative_error);
    Minimizer1D* unset_relative_error();

    void run();

    virtual double x_minimum() const {return x_minimum_;}
    virtual double f_minimum() const {return gsl_min_fminimizer_f_minimum(minimizer_);}

protected:
    const gsl_min_fminimizer_type* type() const {return type_;}

    void require_ordered_x() const;
    void initialize(const gsl_min_fminimizer_type* a_minimizer_type);

    virtual double relative_error() const;
    virtual bool has_reached_tolerance() const;

private:
    const gsl_min_fminimizer_type* type_;
    gsl_min_fminimizer* minimizer_;

    gsl_function function_;

    double x_minimum_;
    double x_lower_;
    double x_upper_;

    std::optional<double> relative_error_;
};


class GoldenSectionMinimizer1D : public Minimizer1D
{
public:
    GoldenSectionMinimizer1D() = delete;
    GoldenSectionMinimizer1D(const gsl_function& a_function, double x_minimum, double x_lower, double x_upper);
    GoldenSectionMinimizer1D(const GoldenSectionMinimizer1D& a_minimizer);
};


class BrentMinimizer1D : public Minimizer1D
{
public:
    BrentMinimizer1D() = delete;
    BrentMinimizer1D(const gsl_function& a_function, double x_minimum, double x_lower, double x_upper);
    BrentMinimizer1D(const BrentMinimizer1D& a_minimizer);
};


class GillMurrayMinimizer1D : public Minimizer1D
{
public:
    GillMurrayMinimizer1D() = delete;
    GillMurrayMinimizer1D(const gsl_function& a_function, double x_minimum, double x_lower, double x_upper);
    GillMurrayMinimizer1D(const GillMurrayMinimizer1D& a_minimizer);
};


////////////////////////////////////////////////////////////////////////


template <class input_iterator>
inline static void
copy_to_gsl_vector(input_iterator first, input_iterator last, gsl_vector* a_vector)
{
    unsigned i = 0U;

    while (first != last)
    {
        gsl_vector_set(a_vector, i, *first);
        ++i;
        ++first;
    }
}


template <class output_iterator>
inline static output_iterator
copy_from_gsl_vector(gsl_vector* a_vector, output_iterator a_result)
{
    for (unsigned i = 0U; i != a_vector->size; ++i)
    {
        *a_result = gsl_vector_get(a_vector, i);
        ++a_result;
    }

    return a_result;
}


////////////////////////////////////////////////////////////////////////


class MinimizerMultiDimensionNoDerivative : public Minimizer
{
public:
    typedef std::vector<double> array_type;

    MinimizerMultiDimensionNoDerivative() = delete;
    MinimizerMultiDimensionNoDerivative(const gsl_multimin_function& a_function,
                                        const array_type& a_start, const array_type& some_step_sizes);
    MinimizerMultiDimensionNoDerivative(const gsl_multimin_function& a_function, const array_type& a_start);
    MinimizerMultiDimensionNoDerivative(const MinimizerMultiDimensionNoDerivative& a_minimizer);
    MinimizerMultiDimensionNoDerivative& operator=(const MinimizerMultiDimensionNoDerivative& a_minimizer);
    virtual ~MinimizerMultiDimensionNoDerivative();

    void set_start(const array_type& a_start);
    void set_step_sizes(const array_type& some_step_sizes);

    template <class output_iterator>
    output_iterator get_step_sizes(output_iterator a_result) const
    {
        return copy_from_gsl_vector(step_sizes_, a_result);
    }

    virtual std::string proper_name() const;

    void run();

    template <class output_iterator>
    output_iterator x_minimum(output_iterator a_result) const
    {
        return copy_from_gsl_vector(gsl_multimin_fminimizer_x(minimizer_), a_result);
    }

    virtual double f_minimum() const;
    double characteristic_size() const {return characteristic_size_;}

protected:
    const gsl_multimin_fminimizer_type* type() const {return type_;}

    void initialize_step_sizes(const array_type& some_step_sizes);
    void set();
    void initialize(const gsl_multimin_fminimizer_type* a_fminimizer_type);

private:
    const gsl_multimin_fminimizer_type* type_;
    gsl_multimin_fminimizer* minimizer_;

    gsl_multimin_function function_;
    gsl_vector* xs_;
    gsl_vector* step_sizes_;

    double characteristic_size_;
};


class MinimizerMultiDimensionSimplex : public MinimizerMultiDimensionNoDerivative
{
public:
    MinimizerMultiDimensionSimplex() = delete;
    MinimizerMultiDimensionSimplex(const gsl_multimin_function& a_function,
                                   const array_type& a_start, const array_type& some_step_sizes);
    MinimizerMultiDimensionSimplex(const gsl_multimin_function& a_function, const array_type& a_start);
    MinimizerMultiDimensionSimplex(const MinimizerMultiDimensionSimplex& a_minimizer);
};


class MinimizerMultiDimensionSimplex2 : public MinimizerMultiDimensionNoDerivative
{
public:
    MinimizerMultiDimensionSimplex2() = delete;

    MinimizerMultiDimensionSimplex2(const gsl_multimin_function& a_function,
                                    const array_type& a_start, const array_type& some_step_sizes);
    MinimizerMultiDimensionSimplex2(const gsl_multimin_function& a_function, const array_type& a_start);
    MinimizerMultiDimensionSimplex2(const MinimizerMultiDimensionSimplex2& a_minimizer);
};


class MinimizerMultiDimensionSimplex2Randomized : public MinimizerMultiDimensionNoDerivative
{
public:
    MinimizerMultiDimensionSimplex2Randomized() = delete;
    MinimizerMultiDimensionSimplex2Randomized(const gsl_multimin_function& a_function,
                                              const array_type& a_start, const array_type& some_step_sizes);
    MinimizerMultiDimensionSimplex2Randomized(const gsl_multimin_function& a_function,
                                              const array_type& a_start);
    MinimizerMultiDimensionSimplex2Randomized(const MinimizerMultiDimensionSimplex2Randomized& a_minimizer);
};


#endif // MINIMIZER_H_INCLUDED

// Local Variables:
// mode: c++
// End:
