// Copyright (C) 2014-2015 Christoph L. Spiel
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
#include <cassert>                // macro assert()
#include <cctype>                 // std::isspace()
#include <fstream>                // std::ifstream
#include <iostream>               // std::cerr
#include <sstream>                // std::istringstream
#include <string>                 // std::string
#include <unordered_map>          // std::unordered_map
#include <vector>                 // std::vector<>

#include "exposure_weight_base.h" // macro FWHM_GAUSSIAN, class ExposureWeight
#include "interpolator.hh"        // class Interpolator


typedef std::pair<std::vector<double>, std::vector<double> > value_table;


static void
dump_value_table(value_table a_value_table)
{
    const value_table::first_type& ys = a_value_table.first;
    const value_table::second_type& ws = a_value_table.second;

    assert(ys.size() == ws.size());

    for (size_t i = 0U; i != ys.size(); ++i)
    {
        std::cout << '[' << i << "]    " << ys[i] << '\t' << ws[i] << '\n';
    }
}


// Answer whether a_string consists *only* of syntactic white-space
// which either is plain white-space or white-space trailed by a
// comment character ('#') and an arbitrary string after in until the
// end of a_string.
static bool
is_syntactic_whitespace(const std::string& a_string)
{
    for (auto c : a_string)
    {
        if (c == '#')
        {
            return true;
        }
        if (!isspace(c))
        {
            return false;
        }
    }

    return true;
}


static value_table
read_data_file(const std::string& a_datafile_name)
{
    std::ifstream data_file(a_datafile_name.c_str());
    unsigned line_number = 1U;
    std::string line;

    std::vector<double> xs;
    std::vector<double> ys;

    while (std::getline(data_file, line))
    {
        if (!is_syntactic_whitespace(line))
        {
            std::istringstream iss(line);
            double x;
            double y;

            if (!(iss >> x >> y))
            {
                std::ostringstream error_message;

                error_message <<
                    "data file: \"" << a_datafile_name << "\", line: " << line_number <<
                    " - failed to parse (luminance, weight) pair";

                std::cerr << error_message.str() << std::endl;
            }
            else
            {
                xs.push_back(x);
                ys.push_back(y);
            }
        }

        ++line_number;
    }

    return std::make_pair(xs, ys);
}


static value_table
read_data_arguments(ExposureWeight::argument_const_iterator arguments_begin,
                    ExposureWeight::argument_const_iterator arguments_end)
{
    int argument_number = 1;
    ExposureWeight::argument_const_iterator argument = arguments_begin;

    std::vector<double> xs;
    std::vector<double> ys;

    while (argument != arguments_end) // for (const auto& a : argument_list)
    {
        double x;
        double y;
        char delimiter;

        if (sscanf(argument->c_str(), "%lf%lf", &x, &y) == 2 ||
            (sscanf(argument->c_str(), "%lf%c%lf", &x, &delimiter, &y) == 3 && delimiter == '/'))
        {
            xs.push_back(x);
            ys.push_back(y);
        }
        else
        {
            std::ostringstream error_message;

            error_message <<
                "immediate data #" << argument_number <<
                " - failed to parse (luminance, weight) pair";

            std::cerr << error_message.str() << std::endl;
        }

        ++argument;
        ++argument_number;
    }

    return std::make_pair(xs, ys);
}


class Tabular : public ExposureWeight
{
    typedef ExposureWeight super;

public:
    ~Tabular() {delete interpolator;}

    void initialize(double y_optimum, double width_parameter,
                    super::argument_const_iterator arguments_begin,
                    super::argument_const_iterator arguments_end) override
    {
        if (y_optimum != 0.5)
        {
            std::cerr <<
                "Tabular::initialize: warning: ignoring parameter \"--exposure-optimum\"" <<
                std::endl;
        }
        if (width_parameter != 0.2)
        {
            std::cerr <<
                "Tabular::initialize: warning: ignoring parameter \"--exposure-width\"" <<
                std::endl;
        }

        if (arguments_begin == arguments_end)
        {
            throw super::error("missing data-source parameter");
        }

        const std::string& data_source_tag = *arguments_begin;

        if (data_source_tag == "immediate")
        {
            super::argument_const_iterator first_data_point = arguments_begin + 1;

            if (first_data_point == arguments_end)
            {
                throw super::error("no data supplied");
            }

            value_table data {read_data_arguments(first_data_point, arguments_end)};
            dump_value_table(data);
            interpolator = new Interpolator(data.first, data.second);
        }
        else if (data_source_tag == "file")
        {
            super::argument_const_iterator filename = arguments_begin + 1;

            if (filename == arguments_end)
            {
                throw super::error("missing data-file name parameter");
            }

            value_table data {read_data_file(*filename)};
            dump_value_table(data);
            interpolator = new Interpolator(data.first, data.second);
        }
        else
        {
            throw super::error("unrecognized data-source");
        }
    }


    double weight(double y) override
    {
        // Implementation Note: Guard against messing up
        // `gsl_interp_accel' inside of Interpolator with a critical
        // section.
#ifdef _OPENMP
#pragma omp critical tabular_evaluate_weight
#endif
        {
            // Implementation Note: We play it safe by using std::max() as
            // an interpolated function may overshoot to the negative
            // side.
            return std::max(interpolator->evaluate(y), 0.0);
        }
    }

private:
    Interpolator* interpolator;
}; // Tabular


class HashedTabular : public Tabular
{
    typedef Tabular super;

public:
    void initialize(double y_optimum, double width_parameter,
                    super::argument_const_iterator arguments_begin,
                    super::argument_const_iterator arguments_end) override
    {
        super::initialize(y_optimum, width_parameter, arguments_begin, arguments_end);
    }

    double weight(double y) override
    {
        return super::weight(y);
    }

private:
    std::unordered_map<double, double> cache;
}; // HashedTabular


Tabular tabular;
HashedTabular htabular;
