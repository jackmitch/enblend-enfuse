/*
 * Copyright (C) 2015 Christoph L. Spiel
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
#include <iostream>

#include "global.h"
#include "openmp_def.h"         // omp::atomic_t

#include "dynamic_loader.h"     // HAVE_DYNAMICLOADER_IMPL
#include "opencl.h"             // macro OPENCL
#include "opencl_exposure_weight.h"

#include "exposure_weight.h"



extern const std::string command;
extern ExposureWeight* ExposureWeightFunction;



namespace exposure_weight
{
#ifdef HAVE_DYNAMICLOADER_IMPL
    class DynamicExposureWeight : public ExposureWeight
    {
    public:
        DynamicExposureWeight() = delete;

        DynamicExposureWeight(const std::string& library_name, const std::string& symbol_name,
                              double y_optimum = 0.5, double width = 0.2) :
            ExposureWeight(y_optimum, width),
            library_(library_name), symbol_(symbol_name),
            dynamic_loader_(DynamicLoader(library_name)),
            function_(dynamic_loader_.resolve<ExposureWeight*>(symbol_name))
        {
            assert(function_);
            if (function_->interface_version() != EXPOSURE_WEIGHT_INTERFACE_VERSION)
            {
                std::cerr <<
                    command << ": user-defined weight function \"" << symbol_name << "\"\n" <<
                    command << ": defined in shared object \"" << library_name << "\"\n" <<
                    command << ": matches interface version " << function_->interface_version() <<
                    ", but " << command << " requires version " << EXPOSURE_WEIGHT_INTERFACE_VERSION <<
                    std::endl;
                exit(1);
            }
        }

        void initialize(double y_optimum, double width_parameter,
                        ExposureWeight::argument_const_iterator user_arguments_begin,
                        ExposureWeight::argument_const_iterator user_arguments_end) override
        {
#ifdef DEBUG
            std::cout << "+ DynamicExposureWeight::initialize\n";
            std::for_each(user_arguments_begin, user_arguments_end,
                          [](const std::string& x)
                          {std::cout << "+ DynamicExposureWeight::initialize: <" << x << ">\n";});
#endif
            function_->initialize(y_optimum, width_parameter, user_arguments_begin, user_arguments_end);
        }

        double weight(double y) override {return function_->weight(y);}

    private:
        std::string library_;
        std::string symbol_;
        DynamicLoader dynamic_loader_;
        ExposureWeight* function_;
    };


    static ExposureWeight*
    make_dynamic_weight_function(const std::string& name,
                                 ExposureWeight::argument_const_iterator arguments_begin,
                                 ExposureWeight::argument_const_iterator arguments_end,
                                 double y_optimum, double width)
    {
        if (arguments_begin == arguments_end)
        {
            // Remember that built-in exposure-weight functions never
            // take any arguments.
            std::cerr <<
                command << ": unknown built-in exposure weight function \"" << name << "\"" << std::endl;
            exit(1);
        }
        else
        {
            const std::string symbol_name = *arguments_begin;
            ExposureWeight* weight_object;

            try
            {
                weight_object = new DynamicExposureWeight(name, symbol_name);
                weight_object->initialize(y_optimum, width, std::next(arguments_begin), arguments_end);
            }
            catch (ExposureWeight::error& exception)
            {
                std::cerr <<
                    command << ": user-defined weight function \"" << symbol_name << "\"\n" <<
                    command << ": defined in shared object \"" << name << "\"\n" <<
                    command << ": raised exception: " << exception.what() << std::endl;
                exit(1);
            }

            return weight_object;
        }
    }
#endif // HAVE_DYNAMICLOADER_IMPL


    ExposureWeight*
    make_weight_function(const std::string& name,
                         ExposureWeight::argument_const_iterator arguments_begin,
                         ExposureWeight::argument_const_iterator arguments_end,
                         double y_optimum, double width)
    {
        delete ExposureWeightFunction;

        std::string possible_built_in(name);
        enblend::to_lower(possible_built_in);

        if (possible_built_in == "gauss" || possible_built_in == "gaussian")
        {
            return new Gaussian(y_optimum, width);
        }
        else if (possible_built_in == "lorentz" || possible_built_in == "lorentzian")
        {
            return new Lorentzian(y_optimum, width);
        }
        else if (possible_built_in == "halfsine" || possible_built_in == "half-sine")
        {
            return new HalfSinusodial(y_optimum, width);
        }
        else if (possible_built_in == "fullsine" || possible_built_in == "full-sine")
        {
            return new FullSinusodial(y_optimum, width);
        }
        else if (possible_built_in == "bisquare" || possible_built_in == "bi-square")
        {
            return new Bisquare(y_optimum, width);
        }
        else
        {
#if defined(HAVE_DYNAMICLOADER_IMPL) && defined(OPENCL)
#ifdef DEBUG
            std::cerr << "+ make_weight_function: HAVE_DYNAMICLOADER_IMPL && OPENCL" << std::endl;
#endif
            if (opencl_exposure_weight::is_opencl_file(name))
            {
                return opencl_exposure_weight::make_weight_function(name,
                                                                    arguments_begin, arguments_end,
                                                                    y_optimum, width);
            }
            else
            {
                return make_dynamic_weight_function(name,
                                                    arguments_begin, arguments_end,
                                                    y_optimum, width);
            }
#elif defined(OPENCL)
#ifdef DEBUG
            std::cerr << "+ make_weight_function: OPENCL only" << std::endl;
#endif
            if (opencl_exposure_weight::is_opencl_file(name))
            {
                return opencl_exposure_weight::make_weight_function(name,
                                                                    arguments_begin, arguments_end,
                                                                    y_optimum, width);
            }
            else
            {
                std::cerr << command << ": OpenCL source file required" << << std::endl;
                exit(1);
            }
#elif defined(HAVE_DYNAMICLOADER_IMPL)
#ifdef DEBUG
            std::cerr << "+ make_weight_function: HAVE_DYNAMICLOADER_IMPL only" << std::endl;
#endif
            if (opencl_exposure_weight::is_opencl_file(name))
            {
                std::cerr << command << ": shared-object file (aka dynamic library) required" << std::endl;
                exit(1);
            }
            else
            {
                return make_dynamic_weight_function(name,
                                                    arguments_begin, arguments_end,
                                                    y_optimum, width);
            }
#else
            std::cerr <<
                command << ": unknown built-in exposure weight function \"" << name << "\"\n" <<
                command << ": note: this binary has no support for dynamic loading of\n" <<
                command << ": note: exposure weight functions" << std::endl;
            exit(1);
#endif
        }
    }


    void
    dump_weight_function(ExposureWeight* weight_function, int n)
    {
        assert(n >= 2);

        for (int i = 0; i < n; ++i)
        {
            const double x = static_cast<double>(i) / static_cast<double>(n - 1);
            const double w = weight_function->weight(x);

            std::cout << i << ' ' << x << ' ' << w << '\n';
        }
    }


    bool
    check_weight_function(ExposureWeight* weight_function, int n)
    {
        assert(n >= 2);

        omp::atomic_t number_of_faults = omp::atomic_t();

#ifdef OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; ++i)
        {
            const double y = static_cast<double>(i) / static_cast<double>(n - 1);
            const double w = weight_function->weight(y);

            if (w < 0.0 || w >= 1.0)
            {
                ++number_of_faults;
            }
        }

        return number_of_faults == omp::atomic_t();
    }
} // namespace exposure_weight
