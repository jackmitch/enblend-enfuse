/*
 * Copyright (C) 2013 Christoph L. Spiel
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
#ifndef OPENCL_H_INCLUDED
#define OPENCL_H_INCLUDED


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>

#define __CL_ENABLE_EXCEPTIONS

#if defined(HAVE_CL_CL_HPP)
#include <CL/cl.hpp>
#elif defined(HAVE_OPENCL_CL_HPP)
#include <OpenCL/cl.hpp>
#endif


namespace ocl
{
#if defined(_OPENCL) || defined(CL_HPP_)

#define OPENCL

    typedef std::vector<cl::Platform> platform_list_t;
    typedef std::vector<cl::Device> device_list_t;


    class runtime_error : public std::runtime_error
    {
    public:
        runtime_error() = delete;
        explicit runtime_error(const std::string& a_message) : std::runtime_error(a_message) {};
        virtual ~runtime_error() throw() {}
    }; // class runtime_error


    std::string string_of_error_code(cl_int error_code);

    void print_opencl_information(bool all_devices = false);

    cl::Platform find_platform(/* input/output */ size_t& a_preferred_platform_id);
    void prefer_device(const cl::Platform& a_platform,
                       size_t a_preferred_platform_id,
                       size_t a_preferred_device_id,
                       /* output */ device_list_t& some_devices);

    // Create a new context given a_platform and some_devices.  This
    // pointer can be deleted as usual.
    cl::Context* create_context(const cl::Platform& a_platform, const device_list_t& some_devices);

#else

#endif // _OPENCL
} // namespace ocl


#endif // OPENCL_H_INCLUDED

// Local Variables:
// mode: c++
// End:
