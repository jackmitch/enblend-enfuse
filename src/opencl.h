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


    std::string string_of_error_code(cl_int error_code);

    void print_opencl_information(bool all_devices = false);

#else

#endif // _OPENCL
} // namespace ocl


#endif // OPENCL_H_INCLUDED

// Local Variables:
// mode: c++
// End:
