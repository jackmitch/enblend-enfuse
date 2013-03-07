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

#include <iostream>
#include <vector>

#if defined(HAVE_CL_CL_HPP)
#include <CL/cl.hpp>
#elif defined(HAVE_OPENCL_CL_HPP)
#include <OpenCL/cl.hpp>
#endif


#if defined(_OPENCL) || defined(__OPENCL_CL_HPP)

#define OPENCL

typedef std::vector<cl::Platform> platform_list_t;
typedef std::vector<cl::Device> device_list_t;


void
print_opencl_information()
{
    platform_list_t platforms;
    cl::Platform::get(&platforms);

    if (platforms.empty())
    {
        std::cout << "  - no platform found\n";
    }
    else
    {
        unsigned platform_index = 0U;
        for (platform_list_t::const_iterator p = platforms.begin(); p != platforms.end(); ++p, ++platform_index)
        {
            std::string info;
            p->getInfo(CL_PLATFORM_VENDOR, &info);
            std::cout << "  - Platform #" << platform_index << ": " << info;
            p->getInfo(CL_PLATFORM_NAME, &info);
            std::cout << ", " << info;
            p->getInfo(CL_PLATFORM_VERSION, &info);
            std::cout << ", " << info << "\n";

            device_list_t devices;
            p->getDevices(CL_DEVICE_TYPE_GPU, &devices);

            if (devices.empty())
            {
                std::cout << "  - no device found\n";
            }
            else
            {
                unsigned device_index = 0U;
                for (device_list_t::const_iterator d = devices.begin(); d != devices.end(); ++d, ++device_index)
                {
                    std::cout << "    * Device #" << device_index << ": " <<
                        d->getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << " cores\n" <<
                        "                 " << d->getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() / 1024UL << " KB global memory\n" <<
                        "                 " << d->getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() / 1024UL << " KB local memory\n" <<
                        "                 " << d->getInfo<CL_DEVICE_IMAGE2D_MAX_WIDTH>() << 'x' <<
                        d->getInfo<CL_DEVICE_IMAGE2D_MAX_HEIGHT>() << " maximum image size\n";
                }
            }
        }
    }
}

#else

#endif // _OPENCL


#endif // OPENCL_H_INCLUDED

// Local Variables:
// mode: c++
// End:
