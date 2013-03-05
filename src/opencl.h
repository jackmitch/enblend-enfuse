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

#if defined(HAVE_CL_OPENCL_H)
#include <CL/opencl.h>
#elif defined(HAVE_OPENCL_OPENCL_H)
#include <OpenCL/opencl.h>
#endif


#if defined(_OPENCL) || defined(__OPENCL_CL_H)

#define OPENCL

void
print_opencl_information()
{
    cl_uint number_of_platforms = 0U;

    if (clGetPlatformIDs(0U, NULL, &number_of_platforms) == CL_SUCCESS)
    {
        cl_platform_id *platform_ids = new cl_platform_id[number_of_platforms];
        clGetPlatformIDs(number_of_platforms, platform_ids, NULL);

        const size_t maximum_info_size = 256U;
        char* info = new char[maximum_info_size];
        size_t size;

        for (cl_uint p = 0U; p != number_of_platforms; ++p)
        {
            clGetPlatformInfo(platform_ids[p], CL_PLATFORM_VENDOR, maximum_info_size, info, &size);
            std::cout << "  - " << info;
            clGetPlatformInfo(platform_ids[p], CL_PLATFORM_NAME, maximum_info_size, info, &size);
            std::cout << ", " << info;
            clGetPlatformInfo(platform_ids[p], CL_PLATFORM_VERSION, maximum_info_size, info, &size);
            std::cout << ", " << info << "\n";

            cl_uint number_of_devices = 0U;

            if (clGetDeviceIDs(platform_ids[p], CL_DEVICE_TYPE_GPU, 0U, NULL, &number_of_devices) == CL_SUCCESS)
            {
                cl_device_id *device_ids = new cl_device_id [number_of_devices];
                clGetDeviceIDs(platform_ids[p], CL_DEVICE_TYPE_GPU, number_of_devices, device_ids, NULL);

                for (cl_uint d = 0U; d != number_of_devices; ++d)
                {
                    cl_uint number_of_cores;
                    std::cout << "    * ";
                    if (clGetDeviceInfo(device_ids[d], CL_DEVICE_MAX_COMPUTE_UNITS,
                                        sizeof(cl_uint), &number_of_cores, NULL) == CL_SUCCESS)
                    {
                        std::cout << number_of_cores;
                    }
                    else
                    {
                        std::cout << "???";
                    }
                    std::cout << " cores\n";

                    cl_ulong memory_size;
                    std::cout << "    * ";
                    if (clGetDeviceInfo(device_ids[d], CL_DEVICE_GLOBAL_MEM_SIZE,
                                        sizeof(cl_ulong), &memory_size, NULL) == CL_SUCCESS)
                    {
                        std::cout << memory_size / 1024UL;
                    }
                    else
                    {
                        std::cout << "???";
                    }
                    std::cout << " KB global memory\n";

                    std::cout << "    * ";
                    if (clGetDeviceInfo(device_ids[d], CL_DEVICE_LOCAL_MEM_SIZE,
                                        sizeof(cl_ulong), &memory_size, NULL) == CL_SUCCESS)
                    {
                        std::cout << memory_size / 1024UL;
                    }
                    else
                    {
                        std::cout << "???";
                    }
                    std::cout << " KB local memory\n";

                    size_t max_width;
                    size_t max_height;
                    if (clGetDeviceInfo(device_ids[d], CL_DEVICE_IMAGE2D_MAX_WIDTH,
                                        sizeof(size_t), &max_width, NULL) == CL_SUCCESS &&
                        clGetDeviceInfo(device_ids[d], CL_DEVICE_IMAGE2D_MAX_HEIGHT,
                                        sizeof(size_t), &max_height, NULL) == CL_SUCCESS)
                    {
                        std::cout << "    * " << max_width << 'x' << max_height;
                    }
                    else
                    {
                        std::cout << "unknown";
                    }
                    std::cout << " maximum image size\n";
                }

                delete [] device_ids;
            }
            else
            {
                std::cout << "  - no GPU devices found\n";
            }
        }

        delete [] info;
        delete [] platform_ids;
    }
    else
    {
        std::cout << "  - no platform found\n";
    }
}

#else

#endif // _OPENCL


#endif // OPENCL_H_INCLUDED

// Local Variables:
// mode: c++
// End:
