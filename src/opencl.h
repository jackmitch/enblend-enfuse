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
#include <sstream>
#include <string>
#include <vector>

#if defined(HAVE_CL_CL_HPP)
#include <CL/cl.hpp>
#elif defined(HAVE_OPENCL_CL_HPP)
#include <OpenCL/cl.hpp>
#endif


namespace ocl
{
#if defined(_OPENCL) || defined(__OPENCL_CL_HPP)

#define OPENCL

    typedef std::vector<cl::Platform> platform_list_t;
    typedef std::vector<cl::Device> device_list_t;


    std::string
    string_of_error_code(cl_int error_code)
    {
        switch (error_code)
        {
        case CL_SUCCESS: return "success"; // not an error

        case CL_DEVICE_NOT_FOUND: return "device not found";
        case CL_DEVICE_NOT_AVAILABLE: return "device not available";
        case CL_COMPILER_NOT_AVAILABLE: return "compiler not available";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE: return "memory object allocation failure";
        case CL_OUT_OF_RESOURCES: return "out of resources";
        case CL_OUT_OF_HOST_MEMORY: return "out of host memory";
        case CL_PROFILING_INFO_NOT_AVAILABLE: return "profiling information not available";
        case CL_MEM_COPY_OVERLAP: return "memory copy overlap";
        case CL_IMAGE_FORMAT_MISMATCH: return "image format mismatch";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED: return "image format not supported";
        case CL_BUILD_PROGRAM_FAILURE: return "build program failure";
        case CL_MAP_FAILURE: return "map failure";
        case CL_MISALIGNED_SUB_BUFFER_OFFSET: return "misaligned sub buffer offset";
        case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST: return "exec status error for events in wait list";

        case CL_INVALID_VALUE: return "invalid value";
        case CL_INVALID_DEVICE_TYPE: return "invalid device type";
        case CL_INVALID_PLATFORM: return "invalid platform";
        case CL_INVALID_DEVICE: return "invalid device";
        case CL_INVALID_CONTEXT: return "invalid context";
        case CL_INVALID_QUEUE_PROPERTIES: return "invalid queue properties";
        case CL_INVALID_COMMAND_QUEUE: return "invalid command queue";
        case CL_INVALID_HOST_PTR: return "invalid host pointer";
        case CL_INVALID_MEM_OBJECT: return "invalid memory object";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR: return "invalid image format descriptor";
        case CL_INVALID_IMAGE_SIZE: return "invalid image size";
        case CL_INVALID_SAMPLER: return "invalid sampler";
        case CL_INVALID_BINARY: return "invalid binary";
        case CL_INVALID_BUILD_OPTIONS: return "invalid build options";
        case CL_INVALID_PROGRAM: return "invalid program";
        case CL_INVALID_PROGRAM_EXECUTABLE: return "invalid program executable";
        case CL_INVALID_KERNEL_NAME: return "invalid kernel name";
        case CL_INVALID_KERNEL_DEFINITION: return "invalid kernel definition";
        case CL_INVALID_KERNEL: return "invalid kernel";
        case CL_INVALID_ARG_INDEX: return "invalid argument index";
        case CL_INVALID_ARG_VALUE: return "invalid argument value";
        case CL_INVALID_ARG_SIZE: return "invalid argument size";
        case CL_INVALID_KERNEL_ARGS: return "invalid kernel arguments";
        case CL_INVALID_WORK_DIMENSION: return "invalid work dimension";
        case CL_INVALID_WORK_GROUP_SIZE: return "invalid work group size";
        case CL_INVALID_WORK_ITEM_SIZE: return "invalid work item size";
        case CL_INVALID_GLOBAL_OFFSET: return "invalid global offset";
        case CL_INVALID_EVENT_WAIT_LIST: return "invalid event wait list";
        case CL_INVALID_EVENT: return "invalid event";
        case CL_INVALID_OPERATION: return "invalid operation";
        case CL_INVALID_GL_OBJECT: return "invalid GL object";
        case CL_INVALID_BUFFER_SIZE: return "invalid buffer size";
        case CL_INVALID_MIP_LEVEL: return "invalid MIP level";
        case CL_INVALID_GLOBAL_WORK_SIZE: return "invalid global work size";
        case CL_INVALID_PROPERTY: return "invalid property";

        default:
            std::ostringstream oss;
            oss << "unknown error code " << error_code;
            return oss.str();
        }
    }


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
            unsigned platform_index = 1U; // We start enumerating at 1 for user convenience.
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
                    unsigned device_index = 1U; // We start enumerating at 1 for user convenience.
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
} // namespace ocl


#endif // OPENCL_H_INCLUDED

// Local Variables:
// mode: c++
// End:
