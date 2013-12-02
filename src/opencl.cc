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


#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "opencl.h"


namespace ocl
{
#if defined(OPENCL)

    std::string
    string_of_error_code(cl_int an_error_code)
    {
        switch (an_error_code)
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
            std::ostringstream error_code;

            error_code << "unknown error code " << an_error_code;

            return error_code.str();
        }
    }


    static void
    print_platform_info(platform_list_t::const_iterator a_platform, unsigned a_platform_index)
    {
        std::string info;

        a_platform->getInfo(CL_PLATFORM_VENDOR, &info);
        std::cout << "  - Platform #" << a_platform_index << ": " << info;
        a_platform->getInfo(CL_PLATFORM_NAME, &info);
        std::cout << ", " << info;
        a_platform->getInfo(CL_PLATFORM_VERSION, &info);
        std::cout << ", " << info << "\n";
    }


    static void
    print_device_info(device_list_t::const_iterator a_device, unsigned a_device_index)
    {
        std::cout <<
            "    * Device #" << a_device_index << ": max. " <<
            a_device->getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << " work-items\n" <<
            "                 " << a_device->getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() / 1024UL <<
            " KB global memory ";

        switch (a_device->getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_TYPE>())
        {
        case CL_NONE: std::cout << "without associated cache";  break;
        case CL_READ_ONLY_CACHE:
            std::cout <<
                "with " <<
                a_device->getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>() / 1024UL <<
                " KB read cache";
            break;
        case CL_READ_WRITE_CACHE:
            std::cout <<
                "with " <<
                a_device->getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>() / 1024UL <<
                " KB read/write cache";
        }

        std::cout << "\n" <<
            "                 " << a_device->getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() / 1024UL << " KB " <<
            (a_device->getInfo<CL_DEVICE_LOCAL_MEM_TYPE>() == CL_LOCAL ? "dedicated " : "") <<
            "local memory\n" <<
            "                 " << a_device->getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>() / 1024UL <<
            " KB maximum constant memory\n";
    }


    struct no_platform : runtime_error
    {
        no_platform() : runtime_error("no OpenCL platform found") {}
        explicit no_platform(const std::string& a_message) : runtime_error(a_message) {}
    };


    struct no_device : runtime_error
    {
        no_device() : runtime_error("no OpenCL device found") {}
        explicit no_device(const std::string& a_message) : runtime_error(a_message) {}
    };


    void
    print_opencl_information(bool all_devices)
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
            for (platform_list_t::const_iterator p = platforms.begin(); p != platforms.end();
                 ++p, ++platform_index)
            {
                print_platform_info(p, platform_index);

                device_list_t devices;
                try
                {
                    p->getDevices(all_devices ? CL_DEVICE_TYPE_ALL : CL_DEVICE_TYPE_GPU, &devices);
                }
                catch (cl::Error& an_error)
                {
                    // CL_DEVICE_NOT_FOUND is a possible error that
                    // does not hurt as the variable devices will be
                    // empty() then, which is checked below.
                    if (an_error.err() != CL_DEVICE_NOT_FOUND)
                    {
                        throw an_error;
                    }
                }

                if (devices.empty())
                {
                    std::cout <<
                        "    * no " << (all_devices ? "" : "GPU ") <<
                        "devices found on this platform\n";
                }
                else
                {
                    unsigned device_index = 1U; // Again, we start enumerating at 1 for user convenience.
                    for (device_list_t::const_iterator d = devices.begin(); d != devices.end();
                         ++d, ++device_index)
                    {
                        print_device_info(d, device_index);
                    }
                }
            }
        }
    }


    void
    print_gpu_preference(size_t a_preferred_platform_id, size_t a_preferred_device_id)
    {
        try
        {
            size_t platform_id = a_preferred_platform_id;
            cl::Platform platform;
            device_list_t some_devices;

            platform = find_platform(platform_id);
            prefer_device(platform, a_preferred_platform_id, a_preferred_device_id, some_devices);

            std::cout <<
                "Currently preferred GPU is device #" << a_preferred_device_id <<
                " on platform #" << platform_id <<
                (a_preferred_platform_id == 0U ? " (autodetected)" : "") << ".\n";
        }
        catch (no_platform&)
        {
            std::cout << "No OpenCL platforms found.\n";
        }
        catch (no_device&)
        {
            std::cout << "No OpenCL (GPU) devices found on any platform.\n";
        }
        catch (runtime_error& an_error)
        {
            std::cout <<
                "Platform number #" << a_preferred_platform_id <<
                (a_preferred_platform_id == 0U ? " (autodetected)" : "") <<
                "/device number #" <<
                a_preferred_device_id << " combination is invalid for this system.\n" <<
                an_error.what() << "\n";
        }
    }


    cl::Platform
    find_platform(size_t& a_preferred_platform_id)
    {
        std::ostringstream message;

        platform_list_t platforms;
        try
        {
            cl::Platform::get(&platforms);
        }
        catch (cl::Error& an_error)
        {
            message << "query for OpenCL platforms failed: " << ocl::string_of_error_code(an_error.err());
            throw runtime_error(message.str());
        }

        if (platforms.empty())
        {
            throw no_platform();
        }
        else
        {
            if (a_preferred_platform_id == 0U)
            {
                ocl::platform_list_t::const_iterator p =
                    std::find_if(platforms.begin(), platforms.end(),
                                 [](const cl::Platform& a_platform)
                                 {
                                     ocl::device_list_t devices;
                                     try {a_platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);}
                                     catch (cl::Error&) {return false;}
                                     return devices.size() >= 1U;
                                 });
                if (p == platforms.end())
                {
                    throw no_device();
                }
                else
                {
                    a_preferred_platform_id = p - platforms.begin() + 1U;
                    return *p;
                }
            }
            else if (a_preferred_platform_id <= platforms.size())
            {
                return platforms[a_preferred_platform_id - 1U];
            }
            else
            {
                message <<
                    "OpenCL platform #" << a_preferred_platform_id <<
                    " is not available; largest OpenCL platform number is " << platforms.size();
                throw runtime_error(message.str());
            }
        }
    }


    void
    prefer_device(const cl::Platform& a_platform, size_t a_preferred_platform_id,
                  size_t a_preferred_device_id, device_list_t& some_devices)
    {
        std::ostringstream message;

        try
        {
            a_platform.getDevices(CL_DEVICE_TYPE_GPU, &some_devices);
        }
        catch (cl::Error& an_error)
        {
            message <<
                "query for OpenCL GPU devices on platform #" << a_preferred_platform_id + 1U << " failed: " <<
                ocl::string_of_error_code(an_error.err());
            throw runtime_error(message.str());
        }

        if (some_devices.empty())
        {
            message << "no OpenCL GPU device found on platform #" << a_preferred_platform_id;
            throw no_device(message.str());
        }
        else
        {
            if (a_preferred_device_id <= some_devices.size())
            {
                // move the preferred device in front
                some_devices.insert(some_devices.begin(), some_devices[a_preferred_device_id - 1U]);
                some_devices.erase(some_devices.begin() + a_preferred_device_id);
            }
            else
            {
                message <<
                    "OpenCL device #" << a_preferred_device_id <<
                    " is not available on platform #" << a_preferred_platform_id + 1U <<
                    ", largest device number there is " << some_devices.size();
                throw runtime_error(message.str());
            }
        }
    }


    static void
    run_self_tests(cl::Context* a_context)
    {
        std::ostringstream message;

        // a_context must be usable.
        std::vector<cl_context_properties> context_properties;
        try
        {
            a_context->getInfo(CL_CONTEXT_PROPERTIES, &context_properties);
        }
        catch (cl::Error& an_error)
        {
            message <<
                "self test failed: cannot query properties of context: " <<
                ocl::string_of_error_code(an_error.err());
            throw runtime_error(message.str());
        }

        // We need at least one device.
        std::vector<cl::Device> devices;
        try
        {
            a_context->getInfo(CL_CONTEXT_DEVICES, &devices);
        }
        catch (cl::Error& an_error)
        {
            message << "self test failed: cannot query devices in context: " <<
                ocl::string_of_error_code(an_error.err());
            throw runtime_error(message.str());
        }

        if (devices.empty())
        {
            throw no_device();
        }
    }


    cl::Context*
    create_context(const cl::Platform& a_platform, const device_list_t& some_devices)
    {
        cl_context_properties context_properties[] = {
            CL_CONTEXT_PLATFORM,
            (cl_context_properties) (a_platform)(),
            0
        };
        cl::Context* context = nullptr;

        try
        {
            context = new cl::Context(some_devices, context_properties, nullptr, nullptr);
        }
        catch (cl::Error& an_error)
        {
            std::ostringstream message;
            message << "failed to create OpenCL context: " << ocl::string_of_error_code(an_error.err());
            throw runtime_error(message.str());
        }

        run_self_tests(context);

        return context;
    }

#else

#endif // OPENCL
} // namespace ocl
