/*
 * Copyright (C) 2009-2014 Dr. Christoph L. Spiel
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


// No test, no bug.  -- Chris Spiel


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <getopt.h>
#include <stdlib.h>
#include <string.h>

#include <cassert>
#include <cmath>                // fabsf
#include <iostream>
#include <string>

#include "self_test.h"


#define lengthof(m_array) (sizeof(m_array) / sizeof(m_array[0]))


extern const std::string command;


////////////////////////////////////////////////////////////////////////


// Number of arguments we pass to getopt_long in each of our tests.
#define ARG_COUNT 4


enum {a_short, b_short, c_short, a_long, b_long, c_long, FLAG_COUNT};


struct test_case {
    const char* arguments[ARG_COUNT];
    int flags[FLAG_COUNT];
};


inline static int
int_of_string(const char* s)
{
    return static_cast<int>(strtol(s, nullptr, 10));
}


static void
reset_getopt_globals()
{
    opterr = 0;                 // silence getopt_long(3)
    optopt = -1;                // reset "unknown option" character
    optind = 1;                 // reset parsing index
    optarg = nullptr;              // reset pointer to value of option argument
}


static int
try_out_getopt_long(int arg_count, const char* arguments[], int* flags)
{
    const char* short_options = "ab:c::";
    const struct option long_options[] = {
        {"long-a", no_argument,       nullptr, 1},
        {"long-b", required_argument, nullptr, 2},
        {"long-c", optional_argument, nullptr, 3},
        {nullptr, 0, nullptr, 0}
    };

    reset_getopt_globals();
    while (true)
    {
        int option_index = 0;
        int code = getopt_long(arg_count, const_cast<char* const*>(arguments),
                               short_options, long_options,
                               &option_index);

        if (code == -1)
        {
            break;
        }

        switch (code)
        {
        case 1:
            flags[a_long] = 1;
            break;
        case 2:
            flags[b_long] = int_of_string(optarg);
            break;
        case 3:
            flags[c_long] = optarg == nullptr ? 1 : int_of_string(optarg);
            break;
        case 'a':
            flags[a_short] = 1;
            break;
        case 'b':
            flags[b_short] = int_of_string(optarg);
            break;
        case 'c':
            flags[c_short] = optarg == nullptr ? 1 : int_of_string(optarg);
            break;
        default:
            return -1;
        }
    }

    return optind;
}


// Write a list of elements separated by spaces to stream out.
template <typename T>
static void
write_list(std::ostream& out, unsigned size, const T list)
{
    for (unsigned i = 0U; i != size; ++i)
    {
        out << list[i];
        if (i != size - 1U)
        {
            out << ' ';
        }
    }
}


// Name of the first argument, i.e. the first non-option in the list
// of arguments.  We need to know its name so that we can check
// whether getopt_long(3) really parsed all options.
#define ARG1 "1"


// Test whether the library function getopt_long(3) works as required.
bool
getopt_long_works_ok()
{
    bool has_passed_test = true;
    struct test_case tests[] = {
        {{"p", ARG1, "2", "3"},                    {0, 0, 0, 0, 0, 0}},

        {{"p", "-a", ARG1, "2"},                   {1, 0, 0, 0, 0, 0}},
        {{"p", "-b2", ARG1, "2"},                  {0, 2, 0, 0, 0, 0}},
        {{"p", "-c", ARG1, "2"},                   {0, 0, 1, 0, 0, 0}},
        {{"p", "-c2", ARG1, "2"},                  {0, 0, 2, 0, 0, 0}},

        {{"p", "--long-a", ARG1, "2"},             {0, 0, 0, 1, 0, 0}},
        {{"p", "--long-b=2", ARG1, "2"},           {0, 0, 0, 0, 2, 0}},
        {{"p", "--long-c", ARG1, "2"},             {0, 0, 0, 0, 0, 1}},
        {{"p", "--long-c=2", ARG1, "2"},           {0, 0, 0, 0, 0, 2}},

        {{"p", "-a", "-b2", ARG1},                 {1, 2, 0, 0, 0, 0}},
        {{"p", "-a", "-b2", ARG1},                 {1, 2, 0, 0, 0, 0}},
        {{"p", "-ab2", "-c", ARG1},                {1, 2, 1, 0, 0, 0}},
        {{"p", "-ab2", "-c3", ARG1},               {1, 2, 3, 0, 0, 0}},

        {{"p", "--long-a", "--long-b=2", ARG1},    {0, 0, 0, 1, 2, 0}},
        {{"p", "--long-a", "--long-b=2", ARG1},    {0, 0, 0, 1, 2, 0}},
        {{"p", "--long-b=2", "--long-c", ARG1},    {0, 0, 0, 0, 2, 1}},
        {{"p", "--long-b=2", "--long-c=3", ARG1},  {0, 0, 0, 0, 2, 3}},

        {{"p", "-a", "--long-a", ARG1},            {1, 0, 0, 1, 0, 0}},
        {{"p", "-b2", "--long-a", ARG1},           {0, 2, 0, 1, 0, 0}},
        {{"p", "-a", "--long-b=2", ARG1},          {1, 0, 0, 0, 2, 0}},
        {{"p", "-b2", "--long-b=2", ARG1},         {0, 2, 0, 0, 2, 0}},

        {{nullptr, nullptr, nullptr}, {0, 0, 0, 0, 0, 0}}
    };
    const unsigned arg_count = lengthof(tests->arguments);
    const unsigned flag_count = lengthof(tests->flags);

    for (struct test_case* t = tests; t->arguments[0] != nullptr; ++t)
    {
        int flags[] = {0, 0, 0, 0, 0, 0};
        assert(lengthof(tests->flags) == lengthof(flags));
        const int index = try_out_getopt_long(arg_count, t->arguments, flags);

        if (index < 0 || index >= static_cast<int>(arg_count) ||
            strcmp(t->arguments[index], ARG1) != 0)
        {
            std::cerr <<
                command <<
                ": failed self test: getopt_long(3) did not parse argument list \"";
            write_list(std::cerr, arg_count, t->arguments);
            std::cerr << "\"\n";

            has_passed_test = false;
        }

        for (unsigned i = 0U; i != flag_count; ++i)
        {
            if (flags[i] != t->flags[i])
            {
                std::cerr <<
                    command <<
                    ": failed self test: getopt_long(3) incorrectly parses argument list \"";
                write_list(std::cerr, arg_count, t->arguments);
                std::cerr << "\";\n";

                std::cerr <<
                    command <<
                    ": failed self test:     expected {";
                write_list(std::cerr, flag_count, t->flags);
                std::cerr << "}, but got {";
                write_list(std::cerr, flag_count, flags);
                std::cerr << "}\n";

                has_passed_test = false;
            }
        }
    }

    reset_getopt_globals();

    return has_passed_test;
}


// Run a kernel, if we have OpenCL support.
#ifdef OPENCL

typedef std::vector<float> float_vector;


static const std::string
axpy_source("kernel void\n"
            "axpy(const float alpha,\n"
            "     global const float *restrict x,\n"
            "     global const float *restrict y,\n"
            "     const int n,\n"
            "     global float *restrict z)\n"
            "{\n"
            "    const int i = get_global_id(0);\n"
            "\n"
            "    if (i >= n)\n"
            "    {\n"
            "        return;\n"
            "    }\n"
            "\n"
            "    z[i] = alpha * x[i] + y[i];\n"
            "}\n");


class Alpha_times_x_plus_y : public ocl::BuildableFunction
{
public:
    Alpha_times_x_plus_y() = delete;
    explicit Alpha_times_x_plus_y(const cl::Context& a_context) : f_(a_context, axpy_source) {}

    void build(const std::string& a_build_option)
    {
        std::cerr << "\n+ Alpha_times_x_plus_y::build: by request\n\n";
        f_.build(a_build_option);
        std::cerr <<
            "+ Alpha_times_x_plus_y::build: log begin ================\n" <<
            f_.build_log() <<
            "\n+ Alpha_times_x_plus_y::build: log end   ================\n";
    }

    void wait()
    {
        f_.wait();
        std::cerr << "\n+ Alpha_times_x_plus_y::wait: carry on...\n\n";
        initialize();
    }

    void run(float alpha, const float_vector& x, const float_vector& y, float_vector& z)
    {
        const size_t n = x.size();
        assert(n == y.size());
        assert(n <= z.size());

        const size_t buffer_size = n * sizeof(float);

        cl::Buffer x_buffer(f_.context(), CL_MEM_READ_ONLY, buffer_size);
        cl::Buffer y_buffer(f_.context(), CL_MEM_READ_ONLY, buffer_size);
        cl::Buffer z_buffer(f_.context(), CL_MEM_WRITE_ONLY, buffer_size);

        kernel_.setArg(0U, static_cast<cl_float>(alpha));
        kernel_.setArg(1U, x_buffer);
        kernel_.setArg(2U, y_buffer);
        kernel_.setArg(3U, static_cast<cl_int>(n));
        kernel_.setArg(4U, z_buffer);

        f_.queue().enqueueWriteBuffer(x_buffer, CL_TRUE, 0U, buffer_size, &x[0]);
        f_.queue().enqueueWriteBuffer(y_buffer, CL_TRUE, 0U, buffer_size, &y[0]);
        f_.queue().enqueueNDRangeKernel(kernel_, cl::NullRange, cl::NDRange(n), cl::NullRange);
        f_.queue().enqueueReadBuffer(z_buffer, CL_TRUE, 0U, buffer_size, &z[0]);
    }

private:
    void initialize()
    {
        kernel_ = f_.create_kernel("axpy");
    }

    ocl::LazyFunctionCXXOfString f_;
    cl::Kernel kernel_;
}; // Alpha_times_x_plus_y


static void
alpha_times_x_plus_y(float alpha, const float_vector& x, const float_vector& y, float_vector& z)
{
    const size_t n = x.size();

    assert(n == y.size());
    assert(n <= z.size());

    for (size_t i = 0U; i != n; ++i)
    {
        z[i] = alpha * x[i] + y[i];
    }
}


template <class forward_iterator, class t>
inline static void
iota(forward_iterator first, forward_iterator last, t a_value)
{
    while (first != last)
    {
        *first++ = a_value++;
    }
}


static bool
test_axpy_on_gpu(cl::Context* a_context)
{
    const float alpha = 2.5F;
    const size_t n = 20U;

    float_vector x(n);
    float_vector y(n);
    float_vector z0(n);
    float_vector z(n);

    iota(x.begin(), x.end(), 1.0F);
    iota(y.begin(), y.end(), 1.0F);

    alpha_times_x_plus_y(alpha, x, y, z0);

    try
    {
        Alpha_times_x_plus_y axpy(*a_context);

        axpy.build("-Werror");
        axpy.wait();
        axpy.run(alpha, x, y, z);
    }
    catch (cl::Error& a_cl_error)
    {
        std::cerr <<
            command << ": warning: plain cl error: " << a_cl_error.what() << "\n" <<
            command << ": warning:     reason: " << ocl::string_of_error_code(a_cl_error.err()) <<
            std::endl;
        return false;
    }
    catch (ocl::runtime_error& an_opencl_runtime_error)
    {
        std::cerr <<
            command << ": warning: ocl error: " << an_opencl_runtime_error.what() << "\n" <<
            command << ": warning:     reason: " <<
            ocl::string_of_error_code(an_opencl_runtime_error.error().err()) << "\n" <<
            command << ": warning:     message: " << an_opencl_runtime_error.additional_message() <<
            std::endl;
        return false;
    }
    catch (...)
    {
        std::cerr <<
            command << ": warning: unknown exception thrown during self test \"test_axpy_on_gpu\"" <<
            std::endl;
        return false;
    }

    for (size_t i = 0U; i != n; ++i)
    {
#ifdef DEBUG
        std::cout <<
            "+ test_axpy_on_gpu: [" << i << "]  " <<
            alpha << " * " << x[i] << " + " << y[i] << " = " << z[i] << " (reference: " << z0[i] << "), " <<
            "delta: " << std::scientific << fabsf(z[i] - z0[i]) << std::fixed << "\n";
#endif // DEBUG
        if (fabsf(z[i] - z0[i]) > std::numeric_limits<float>::epsilon())
        {
            std::cerr <<
                "+ test_axpy_on_gpu: failure at index " << i <<
                ", expected " << z0[i] << ", but got " << z[i] << "\n";
            return false;
        }
    }

    return true;
}


bool
gpu_is_ok(cl::Context* a_context)
{
    return test_axpy_on_gpu(a_context);
}

#endif // OPENCL
