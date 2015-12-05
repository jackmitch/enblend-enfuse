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


#include <algorithm>            // std::for_each
#include <cassert>
#include <cmath>                // std::lround
#include <ctime>                // std::time_t, std::time()
#include <fstream>              // std::ifstream
#include <iostream>
#include <string>

#include "opencl.h"
#include "parameter.h"

#include "opencl_exposure_weight.h"


extern const std::string command;

#ifdef OPENCL
extern cl::Context* GPUContext;
#endif


namespace opencl_exposure_weight
{
#ifdef OPENCL
    class UserWeightFunction : public ocl::BuildableFunction
    {
        static const std::string header_code;
        static const std::string kernel_code;

    public:
        UserWeightFunction() = delete;

        explicit UserWeightFunction(const cl::Context& a_context) :
            context_(a_context), f_(nullptr),
            preferred_work_group_size_multiple_(0U), work_group_size_(0U)
        {}

        virtual ~UserWeightFunction()
        {
            delete f_;
        }

        std::string user_code;

        void build(const std::string& a_build_option)
        {
            delete f_;
            f_ = new ocl::LazyFunctionCXXOfString(context_, header_code + user_code + kernel_code);

            f_->add_build_option("-cl-single-precision-constant");
            f_->add_build_option("-DENFUSE_FWHM_GAUSSIAN=%.6ff", FWHM_GAUSSIAN);

            switch (f_->vendor_id())
            {
            case ocl::vendor::amd:
                f_->add_build_option("-g");
                break;
            case ocl::vendor::nvidia:
                f_->add_build_option("-cl-nv-verbose");
                break;
            case ocl::vendor::unknown:
                break;
            }

            try
            {
#ifdef DEBUG_OPENCL
                std::cout << "\n+ UserWeightFunction::build: by request\n\n";
#endif

                f_->build(a_build_option);

#ifdef DEBUG_OPENCL
                std::cout <<
                    "+ UserWeightFunction::build: log begin ================\n" <<
                    f_->build_log() <<
                    "\n+ UserWeightFunction::build: log end   ================\n" <<
                    std::endl;
#endif
            }
            catch (ocl::runtime_error& an_error)
            {
                std::cerr << command << ": " << ocl::string_of_error_code(an_error.error().err()) << "\n";

                std::vector<std::string> messages =
                    ocl::split_string(an_error.additional_message(), '\n', true);
                for (auto m : messages)
                {
                    std::cerr << command << ": note: " << m << "\n";
                }

                exit(1);
            }

            initialize();
        }

        void wait()
        {
            f_->wait();

#ifdef DEBUG_OPENCL
            std::cerr << "\n+ UserWeightFunction::wait: carry on...\n\n";
#endif

            initialize();
        }

        void initialize()
        {
            if (work_group_size_ > 0U) // initialize only once
            {
                return;
            }

            weight_kernel_ = f_->create_kernel("enfuse_user_weight_kernel");

            weight_kernel_.getWorkGroupInfo(f_->device(),
                                            CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                            &preferred_work_group_size_multiple_);

            size_t device_max_group_size;
            f_->device().getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &device_max_group_size);
            size_t kernel_group_size;
            weight_kernel_.getWorkGroupInfo(f_->device(), CL_KERNEL_WORK_GROUP_SIZE, &kernel_group_size);
            work_group_size_ = std::min(device_max_group_size, kernel_group_size);
        }

        size_t work_group_size() const
        {
            assert(work_group_size_ > 0U);
            return work_group_size_;
        }

        size_t work_group_size(size_t a_suggested_work_group_size) const
        {
            assert(work_group_size_ > 0U);
            return ocl::round_up_to_next_multiple(a_suggested_work_group_size,
                                                  preferred_work_group_size_multiple_);
        }

        template <typename random_iterator, typename output_iterator>
        void evaluate(random_iterator luminances_begin, random_iterator luminances_end,
                      output_iterator some_weights)
        {
            const size_t size = luminances_end - luminances_begin;
            const size_t raw_size = size * sizeof(float);

            cl::Buffer luminance_buffer(context_, CL_MEM_READ_ONLY, raw_size);
            cl::Buffer weight_buffer(context_, CL_MEM_READ_WRITE, raw_size);

            const size_t group_size = work_group_size(32U);
            const size_t m = size % group_size;
            const size_t n = size / group_size;

#ifdef DEBUG_OPENCL
            std::cout <<
                "+ UserWeightFunction::evaluate: size = " << size << " -> raw-size = " << raw_size << "\n" <<
                "+ UserWeightFunction::evaluate: work_group_size_ = " << work_group_size_ << "\n" <<
                "+ UserWeightFunction::evaluate: group_size = " << group_size << "\n" <<
                "+ UserWeightFunction::evaluate: n = " << n << ", m = " << m <<
                std::endl;
#endif

            weight_kernel_.setArg(0U, static_cast<cl_int>(size));
            weight_kernel_.setArg(1U, luminance_buffer);
            weight_kernel_.setArg(2U, weight_buffer);

            cl::CommandQueue& queue = f_->queue();

            {
                ocl::ScopedWriteMap luminance_map(queue, luminance_buffer, CL_TRUE, 0U, raw_size);

                float* ys = luminance_map.get<float*>();
                for (random_iterator lumi = luminances_begin; lumi != luminances_end; ++lumi)
                {
                    *ys++ = static_cast<float>(*lumi);
                }
            }

            std::vector<cl::Event> done_with_kernels(n);
            for (size_t i = 0U; i != n; ++i)
            {
                queue.enqueueNDRangeKernel(weight_kernel_,
                                           cl::NDRange(i * group_size),
                                           cl::NDRange(group_size),
                                           cl::NullRange,
                                           nullptr,
                                           &done_with_kernels[i]);
                DEBUG_CHECK_OPENCL_EVENT(done_with_kernels[i]);
            }
            cl::Event::waitForEvents(done_with_kernels);

            {
                ocl::ScopedReadMap weight_map(queue, weight_buffer, CL_TRUE, 0U, raw_size);

                const float* const weights = weight_map.get<float*>();
                for (size_t j = 0U; j != size; ++j)
                {
                    *some_weights++ = static_cast<double>(weights[j]);
                }
            }
        }

    private:
        cl::Context context_;
        ocl::LazyFunctionCXXOfString* f_;
        cl::Kernel weight_kernel_;
        size_t preferred_work_group_size_multiple_;
        size_t work_group_size_;
    }; // class UserWeightFunction


    const std::string
    UserWeightFunction::header_code =
        "#line 1 \"ENFUSE-INTERNAL-HEADER\"\n"
        "float\n"
        "enfuse_normalized_luminance(float y)\n"
        "{\n"
        "    return (y - ENFUSE_OPTIMUM_Y) / ENFUSE_WIDTH;\n"
        "}\n"
        "\n"
        "\n"
        "float ENFUSE_USER_WEIGHT_FUNCTION(float);\n"
        "\n"
        "\n";


    const std::string
    UserWeightFunction::kernel_code =
        "\n"
        "\n"
        "#line 1 \"ENFUSE-INTERNAL-KERNEL\"\n"
        "kernel void\n"
        "enfuse_user_weight_kernel(const int n,\n"
        "                          global const float *y,\n"
        "                          global       float *w)\n"
        "{\n"
        "    const int i = get_global_id(0);\n"
        "\n"
        "    if (i < n)\n"
        "    {\n"
        "        w[i] = ENFUSE_USER_WEIGHT_FUNCTION(y[i]);\n"
        "    }\n"
        "}\n"
        "\n"
        "\n";


    class OpenCLUserExposureWeight : public ExposureWeight
    {
        enum {OPENCL_USER_WEIGHT_SAMPLES = 32U};

    public:
        OpenCLUserExposureWeight() = delete;

        OpenCLUserExposureWeight(const std::string& source_file_name,
                                 const std::string& weight_function_name) :
            ExposureWeight(0.5, 0.25),
            source_file_name_(source_file_name),
            weight_function_name_(weight_function_name),
            user_weight_function_(ocl::create_function<UserWeightFunction>(GPUContext)),
            number_of_samples_(parameter::as_unsigned("opencl-user-weight-samples", OPENCL_USER_WEIGHT_SAMPLES))
        {}

        OpenCLUserExposureWeight(const std::string& source_file_name,
                                 const std::string& weight_function_name,
                                 double y_optimum, double width) :
            ExposureWeight(y_optimum, width),
            source_file_name_(source_file_name),
            weight_function_name_(weight_function_name),
            user_weight_function_(ocl::create_function<UserWeightFunction>(GPUContext)),
            number_of_samples_(parameter::as_unsigned("opencl-user-weight-samples", OPENCL_USER_WEIGHT_SAMPLES))
        {}

        void initialize(double y_optimum, double width_parameter,
                        ExposureWeight::argument_const_iterator arguments_begin,
                        ExposureWeight::argument_const_iterator arguments_end) override
        {
            std::ostringstream build_options;

            for (auto p : ocl::construct_search_path())
            {
#ifdef _WIN32
#define PATHS_DELIMITER ';'
#else
#define PATHS_DELIMITER ':'
#endif
                for (auto d : ocl::split_string(p, PATHS_DELIMITER))
                {
                    build_options << "-I" << d << " ";
                }
            }

            build_options <<
                "-DENFUSE_USER_WEIGHT_FUNCTION=" << weight_function_name_ << " " <<
                "-DENFUSE_OPTIMUM_Y=" << y_optimum << "f " <<
                "-DENFUSE_WIDTH=" << width_parameter << "f ";

            std::for_each(arguments_begin, arguments_end,
                          [&] (const std::string& x) {build_options << "-D" << x << " ";});

            {
                std::ostringstream user_code;
                if (parameter::as_boolean("consult-opencl-user-exposure-weight-file", true))
                {
                    user_code <<
                        "#line 1 \"" << source_file_name_ << "\"\n" <<
                        ocl::consult_file(source_file_name_);
                }
                else
                {
                    // Implementation Note: We put an ever-changing
                    // time-stamp into the source code to trick the
                    // OpenCL build-system into re-compiling.
                    //
                    // Usually the build-system is dependency
                    // agnostic, i.e. it will not recompile if the
                    // contents of `source_file_name_' changes.

                    const std::time_t now(std::time(nullptr));

                    user_code <<
                        "#define ENFUSE_TIME_STAMP " << now << "\n" <<
                        "#include \"" << source_file_name_ << "\"\n";

                }
                user_weight_function_->user_code = user_code.str();
            }

#ifdef DEBUG_OPENCL
            std::cout <<
                "+ OpenCLUserExposureWeight::initialize\n" <<
                "+     source_file_name_ = <" << source_file_name_ << ">\n" <<
                "+     weight_function_name_ = <" << weight_function_name_ << ">\n" <<
                "+     build-options = <" << build_options.str() << ">\n";
            std::for_each(arguments_begin, arguments_end,
                          [&] (const std::string& x) { std::cout << "+     argument: " << x << "\n";});
#endif

            user_weight_function_->build(build_options.str());
            user_weight_function_->wait();

            evaluate();
        }

        template <class random_iterator>
        void generate_interval_decomposition(random_iterator first, random_iterator last)
        {
            const size_t n = last - first - 1U;
            const double delta = 1.0 / static_cast<double>(n);
            size_t i = 0U;

            while (first != last)
            {
                *first++ = static_cast<double>(i++) * delta;
            }
        }

        void evaluate()
        {
            number_of_samples_ = user_weight_function_->work_group_size(number_of_samples_);

            std::vector<double> luminances(number_of_samples_);
            generate_interval_decomposition(luminances.begin(), luminances.end());

            weight_samples_.clear();
            user_weight_function_->evaluate(luminances.begin(), luminances.end(),
                                            std::back_inserter(weight_samples_));

#if 1
            for (size_t i = 0U; i != weight_samples_.size(); ++i)
            {
                std::cout <<
                    "+ OpenCLUserExposureWeight::evaluate: Y[" << i << "] = " << luminances[i] << "\t" <<
                    "w[" << i << "] = " << weight_samples_[i] << "\n";
            }
#endif
        }

        double weight(double y) override
        {
            assert(y >= 0.0);
            assert(y <= 1.0);

            return weight_samples_[std::lround(y * static_cast<double>(number_of_samples_))];
        }

    private:
        std::string source_file_name_;
        std::string weight_function_name_;
        std::unique_ptr<UserWeightFunction> user_weight_function_;

        size_t number_of_samples_;
        std::vector<double> weight_samples_;
    }; // class OpenCLUserExposureWeight


    // If the beginning of `a_source_file_name' do not contain any
    // null-character we assume it is an OpenCL-source file.
    bool
    is_opencl_file(const std::string& a_source_file_name)
    {
        const std::streamsize buffer_size = 256U;
        std::ifstream source(a_source_file_name, std::ios::binary);
        std::unique_ptr<char[]> buffer(new char[buffer_size]);

        source.read(buffer.get(), buffer_size);
        const std::streamsize actual_size = source.gcount();

        for (std::streamsize i = 0U; i != actual_size; ++i)
        {
            if (buffer[i] == '\0')
            {
                return false;
            }
        }

        return true;
    }


    ExposureWeight*
    make_weight_function(const std::string& a_source_file_name,
                         ExposureWeight::argument_const_iterator some_arguments_begin,
                         ExposureWeight::argument_const_iterator some_arguments_end,
                         double a_y_optimum, double a_width)
    {
        if (some_arguments_begin == some_arguments_end)
        {
            // Remember that built-in exposure-weight functions never
            // take any arguments.
            std::cerr <<
                command << ": unknown built-in exposure weight function \"" << a_source_file_name << "\"" <<
                std::endl;
            exit(1);
        }
        else
        {
            ExposureWeight* weight_object = nullptr;

            try
            {
                const std::string weight_function_name = *some_arguments_begin;

                if (!GPUContext)
                {
                    throw std::runtime_error("cannot compile OpenCL user-weight function: GPU not enabled");
                }

                weight_object = new OpenCLUserExposureWeight(a_source_file_name, weight_function_name);
                weight_object->initialize(a_y_optimum, a_width,
                                          std::next(some_arguments_begin), some_arguments_end);
            }
            catch (cl::Error& an_error)
            {
                std::cerr << command << ": " << ocl::string_of_error_code(an_error.err()) << "\n";
                exit(1);
            }
            catch (ocl::runtime_error& an_error)
            {
                std::cerr << command << ": " << ocl::string_of_error_code(an_error.error().err()) << "\n";

                std::vector<std::string> messages = ocl::split_string(an_error.additional_message(), '\n', true);
                for (auto m : messages)
                {
                    std::cerr << command << ": note: " << m << "\n";
                }
                exit(1);
            }
            catch (ExposureWeight::error& an_exception)
            {
                std::cerr <<
                    command << ": user-defined OpenCL weight function \"" << "???" <<
                    "\" defined in source file \"" << a_source_file_name <<
                    "\" raised exception: " << an_exception.what() << std::endl;
                exit(1);
            }

            return weight_object;
        }
    }
#endif // OPENCL
} // namespace opencl_exposure_weight
