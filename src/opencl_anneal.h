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
#ifndef OPENCL_ANNEAL_H_INCLUDED
#define OPENCL_ANNEAL_H_INCLUDED


/*
    Test

    cd ~/photos/multi-image-techniques/hugin-multi-row
    rm -f vis-?.tif;  env ENBLEND_OPENCL_PATH=/home/cspiel/src/enblend/src /home/cspiel/src/enblend/BUILD-GCC-O0/src/enblend --gpu --primary-seam-generator=nft --fine-mask --visualize --parameter='gpu-kernel-dt=false' --parameter='time-anneal-snake=true:time-state-probabilities=true:profile-state-probabilities=true' remapped-000?-08bit.tif

*/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "timer.h"

#include "opencl.h"


namespace ocl
{
#ifdef OPENCL

    template <typename floating_point_t>
    inline static void
    let_host_calculate_state_probabilities(int local_k,
                                           floating_point_t * RESTRICT state_probabilities,
                                           float * RESTRICT e, float * RESTRICT pi)
    {
        for (int j = 0; j < local_k; j++)
        {
            const floating_point_t pi_t_j = state_probabilities[j];
            pi[j] += pi_t_j;

            for (int i = j + 1; i < local_k; i++)
            {
                const floating_point_t pi_t = state_probabilities[i] + pi_t_j;
                floating_point_t pi_t_an = pi_t / (static_cast<floating_point_t>(1) + std::exp(e[j] - e[i]));
                if (std::isnan(pi_t_an))
                {
                    pi_t_an = e[j] > e[i] ? floating_point_t() : pi_t;
                }
                pi[j] += pi_t_an;
                pi[i] += pi_t - pi_t_an;
            }

            state_probabilities[j] = pi[j] / static_cast<floating_point_t>(local_k);
        }
    }


#ifndef PREFER_SEPARATE_OPENCL_SOURCE
#include "calculate_state_probabilities.icl"
#endif

    class CalculateStateProbabilities : public ::ocl::BuildableFunction
    {
        enum write_prereq_indexes
        {
            E_MAPPED, PI_MAPPED,
            N_MAPPED_
        };

        enum kernel_prereq_indexes
        {
            E_BUFFER_WRITTEN, PI_BUFFER_WRITTEN, STATE_PROBABILITIES_BUFFER_WRITTEN,
            N_WRITTEN_
        };

        enum unmap_prereq_indexes
        {
            PI_BUFFER_UPDATED, STATE_PROBABILITIES_BUFFER_UPDATED,
            N_UPDATED_
        };

    public:
        CalculateStateProbabilities() = delete;

        CalculateStateProbabilities(const cl::Context& a_context) :
            immediately_fallback_(false),
#ifdef PREFER_SEPARATE_OPENCL_SOURCE
            f_(a_context, std::string("calculate_state_probabilities.cl")),
#else
            f_(a_context, calculate_state_probabilities_source_code),
#endif
            preferred_work_group_size_multiple_(0U), work_group_size_(0U),
            e_begin_(nullptr), pi_begin_(nullptr),
            map_complete_(N_MAPPED_), kernel_prereq_(N_WRITTEN_),
            read_buffer_prereq_(1U), unmap_buffer_prereq_(N_UPDATED_)
        {
            std::string device_extensions;
            f_.device().getInfo(CL_DEVICE_EXTENSIONS, &device_extensions);
            extensions_ = split_string(device_extensions, ' ');
            std::sort(extensions_.begin(), extensions_.end());

            f_.add_build_option("-cl-fast-relaxed-math");
            f_.add_build_option("-cl-unsafe-math-optimizations");

            switch (f_.vendor_id())
            {
            case ::ocl::vendor::amd:
                // f_.add_build_option("...");
                break;
            case ::ocl::vendor::nvidia:
                f_.add_build_option("-cl-nv-verbose");
                break;
            case ::ocl::vendor::unknown:
                break;
            }

#ifdef BUILD_EAGERLY
            f_.build();

            initialize();
#endif
        }

        bool has_extension(const std::string& an_extension)
        {
            return std::binary_search(extensions_.begin(), extensions_.end(), an_extension);
        }

        void build(const std::string& a_build_option)
        {
            std::cerr << "\n+ CalculateStateProbabilities::build: by request\n\n";
            f_.build(a_build_option);
            std::cerr <<
                "+ CalculateStateProbabilities::build: log begin ================\n" <<
                f_.build_log() <<
                "\n+ CalculateStateProbabilities::build: log end   ================\n";
        }

        void wait()
        {
            f_.wait();
#ifndef BUILD_EAGERLY
            std::cerr << "\n+ CalculateStateProbabilities::wait: carry on...\n\n";
            initialize();
#endif
        }

        void run(int local_k, std::vector<double>* state_probabilities, int k_max, float* e, float* pi)
        {
            if (EXPECT_RESULT(!immediately_fallback_, true))
            {
                try
                {
                    run0(local_k, state_probabilities, k_max, e, pi);
                    return;
                }
                catch (cl::Error& a_cl_error)
                {
                    std::cerr <<
                        command << ": warning: falling back from OpenCL to CPU path because of\n" <<
                        command << ": warning: plain cl error in function: " << a_cl_error.what() << "\n" <<
                        command << ": note: " << ::ocl::string_of_error_code(a_cl_error.err()) <<
                        std::endl;
                }
                catch (::ocl::runtime_error& an_opencl_runtime_error)
                {
                    std::cerr <<
                        command << ": warning: falling back from OpenCL to CPU path because of\n" <<
                        command << ": warning: ocl error in function: " <<
                        an_opencl_runtime_error.what() << "\n" <<
                        command << ": note: " <<
                        ::ocl::string_of_error_code(an_opencl_runtime_error.error().err()) << "\n" <<
                        command << ": note: " << an_opencl_runtime_error.additional_message() <<
                        std::endl;
                }
            }

            let_host_calculate_state_probabilities<double>(local_k, &(*state_probabilities)[0], e, pi);
        }

        // In setup() the parameter `size' is the maximum number of
        // elements of any of the `state_probabilities' vectors.
        void setup(size_t size, size_t k_max, float*& e, float*& pi)
        {
            const size_t local_k = k_max;

            state_probabilities_buffer_ = cl::Buffer(f_.context(), CL_MEM_READ_WRITE, size * sizeof(double));
            e_buffer_ = cl::Buffer(f_.context(), CL_MEM_ALLOC_HOST_PTR, k_max * sizeof(float));
            pi_buffer_ = cl::Buffer(f_.context(), CL_MEM_ALLOC_HOST_PTR, k_max * sizeof(float));

            // kernel argument #0, i.e. k_max, will be set in method run().
            state_probabilities_kernel_.setArg(1U, state_probabilities_buffer_);
            state_probabilities_kernel_.setArg(2U, cl::Local(size * sizeof(float)));
            state_probabilities_kernel_.setArg(3U, e_buffer_);
            state_probabilities_kernel_.setArg(4U, pi_buffer_);
            state_probabilities_kernel_.setArg(5U, cl::Local(k_max * sizeof(float)));

            e_begin_ = static_cast<float*>(f_.queue().enqueueMapBuffer(e_buffer_, CL_FALSE,
                                                                       CL_MEM_READ_ONLY,
                                                                       0U, local_k * sizeof(float),
                                                                       nullptr, &map_complete_[E_MAPPED]));
            pi_begin_ = static_cast<float*>(f_.queue().enqueueMapBuffer(pi_buffer_, CL_FALSE,
                                                                        CL_MEM_READ_WRITE,
                                                                        0U, local_k * sizeof(float),
                                                                        nullptr, &map_complete_[PI_MAPPED]));
            cl::Event::waitForEvents(map_complete_);

            e = e_begin_;
            pi = pi_begin_;
        }

        void teardown()
        {
            f_.queue().enqueueUnmapMemObject(e_buffer_, e_begin_, &unmap_buffer_prereq_);
            f_.queue().enqueueUnmapMemObject(pi_buffer_, pi_begin_, &unmap_buffer_prereq_);

            e_begin_ = nullptr;
            pi_begin_ = nullptr;
        }

    private:
        void run0(int local_k, std::vector<double>* state_probabilities, int k_max, float* e, float* pi)
        {
            if (!has_extension("cl_khr_fp64"))
            {
                immediately_fallback_ = true;
                throw ocl::runtime_error(cl::Error(CL_BUILD_PROGRAM_FAILURE),
                                         "device does not support double-precision floating-point numbers");
            }

            const int size = static_cast<int>(state_probabilities->size());
            const size_t buffer_size = size * sizeof(double);
            double* const state_probabilities_begin = &(*state_probabilities)[0];

            state_probabilities_kernel_.setArg(0U, static_cast<cl_int>(local_k));

            f_.queue().enqueueWriteBuffer(state_probabilities_buffer_, CL_FALSE,
                                          0U, buffer_size,
                                          state_probabilities_begin,
                                          nullptr, // no prerequisite
                                          &kernel_prereq_[STATE_PROBABILITIES_BUFFER_WRITTEN]);
            f_.queue().enqueueWriteBuffer(e_buffer_, CL_FALSE,
                                          0U, local_k * sizeof(float),
                                          e_begin_,
                                          &map_complete_,
                                          &kernel_prereq_[E_BUFFER_WRITTEN]);
            f_.queue().enqueueWriteBuffer(pi_buffer_, CL_FALSE,
                                          0U, local_k * sizeof(float),
                                          pi_begin_,
                                          &map_complete_,
                                          &kernel_prereq_[PI_BUFFER_WRITTEN]);

            f_.queue().enqueueNDRangeKernel(state_probabilities_kernel_,
                                            cl::NullRange,
                                            cl::NDRange(work_group_size(k_max)), // global size
                                            cl::NDRange(work_group_size(k_max)), //cl::NullRange, // local size
                                            &kernel_prereq_,
                                            &read_buffer_prereq_[0]);
            DEBUG_CHECK_OPENCL_EVENT(read_buffer_prereq_[0]);

            f_.queue().enqueueReadBuffer(state_probabilities_buffer_, CL_FALSE,
                                         0U, buffer_size,
                                         state_probabilities_begin,
                                         &read_buffer_prereq_,
                                         &unmap_buffer_prereq_[STATE_PROBABILITIES_BUFFER_UPDATED]);
            f_.queue().enqueueReadBuffer(pi_buffer_, CL_FALSE,
                                         0U, local_k * sizeof(float),
                                         pi_begin_,
                                         &read_buffer_prereq_,
                                         &unmap_buffer_prereq_[PI_BUFFER_UPDATED]);

            cl::Event::waitForEvents(unmap_buffer_prereq_);

            if (parameter::as_boolean("profile-state-probabilities", false))
            {
                show_profile_data(buffer_size, local_k);
            }
        }

        void show_profile_data(size_t buffer_size, int local_k)
        {
            if (f_.queue().getInfo<CL_QUEUE_PROPERTIES>() & CL_QUEUE_PROFILING_ENABLE)
            {
                ShowProfileData show_profile;

                show_profile.set_device(&f_.device());

                show_profile.add_event_latencies("write buffer",
                                                 kernel_prereq_.begin(), kernel_prereq_.end());
                show_profile.add_event_latency("kernel", read_buffer_prereq_[0]);
                show_profile.add_event_latencies("read buffer",
                                                 unmap_buffer_prereq_.begin(), unmap_buffer_prereq_.end());

                const double delta_t = 1e-3 * f_.device().getInfo<CL_DEVICE_PROFILING_TIMER_RESOLUTION>();

                {
                    const double n = static_cast<double>(local_k * (local_k + 1));
                    const double t = show_profile.retrieve_result("kernel", ShowProfileData::MICRO_SECONDS);

                    std::cerr <<
                        "\n" <<
                        command << ": timing: OpenCL latencies of `Calculate State Probabilities' for k = " <<
                        local_k << ".\n" <<
                        command << ": timing: Kernel performance " << n / t <<
                        "±" << n / t * delta_t / t << " probabilities/µs.\n";
                }

                {
                    const double theta =
                        static_cast<double>(2 * buffer_size + 3 * local_k * sizeof(float)) / 1024.0;
                    const double t =
                        show_profile.retrieve_result("write buffer", ShowProfileData::MICRO_SECONDS) +
                        show_profile.retrieve_result("read buffer", ShowProfileData::MICRO_SECONDS);

                    std::cerr <<
                        command << ": timing: Effective bandwidth " << theta / t <<
                        "±" << theta / t * delta_t / t << " kB/µs (≙ GB/s).\n";
                }

                show_profile.show_results(std::cerr, command + ": ", ShowProfileData::MICRO_SECONDS);
            }
        }

        void initialize()
        {
            if (EXPECT_RESULT(work_group_size_ > 0U, true))
            {
                return;
            }

            state_probabilities_kernel_ = f_.create_kernel("calculate_state_probabilities");

            cl::Kernel& k = state_probabilities_kernel_;
            k.getWorkGroupInfo(f_.device(),
                               CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                               &preferred_work_group_size_multiple_);

            size_t device_max_group_size;
            f_.device().getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &device_max_group_size);
            size_t kernel_group_size;
            k.getWorkGroupInfo(f_.device(), CL_KERNEL_WORK_GROUP_SIZE, &kernel_group_size);
            work_group_size_ = std::min(device_max_group_size, kernel_group_size);
        }

        size_t work_group_size(size_t a_suggested_work_group_size) const
        {
            return ::ocl::round_up_to_next_multiple(a_suggested_work_group_size,
                                                    preferred_work_group_size_multiple_);
        }

        bool immediately_fallback_;

#ifdef PREFER_SEPARATE_OPENCL_SOURCE
        ::ocl::LazyFunctionCXXOfFile f_;
#else
        ::ocl::LazyFunctionCXXOfString f_;
#endif

        cl::Kernel state_probabilities_kernel_;

        std::vector<std::string> extensions_;
        size_t preferred_work_group_size_multiple_;
        size_t work_group_size_;

        cl::Buffer state_probabilities_buffer_;
        cl::Buffer e_buffer_;
        cl::Buffer pi_buffer_;

        float* e_begin_;
        float* pi_begin_;

        std::vector<cl::Event> map_complete_;
        std::vector<cl::Event> kernel_prereq_;
        std::vector<cl::Event> read_buffer_prereq_;
        std::vector<cl::Event> unmap_buffer_prereq_;
    }; // class CalculateStateProbabilities

#endif // OPENCL
} // namespace ocl


#endif // OPENCL_ANNEAL_H_INCLUDED


// Local Variables:
// mode: c++
// End:
