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
#ifndef OPENCL_VIGRA_H_INCLUDED
#define OPENCL_VIGRA_H_INCLUDED


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vigra/basicimageview.hxx>

#include "muopt.h"
#include "opencl.h"
#include "openmp_vigra.h"


namespace vigra
{
    namespace ocl
    {
#ifdef OPENCL

#ifndef PREFER_SEPARATE_OPENCL_SOURCE
#include "distance_transform_fh.icl"
#endif

        class DistanceTransformFH : public ::ocl::BuildableFunction
        {
        public:
            DistanceTransformFH() = delete;

            DistanceTransformFH(const cl::Context& a_context) :
#ifdef PREFER_SEPARATE_OPENCL_SOURCE
                f_(a_context, std::string("distance_transform_fh.cl")),
#else
                f_(a_context, distance_transform_fh_source_code),
#endif
                preferred_work_group_size_multiple_(0U), work_group_size_(0U),
                f_scratch_buffer_(nullptr), d_scratch_buffer_(nullptr),
                v_scratch_buffer_(nullptr), z_scratch_buffer_(nullptr),
                write_buffer_prereq_(1U), column_kernel_prereq_(1U), row_kernel_prereq_(1U),
                read_buffer_prereq_(1U), unmap_buffer_prereq_(1U)
            {
                f_.add_build_option("-cl-fast-relaxed-math");
                f_.add_build_option("-cl-strict-aliasing");

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
                build("");
                initialize();
#endif
            }

            void build(const std::string& a_build_option)
            {
                try
                {
                    std::cerr << "\n+ DistanceTransformFH::build: by request\n\n";
                    f_.build(a_build_option);
                    std::cerr <<
                        "+ DistanceTransformFH::build: log begin ================\n" <<
                        f_.build_log() <<
                        "\n+ DistanceTransformFH::build: log end   ================\n";
                }
                catch (::ocl::runtime_error& an_error)
                {
                    std::cerr <<
                        command << ": " << ::ocl::string_of_error_code(an_error.error().err()) << "\n";

                    std::vector<std::string> messages =
                        ::ocl::split_string(an_error.additional_message(), '\n', true);
                    for (auto m : messages)
                    {
                        std::cerr << command << ": note: " << m << "\n";
                    }

                    exit(1);
                }
            }

            void wait()
            {
                f_.wait();
#ifndef BUILD_EAGERLY
                std::cerr << "\n+ DistanceTransformFH::wait: carry on...\n\n";
                initialize();
#endif
            }

            template <class source_iterator, class source_accessor,
                      class destination_iterator, class destination_accessor,
                      class value_type>
            void run(source_iterator a_source_upperleft, source_iterator a_source_lowerright,
                     source_accessor a_source_accessor,
                     destination_iterator a_destination_upperleft, destination_accessor a_destination_accessor,
                     value_type a_background_value,
                     int a_distance_norm)
            {
                try
                {
                    run0(a_source_upperleft, a_source_lowerright, a_source_accessor,
                         a_destination_upperleft, a_destination_accessor,
                         a_background_value,
                         a_distance_norm);
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

                vigra::omp::distanceTransform(a_source_upperleft, a_source_lowerright, a_source_accessor,
                                              a_destination_upperleft, a_destination_accessor,
                                              a_background_value,
                                              a_distance_norm);
            }

            template <class source_iterator, class source_accessor,
                      class destination_iterator, class destination_accessor,
                      class value_type>
            void run(vigra::triple<source_iterator, source_iterator, source_accessor> a_source,
                     vigra::pair<destination_iterator, destination_accessor> a_destination,
                     value_type a_background,
                     int a_distance_norm)
            {
                run(a_source.first, a_source.second, a_source.third,
                    a_destination.first, a_destination.second,
                    a_background,
                    a_distance_norm);
            }

        private:
            template <class source_iterator, class source_accessor,
                      class destination_iterator, class destination_accessor,
                      class value_type>
            void run0(source_iterator a_source_upperleft, source_iterator a_source_lowerright,
                      source_accessor a_source_accessor,
                      destination_iterator a_destination_upperleft, destination_accessor a_destination_accessor,
                      value_type a_background_value,
                      int a_distance_norm)
            {
                wait();         // Ensure that all kernels were built.

                const vigra::Size2D size(a_source_lowerright - a_source_upperleft);
                const size_t buffer_size = static_cast<size_t>(size.area()) * sizeof(float);

                output_buffer_ = cl::Buffer(f_.context(), CL_MEM_ALLOC_HOST_PTR, buffer_size);
                float* const buffer_begin =
                    static_cast<float*>(f_.queue().enqueueMapBuffer(output_buffer_, CL_TRUE,
                                                                    CL_MAP_READ | CL_MAP_WRITE_INVALIDATE_REGION,
                                                                    0U, buffer_size,
                                                                    nullptr, &write_buffer_prereq_[0]));
                vigra::BasicImageView<float> distance_image(buffer_begin, size);

                vigra::omp::transformImage(a_source_upperleft, a_source_lowerright, a_source_accessor,
                                           distance_image.upperLeft(), distance_image.accessor(),
                                           vigra::functor::ifThenElse(vigra::functor::Arg1() >
                                                                      vigra::functor::Param(a_background_value),
                                                                      vigra::functor::Param(0.0f),
                                                                      vigra::functor::Param(FLT_MAX)));

                const cl::NDRange column_global_size(work_group_size(size.width()));
                const cl::NDRange row_global_size(work_group_size(size.height()));
                const cl::NDRange local_size(16U);

                setup(a_distance_norm, size);

                f_.queue().enqueueWriteBuffer(output_buffer_, CL_FALSE, 0U, buffer_size,
                                              buffer_begin,
                                              &write_buffer_prereq_, &column_kernel_prereq_[0]);
                f_.queue().enqueueNDRangeKernel(a_distance_norm >= 2 ?
                                                euclidean_column_kernel_ : manhattan_column_kernel_,
                                                cl::NullRange,
                                                column_global_size, local_size,
                                                &column_kernel_prereq_, &row_kernel_prereq_[0]);
                DEBUG_CHECK_OPENCL_EVENT(row_kernel_prereq_[0]);
                f_.queue().enqueueNDRangeKernel(a_distance_norm >= 2 ?
                                                euclidean_row_kernel_ : manhattan_row_kernel_,
                                                cl::NullRange,
                                                row_global_size, local_size,
                                                &row_kernel_prereq_, &read_buffer_prereq_[0]);
                DEBUG_CHECK_OPENCL_EVENT(read_buffer_prereq_[0]);
                f_.queue().enqueueReadBuffer(output_buffer_, CL_TRUE, 0U, buffer_size,
                                             buffer_begin,
                                             &read_buffer_prereq_, &unmap_buffer_prereq_[0]);

                vigra::omp::copyImage(distance_image.upperLeft(), distance_image.lowerRight(),
                                      distance_image.accessor(),
                                      a_destination_upperleft, a_destination_accessor);

                f_.queue().enqueueUnmapMemObject(output_buffer_, buffer_begin, &unmap_buffer_prereq_, &done_);

                if (parameter::as_boolean("time-distance-transform", false))
                {
                    show_profile_data(size);
                }

                teardown(a_distance_norm);
            }

            void show_profile_data(const vigra::Size2D& a_size)
            {
                if (f_.queue().getInfo<CL_QUEUE_PROPERTIES>() & CL_QUEUE_PROFILING_ENABLE)
                {
                    ::ocl::ShowProfileData show_profile;

                    show_profile.set_device(&f_.device());

                    show_profile.add_event_latency("map buffer", write_buffer_prereq_[0]);
                    show_profile.add_event_latency("write buffer", column_kernel_prereq_[0]);
                    show_profile.add_event_latency("column kernel", row_kernel_prereq_[0]);
                    show_profile.add_event_latency("row kernel", read_buffer_prereq_[0]);
                    show_profile.add_event_latency("read buffer", unmap_buffer_prereq_[0]);
                    show_profile.add_event_latency("unmap buffer", done_);

                    const double n = static_cast<double>(a_size.area());
                    const double delta_t = 1e-3 * f_.device().getInfo<CL_DEVICE_PROFILING_TIMER_RESOLUTION>();

                    {
                        const double t_column =
                            show_profile.retrieve_result("column kernel", ::ocl::ShowProfileData::MICRO_SECONDS);
                        const double t_row =
                            show_profile.retrieve_result("row kernel", ::ocl::ShowProfileData::MICRO_SECONDS);

                        std::cerr <<
                            "\n" <<
                            command << ": timing: OpenCL latencies of `Distance Transform' for\n" <<
                            command << ": timing: " << a_size << " = " << n / 1028196.0 << " Mpixel image.\n" <<
                            command << ": timing: Column-kernel performance " << n / t_column <<
                            "±" << n / t_column * delta_t / t_column << " pixels/µs.\n" <<
                            command << ": timing: Row-kernel performance " << n / t_row <<
                            "±" << n / t_row * delta_t / t_row << " pixels/µs.\n";
                    }

                    {
                        const double theta = 2.0 * n * sizeof(float) / 1024.0;
                        const double t =
                            show_profile.retrieve_result("write buffer", ::ocl::ShowProfileData::MICRO_SECONDS) +
                            show_profile.retrieve_result("read buffer", ::ocl::ShowProfileData::MICRO_SECONDS);

                        std::cerr <<
                            command << ": timing: Effective bandwidth " << theta / t <<
                            "±" << theta / t * delta_t / t << " kB/µs (≙ GB/s).\n";
                    }

                    show_profile.show_results(std::cerr, command + ": ", ::ocl::ShowProfileData::MICRO_SECONDS);
                }
            }

            void initialize()
            {
                if (EXPECT_RESULT(work_group_size_ > 0U, true))
                {
                    return;
                }

                manhattan_row_kernel_ = f_.create_kernel("manhattan_2d_rows");
                manhattan_column_kernel_ = f_.create_kernel("manhattan_2d_columns");
                euclidean_row_kernel_ = f_.create_kernel("euclidean_2d_rows");
                euclidean_column_kernel_ = f_.create_kernel("euclidean_2d_columns");

                cl::Kernel& k = euclidean_row_kernel_;
                k.getWorkGroupInfo(f_.device(),
                                   CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                   &preferred_work_group_size_multiple_);

                size_t device_max_group_size;
                f_.device().getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &device_max_group_size);
                size_t kernel_group_size;
                k.getWorkGroupInfo(f_.device(), CL_KERNEL_WORK_GROUP_SIZE, &kernel_group_size);
                work_group_size_ = std::min(device_max_group_size, kernel_group_size);
            }

            void setup(int a_distance_norm, const vigra::Size2D& a_size)
            {
                const size_t max_length = static_cast<size_t>(std::max(a_size.width(), a_size.height()));
                const size_t float_size = max_length * max_length * sizeof(cl_float);
                const size_t int_size = max_length * max_length * sizeof(cl_int);

                f_scratch_buffer_ = new cl::Buffer(f_.context(), CL_MEM_READ_WRITE, float_size);
                d_scratch_buffer_ = new cl::Buffer(f_.context(), CL_MEM_READ_WRITE, float_size);

                cl::Kernel row_kernel;
                cl::Kernel column_kernel;
                if (a_distance_norm >= 2)
                {
                    row_kernel = euclidean_row_kernel_;
                    column_kernel = euclidean_column_kernel_;
                }
                else
                {
                    row_kernel = manhattan_row_kernel_;
                    column_kernel = manhattan_column_kernel_;
                }

                row_kernel.setArg(0U, output_buffer_);
                row_kernel.setArg(1U, static_cast<cl_int>(a_size.width()));
                row_kernel.setArg(2U, static_cast<cl_int>(a_size.height()));
                row_kernel.setArg(3U, *f_scratch_buffer_);
                row_kernel.setArg(4U, *d_scratch_buffer_);

                column_kernel.setArg(0U, output_buffer_);
                column_kernel.setArg(1U, static_cast<cl_int>(a_size.width()));
                column_kernel.setArg(2U, static_cast<cl_int>(a_size.height()));
                column_kernel.setArg(3U, *f_scratch_buffer_);
                column_kernel.setArg(4U, *d_scratch_buffer_);

                if (a_distance_norm >= 2)
                {
                    v_scratch_buffer_ = new cl::Buffer(f_.context(), CL_MEM_READ_WRITE, int_size);
                    z_scratch_buffer_ = new cl::Buffer(f_.context(), CL_MEM_READ_WRITE, float_size);

                    row_kernel.setArg(5U, *v_scratch_buffer_);
                    row_kernel.setArg(6U, *z_scratch_buffer_);

                    column_kernel.setArg(5U, *v_scratch_buffer_);
                    column_kernel.setArg(6U, *z_scratch_buffer_);
                }
            }

            void teardown(int a_distance_norm __attribute__((unused)))
            {
                delete z_scratch_buffer_;
                delete d_scratch_buffer_;
                delete v_scratch_buffer_;
                delete f_scratch_buffer_;
            }

            size_t work_group_size(size_t a_suggested_work_group_size) const
            {
                return ::ocl::round_up_to_next_multiple(a_suggested_work_group_size,
                                                        preferred_work_group_size_multiple_);
            }

#ifdef PREFER_SEPARATE_OPENCL_SOURCE
            ::ocl::LazyFunctionCXXOfFile f_;
#else
            ::ocl::LazyFunctionCXXOfString f_;
#endif

            cl::Kernel manhattan_row_kernel_;
            cl::Kernel manhattan_column_kernel_;
            cl::Kernel euclidean_row_kernel_;
            cl::Kernel euclidean_column_kernel_;

            size_t preferred_work_group_size_multiple_;
            size_t work_group_size_;

            cl::Buffer output_buffer_;
            cl::Buffer* f_scratch_buffer_;
            cl::Buffer* d_scratch_buffer_;
            cl::Buffer* v_scratch_buffer_;
            cl::Buffer* z_scratch_buffer_;

            std::vector<cl::Event> write_buffer_prereq_;
            std::vector<cl::Event> column_kernel_prereq_;
            std::vector<cl::Event> row_kernel_prereq_;
            std::vector<cl::Event> read_buffer_prereq_;
            std::vector<cl::Event> unmap_buffer_prereq_;
            cl::Event done_;
        }; // class DistanceTransformFH

#endif // OPENCL
    } // namespace vigra
} // namespace ocl


#endif // OPENCL_VIGRA_H_INCLUDED

// Local Variables:
// mode: c++
// End:
