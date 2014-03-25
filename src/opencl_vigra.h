/*
 * Copyright (C) 2014 Christoph L. Spiel
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

#include "opencl.h"

#include "openmp_vigra.h"
#include "vigra/basicimageview.hxx"

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
                f_scratch_buffer_(nullptr), d_scratch_buffer_(nullptr),
                v_scratch_buffer_(nullptr), z_scratch_buffer_(nullptr),
                write_buffer_prereq_(1U), column_kernel_prereq_(1U), row_kernel_prereq_(1U),
                read_buffer_prereq_(1U), unmap_buffer_prereq_(1U)
            {
                f_.add_build_option("-cl-fast-relaxed-math");
                f_.add_build_option("-cl-nv-verbose");

#ifdef BUILD_EAGERLY
                f_.build();

                initialize();
#endif
            }

            void build(const std::string& a_build_option)
            {
                std::cerr << "\n+ DistanceTransformFH::build: by request\n\n";
                f_.build(a_build_option);
                std::cerr <<
                    "+ DistanceTransformFH::build: log begin ================\n" <<
                    f_.build_log() <<
                    "\n+ DistanceTransformFH::build: log end   ================\n";
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
                        command << ": warning: plain cl error: " << a_cl_error.what() << std::endl;
                }
                catch (::ocl::runtime_error& an_opencl_runtime_error)
                {
                    std::cerr <<
                        command << ": warning: falling back from OpenCL to CPU path because of\n" <<
                        command << ": warning: ocl error: " << an_opencl_runtime_error.what() << "\n" <<
                        command << ": warning:            " << an_opencl_runtime_error.additional_message() <<
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
                const vigra::Size2D size(a_source_lowerright - a_source_upperleft);
                const size_t buffer_size = size.area() * sizeof(float);

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
                const cl::NDRange local_size(16);

                setup(a_distance_norm, size);

                f_.queue().enqueueWriteBuffer(output_buffer_, CL_FALSE, 0U, buffer_size,
                                              buffer_begin,
                                              &write_buffer_prereq_, &column_kernel_prereq_[0]);
                f_.queue().enqueueNDRangeKernel(a_distance_norm >= 2 ?
                                                euclidean_column_kernel_ : manhattan_column_kernel_,
                                                cl::NullRange,
                                                column_global_size, local_size,
                                                &column_kernel_prereq_, &row_kernel_prereq_[0]);
                f_.queue().enqueueNDRangeKernel(a_distance_norm >= 2 ?
                                                euclidean_row_kernel_ : manhattan_row_kernel_,
                                                cl::NullRange,
                                                row_global_size, local_size,
                                                &row_kernel_prereq_, &read_buffer_prereq_[0]);
                f_.queue().enqueueReadBuffer(output_buffer_, CL_TRUE, 0U, buffer_size,
                                             buffer_begin,
                                             &read_buffer_prereq_, &unmap_buffer_prereq_[0]);

                vigra::omp::copyImage(distance_image.upperLeft(), distance_image.lowerRight(),
                                      distance_image.accessor(),
                                      a_destination_upperleft, a_destination_accessor);

                f_.queue().enqueueUnmapMemObject(output_buffer_, buffer_begin, &unmap_buffer_prereq_, &done_);
                show_profile_data(size);

                teardown(a_distance_norm);
            }

            void show_profile_data(const vigra::Size2D& a_size)
            {
                if (f_.queue().getInfo<CL_QUEUE_PROPERTIES>() & CL_QUEUE_PROFILING_ENABLE)
                {
                    typedef std::pair<std::string, double> pair_t;
                    std::vector<pair_t> results = {
                        std::make_pair("map buffer", ::ocl::event_latency(write_buffer_prereq_[0])),
                        std::make_pair("write buffer", ::ocl::event_latency(column_kernel_prereq_[0])),
                        std::make_pair("column kernel", ::ocl::event_latency(row_kernel_prereq_[0])),
                        std::make_pair("row kernel", ::ocl::event_latency(read_buffer_prereq_[0])),
                        std::make_pair("read buffer", ::ocl::event_latency(unmap_buffer_prereq_[0])),
                        std::make_pair("unmap buffer", ::ocl::event_latency(done_))
                    };
                    auto partial_sum =
                        [] (double a_partial_sum, const pair_t& a_pair)
                        {return a_partial_sum + a_pair.second;};
                    results.push_back(std::make_pair("TOTAL",
                                                     std::accumulate(results.begin(), results.end(),
                                                                     double(), partial_sum)));

                    const std::ios::fmtflags flags(std::cerr.flags());
                    std::cerr <<
                        command << ": timing: OpenCL latencies of `Distance Transform' for " <<
                        a_size << " = " << a_size.area() << " pixel image\n";
                    for (auto& r : results)
                    {
                        std::cerr <<
                            command << ": timing:         " <<
                            std::setw(16) << r.first <<
                            std::fixed <<
                            std::setw(10) << std::setprecision(3) << 1000.0 * r.second << " ms    " <<
                            std::setw(6) << std::setprecision(1) <<
                            100.0 * r.second / results.back().second << "%\n";
                    }
                    std::cerr <<
                        command << ": timing: device profiling timer resolution " <<
                        std::setprecision(6) <<
                        f_.device().getInfo<CL_DEVICE_PROFILING_TIMER_RESOLUTION>() / 1.0e6 << " ms\n";
                    std::cerr.flags(flags);
                }
            }

            void initialize()
            {
                if (work_group_size_)
                {
                    return;
                }

                manhattan_row_kernel_ = f_.create_kernel("manhattan_2d_rows");
                manhattan_column_kernel_ = f_.create_kernel("manhattan_2d_columns");
                euclidean_row_kernel_ = f_.create_kernel("euclidean_2d_rows");
                euclidean_column_kernel_ = f_.create_kernel("euclidean_2d_columns");

                cl::Kernel& k = euclidean_row_kernel_;
                size_t size_multiple;
                k.getWorkGroupInfo(f_.device(), CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, &size_multiple);
                preferred_work_group_size_multiple_ = static_cast<int>(size_multiple);

                size_t device_max_group_size;
                f_.device().getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &device_max_group_size);
                size_t kernel_group_size;
                k.getWorkGroupInfo(f_.device(), CL_KERNEL_WORK_GROUP_SIZE, &kernel_group_size);
                work_group_size_ = static_cast<int>(std::min(device_max_group_size, kernel_group_size));
            }

            void setup(int a_distance_norm, const vigra::Size2D& a_size)
            {
                const int max_length = std::max(a_size.width(), a_size.height());
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
                row_kernel.setArg(1U, a_size.width());
                row_kernel.setArg(2U, a_size.height());
                row_kernel.setArg(3U, *f_scratch_buffer_);
                row_kernel.setArg(4U, *d_scratch_buffer_);

                column_kernel.setArg(0U, output_buffer_);
                column_kernel.setArg(1U, a_size.width());
                column_kernel.setArg(2U, a_size.height());
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

            int work_group_size(int a_suggested_work_group_size) const
            {
                return ::ocl::round_to_next_multiple(a_suggested_work_group_size,
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

            int preferred_work_group_size_multiple_;
            int work_group_size_;

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
