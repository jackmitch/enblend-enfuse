/*
 * Copyright (C) 2013, 2014 Christoph L. Spiel
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

#include <condition_variable>
#include <cstdint>              // std::uint8_t
#include <deque>
#include <memory>               // std::unique_ptr
#include <mutex>
#include <stdexcept>            // std::runtime_error
#include <string>
#include <thread>
#include <vector>

#define __CL_ENABLE_EXCEPTIONS

#if defined(HAVE_CL_CL_HPP)
#include <CL/cl.hpp>
#elif defined(HAVE_OPENCL_CL_HPP)
#include <OpenCL/cl.hpp>
#endif

#ifdef _MSC_VER
#define NOINLINE __declspec(noinline)
#define UNUSEDVAR
#define strncasecmp(s1, s2, n) _strnicmp(s1, s2, n)
#else
#define NOINLINE __attribute__((noinline))
#define UNUSEDVAR __attribute__((unused))
#endif

namespace ocl
{
#if defined(_OPENCL) || defined(CL_HPP_)

#define OPENCL

    typedef std::vector<cl::Platform> platform_list_t;
    typedef std::vector<cl::Device> device_list_t;


    class runtime_error : public std::runtime_error
    {
    public:
        runtime_error() = delete;
        runtime_error(const std::string& a_message);
        runtime_error(const cl::Error& an_opencl_error,
                      const std::string& an_additional_message);

        virtual ~runtime_error() throw() {}

        const cl::Error& error() const;
        const std::string& additional_message() const;

    private:
        cl::Error opencl_error_;
        std::string additional_message_;
    }; // class runtime_error


    std::string string_of_error_code(cl_int error_code);

    void print_opencl_information(bool all_devices = false);
    void print_gpu_preference(size_t a_preferred_platform_id, size_t a_preferred_device_id);

    cl::Platform find_platform(/* input/output */ size_t& a_preferred_platform_id);
    void prefer_device(const cl::Platform& a_platform,
                       size_t a_preferred_platform_id,
                       size_t a_preferred_device_id,
                       /* output */ device_list_t& some_devices);

    // Create a new context given a_platform and some_devices.  This
    // pointer can be deleted as usual.
    cl::Context* create_context(const cl::Platform& a_platform, const device_list_t& some_devices);


    namespace vendor
    {
        typedef enum
        {
            amd,
            nvidia,
            unknown
        } id_t;
    } // namespace vendor


    // Recover the platform vendor's ID from a_context.
    vendor::id_t derive_vendor_id_from_context(const cl::Context* a_context);


    template <typename t>
    inline static t
    round_to_next_multiple(t a_number, t a_multiple)
    {
        if (a_multiple == t())
        {
            return a_number;
        }
        else
        {
            const t remainder = a_number % a_multiple;

            if (remainder == t())
            {
                return a_number;
            }
            else
            {
                return a_number + a_multiple - remainder;
            }
        }
    }


    double event_latency(cl::Event& an_event); // latency in seconds


    // Macro CHECK_OPENCL_EVENT is useful to debug out-of-bound
    // read/write operations of OpenCL kernels.  To debug add
    // invocations of the macro right after each OpenCL
    // function-call that has a "task-completed event" argument,
    // m_event being the "task-completed event".
    //
    // With CHECK_OPENCL_EVENT out-of-bound reads/writes often
    // show up as exception cl::Error/"resource exhausted" and
    // usually the offending kernel is the one to which m_event is
    // attached.
    void check_opencl_event(cl::Event& an_event, const char* a_filename, int a_linenumber);
#define CHECK_OPENCL_EVENT(m_event) ::ocl::check_opencl_event(m_event, __FILE__, __LINE__)

#ifdef DEBUG
#define DEBUG_CHECK_OPENCL_EVENT(m_event) CHECK_OPENCL_EVENT(m_event)
#else
#define DEBUG_CHECK_OPENCL_EVENT(m_event)
#endif


    ////////////////////////////////////////////////////////////////////////////


    class CodePolicy
    {
    public:
        virtual ~CodePolicy() {}
    }; // class CodePolicy


    class SourcePolicy : public CodePolicy
    {
    public:
        virtual const std::string& text() = 0;
        std::pair<const char*, size_t> source();
    }; // class SourcePolicy


    class BinaryPolicy : public CodePolicy
    {
    public:
        typedef std::vector<std::uint8_t> code_t;

        virtual code_t code() = 0;
        std::pair<const void*, size_t> binary();
    }; // class BinaryPolicy


    class SourceStringPolicy : public SourcePolicy
    {
    public:
        SourceStringPolicy() = delete;
        explicit SourceStringPolicy(const std::string& a_source_text);

        const std::string& text() override;

    private:
        std::string text_;
    }; // class SourceStringPolicy


    class SourceFilePolicy : public SourcePolicy
    {
    public:
        SourceFilePolicy() = delete;
        explicit SourceFilePolicy(const std::string& a_source_filename);

        std::string filename() const {return filename_;}
        const std::string& text() override;

    private:
        NOINLINE void consult();

        std::string filename_;
        std::string text_;
    }; // class SourceFilePolicy


    class BinaryCodePolicy : public BinaryPolicy
    {
    public:
        BinaryCodePolicy() = delete;
        explicit BinaryCodePolicy(const code_t& a_binary_code);

        code_t code() override {return code_;}

    private:
        code_t code_;
    }; // class BinaryCodePolicy


    class BinaryFilePolicy : public BinaryPolicy
    {
    public:
        BinaryFilePolicy() = delete;
        explicit BinaryFilePolicy(const std::string& a_binary_filename);

        std::string filename() const {return filename_;}
        code_t code() override;

    private:
        NOINLINE void consult();

        std::string filename_;
        code_t code_;
    }; // class BinaryFilePolicy


    ////////////////////////////////////////////////////////////////////////////


    class BuildableFunction
    {
    public:
        virtual ~BuildableFunction() {}

        virtual void build(const std::string& a_build_option) = 0;
        virtual void wait() = 0;
    }; // class BuildableFunction


    // Class "Function" is the main helper for constructing
    // OpenCL-based Vigra extensions.
    //     * It supplies access to cl::Context, one or more
    //     * cl::Device objects, as well as one or more associated
    //       cl::CommandQueue objects.
    //     * It frees the developer from caring about the origin
    //       of the source or binary code.
    //     * It assists in getting OpenCL source code compiled.

    template <class actual_code_policy,
              int default_queue_flags = CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE>
    class Function : public BuildableFunction, public actual_code_policy
    {
    public:
        typedef actual_code_policy code_policy;

        Function() = delete;
        Function(const cl::Context& a_context, const std::string& a_string);
        virtual ~Function() {finalize();}

        void clear_build_options();

        Function& add_build_option(const std::string& an_option);
        Function& add_build_option(const char* a_format_string, ...);

        virtual void build(const std::string& an_extra_build_option = std::string());
        virtual void build(const char* a_format_string, ...);

        std::vector<std::string> build_logs() const;
        std::string build_log() const;

        std::vector<BinaryPolicy::code_t> binaries() const;
        BinaryPolicy::code_t binary() const;

        const cl::Context& context() const;

        vendor::id_t vendor_id() const;

        const std::vector<cl::Device>& devices() const;
        const cl::Device& device() const;

        virtual void wait() {}  // Wait for device until source has been compiled.

        virtual const cl::Program& program();

        cl::Kernel create_kernel(const std::string& an_entry_point);

        std::string build_options(const std::string& an_extra_build_option) const;

        const std::vector<cl::CommandQueue>& queues() const {return queues_;}
        cl::CommandQueue queue() const {return queues_.front();}

    protected:
        virtual void update_program_from_source(const cl::Program::Sources& a_source);

    private:
        void initialize();
        void finalize();

        cl::Program program_;
        cl::Context context_;
        vendor::id_t vendor_id_;
        std::vector<cl::Device> devices_;
        std::vector<cl::CommandQueue> queues_;
        std::vector<std::string> build_options_;
    }; // class Function


    class FunctionOfString : public Function<SourceStringPolicy>
    {
    public:
        FunctionOfString() = delete;
        FunctionOfString(const cl::Context& a_context, const std::string& a_string) :
            Function<SourceStringPolicy>(a_context, a_string) {}
    }; // class FunctionOfString


    class FunctionOfFile : public Function<SourceFilePolicy>
    {
    public:
        FunctionOfFile() = delete;
        FunctionOfFile(const cl::Context& a_context, const std::string& a_source_filename) :
            Function<SourceFilePolicy>(a_context, a_source_filename) {}
    }; // class FunctionOfFile


    template <class actual_code_policy>
    class LazyFunction : public Function<actual_code_policy>
    {
        typedef Function<actual_code_policy> super;

    public:
        typedef actual_code_policy code_policy;

        LazyFunction(const cl::Context& a_context, const std::string& a_string);

        void build(const std::string& an_extra_build_option = std::string()) override;

        const cl::Program& program() override
        {
            wait();
            return super::program();
        }

        virtual void notify(cl_program a_program) = 0;

    protected:
        // IMPLEMENTATION NOTE: This method could be 'const' in this
        // class, but derived classes may want to change the instance
        // (*not* the bool itself, of course), e.g. to lock the access
        // in a multi-threaded environment.
        virtual bool build_completed() {return build_completed_;}

        void set_build_completed(bool has_completed) {build_completed_ = has_completed;}

    private:
        static void notify_trampoline(cl_program a_program, void* an_instance);

        void update_hashes(const std::string& an_extra_build_option);
        bool needs_building(const std::string& an_extra_build_option);

        bool build_completed_;
        size_t text_hash_;
        size_t build_option_hash_;
    }; // class LazyFunction


    template <class actual_code_policy>
    class LazyFunctionCXX : public LazyFunction<actual_code_policy>
    {
        typedef LazyFunction<actual_code_policy> super;

    public:
        typedef actual_code_policy code_policy;

        LazyFunctionCXX() = delete;
        LazyFunctionCXX(const cl::Context& a_context, const std::string& a_string);
        LazyFunctionCXX(const LazyFunctionCXX&) = delete;
        LazyFunctionCXX& operator=(const LazyFunctionCXX&) = delete;

        void wait() override;

        bool build_completed() override;

        void notify(cl_program a_program);

    private:
        std::mutex build_completed_mutex_;
        std::condition_variable build_completed_condition_;
    }; //  class LazyFunctionCXX


    class LazyFunctionCXXOfString : public LazyFunctionCXX<SourceStringPolicy>
    {
    public:
        LazyFunctionCXXOfString() = delete;
        LazyFunctionCXXOfString(const cl::Context& a_context, const std::string& a_string) :
            LazyFunctionCXX<SourceStringPolicy>(a_context, a_string) {}
    }; // class LazyFunctionCXXOfString


    class LazyFunctionCXXOfFile : public LazyFunctionCXX<SourceFilePolicy>
    {
    public:
        LazyFunctionCXXOfFile() = delete;
        LazyFunctionCXXOfFile(const cl::Context& a_context, const std::string& a_source_filename) :
            LazyFunctionCXX<SourceFilePolicy>(a_context, a_source_filename) {}
    }; // class LazyFunctionCXXOfFile


    template <class ocl_function>
    std::unique_ptr<ocl_function>
    create_function(cl::Context* a_context)
    {
        try
        {
            return std::unique_ptr<ocl_function>(new ocl_function(*a_context));
        }
        catch (ocl::runtime_error& a_runtime_error)
        {
#ifdef DEBUG
            std::cerr <<
                "+ ocl::create_function: function creation failed with ocl::runtime_error\n" <<
                "+ ocl::create_function:     reason: " << a_runtime_error.what() << "\n" <<
                "+ ocl::create_function:     message: " << a_runtime_error.additional_message() << "\n";
#endif
            return std::unique_ptr<ocl_function>(nullptr);
        }
        catch (cl::Error& an_error)
        {
#ifdef DEBUG
            std::cerr <<
                "+ ocl::create_function: function creation failed with cl::Error\n" <<
                "+ ocl::create_function:     in function " << an_error.what() << "\n" <<
                "+ ocl::create_function:     because of " << string_of_error_code(an_error.err()) << "\n";
#endif
            return std::unique_ptr<ocl_function>(nullptr);
        }
    }


    class BatchBuilder
    {
    public:
        virtual void finalize() {}
        virtual ~BatchBuilder() {}

        typedef BuildableFunction* value_t;
        virtual void submit(value_t a_function, const std::string& a_build_option = std::string()) = 0;
        virtual void submit(value_t a_function, const char *a_format_string, ...);
    }; // class BatchBuilder


    // This is the most basic implementation of an BatchBuilder.
    // It does not perform any parallelization or implements a
    // sophisticated signal/wait logic.  Still the class may be
    // valuable
    //     (1) for debugging -- in particular the BatchBuilder
    //         itself or
    //     (2) for systems that cannot reliably implement any of the
    //         more complicated schemes.
    class SerialBatchBuilder : public BatchBuilder
    {
    public:
        void submit(value_t a_function, const std::string& a_build_option = std::string());
    }; // class SerialBatchBuilder


    class ThreadedBatchBuilder : public BatchBuilder
    {
    public:
        ThreadedBatchBuilder();
        ~ThreadedBatchBuilder();

        void submit(value_t a_function, const std::string& a_build_option = std::string());

        void finalize();

    private:
        static void build_all_trampoline(ThreadedBatchBuilder* self);

        void build();
        void build_all();

        struct BuildCommand
        {
            BuildCommand() = delete;
            BuildCommand(value_t a_function, const std::string& a_build_option) :
                function(a_function), option(a_build_option) {}

            value_t function;
            std::string option;
        }; // class BuildCommand

        bool run_;
        std::deque<BuildCommand> compile_queue_;
        std::recursive_mutex queue_mutex_;
        std::condition_variable_any queue_not_empty_;
    }; // class ThreadedBatchBuilder

#else

#endif // _OPENCL
} // namespace ocl


#endif // OPENCL_H_INCLUDED

// Local Variables:
// mode: c++
// End:
