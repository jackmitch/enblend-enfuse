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


#include "timer.h"


namespace timer
{
    inline static double
    clock_in_seconds(clock_t a_clock_time)
    {
        return static_cast<double>(a_clock_time) / static_cast<double>(CLOCKS_PER_SEC);
    }


    StandardWallClock::StandardWallClock()
    {
        start();
    }


    void
    StandardWallClock::start()
    {
        value_ = clock_t();
        start_ = clock();
        stop_ = start_;
    }


    void
    StandardWallClock::stop()
    {
        stop_ = clock();
        value_ += stop_ - start_;
        start_ = stop_;
    }


    void
    StandardWallClock::restart()
    {
        start_ = clock();
        stop_ = start_;
    }


    double
    StandardWallClock::value() const
    {
        return clock_in_seconds(value_);
    }


#ifdef HAVE_CLOCK_GETTIME
    inline static std::uint64_t
    timespec_in_nano_seconds(const struct timespec& a_timespec)
    {
        return
            static_cast<std::uint64_t>(a_timespec.tv_sec) * 1000000000UL +
            static_cast<std::uint64_t>(a_timespec.tv_nsec);
    }


    RealTimeWallClock::RealTimeWallClock()
    {
        start();
    }


    void
    RealTimeWallClock::start()
    {
        value_ = value_t();
        clock_gettime(CLOCK_REALTIME, &start_);
        stop_.tv_sec = start_.tv_sec;
        stop_.tv_nsec = long();
    }


    void
    RealTimeWallClock::stop()
    {
        clock_gettime(CLOCK_REALTIME, &stop_);
        value_ += timespec_in_nano_seconds(stop_) - timespec_in_nano_seconds(start_);
        start_ = stop_;
    }


    void
    RealTimeWallClock::restart()
    {
        clock_gettime(CLOCK_REALTIME, &start_);
        stop_ = start_;
    }


    double
    RealTimeWallClock::value() const
    {
        return static_cast<double>(value_) * 1e-9;
    }
#endif // HAVE_CLOCK_GETTIME


#ifdef WIN32
    inline static ULONGLONG
    filetime_in_100nanoseconds(const _FILETIME& a_filetime)
    {
        return
            (static_cast<ULONGLONG>(a_filetime.dwHighDateTime) << 32) +
             static_cast<ULONGLONG>(a_filetime.dwLowDateTime);
    }


    ProcessorTime::ProcessorTime()
    {
        start();
    }


    void
    ProcessorTime::start()
    {
        _FILETIME idle;
        _FILETIME kernel;
        _FILETIME user;

        user_value_ = ULONGLONG();
        system_value_ = ULONGLONG();

        GetSystemTimes(&idle, &kernel, &user);
        start_idle_value_ = filetime_in_100nanoseconds(idle);
        start_system_value_ = filetime_in_100nanoseconds(kernel);
        start_user_value_ = filetime_in_100nanoseconds(user);
    }


    void
    ProcessorTime::stop()
    {
        _FILETIME idle;
        _FILETIME kernel;
        _FILETIME user;

        GetSystemTimes(&idle, &kernel, &user);
        const ULONGLONG stop_idle_value = filetime_in_100nanoseconds(idle);
        const ULONGLONG stop_system_value = filetime_in_100nanoseconds(kernel);
        const ULONGLONG stop_user_value = filetime_in_100nanoseconds(user);

        system_value_ = stop_system_value - start_system_value_ - (stop_idle_value - start_idle_value_);
        user_value_ = stop_user_value - start_user_value_;
        start_idle_value_ = stop_idle_value;
        start_system_value_ = stop_system_value;
        start_user_value_ = stop_user_value;
    }


    void
    ProcessorTime::restart()
    {
        _FILETIME idle;
        _FILETIME kernel;
        _FILETIME user;

        GetSystemTimes(&idle, &kernel, &user);
        start_idle_value_ = filetime_in_100nanoseconds(idle);
        start_system_value_ = filetime_in_100nanoseconds(kernel);
        start_user_value_ = filetime_in_100nanoseconds(user);
    }


    double
    UserTime::value() const
    {
        return user_value_ / 1.0E7;
    }


    double
    SystemTime::value() const
    {
        return system_value_ / 1.0E7;
    }

#elif defined(HAVE_SYS_TIMES_H)

    ProcessorTime::ProcessorTime()
    {
        start();
    }


    void
    ProcessorTime::start()
    {
        user_value_ = clock_t();
        system_value_ = clock_t();
        times(&start_);
        stop_.tms_utime = start_.tms_utime;
        stop_.tms_stime = start_.tms_stime;
    }


    void
    ProcessorTime::stop()
    {
        times(&stop_);
        user_value_ = stop_.tms_utime - start_.tms_utime;
        system_value_ = stop_.tms_stime - start_.tms_stime;
        start_.tms_utime = stop_.tms_utime;
        start_.tms_stime = stop_.tms_stime;
    }


    void
    ProcessorTime::restart()
    {
        times(&start_);
        stop_.tms_utime = start_.tms_utime;
        stop_.tms_stime = start_.tms_stime;
    }


    double
    UserTime::value() const
    {
        return clock_in_seconds(user_value_);
    }


    double
    SystemTime::value() const
    {
        return clock_in_seconds(system_value_);
    }

#else

    double
    UserTime::value() const
    {
        return 0.0;
    }


    double
    SystemTime::value() const
    {
        return 0.0;
    }
#endif
} // namespace timer
