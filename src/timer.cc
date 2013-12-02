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
} // namespace timer
