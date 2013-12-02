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
#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdint>

#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>   // times()
#endif

#include <time.h>        // clock_t, clock(); timespec, clock_gettime()

#ifdef WIN32
#include <Windows.h>
#endif


namespace timer
{
    class Interface
    {
    public:
        // The constructor of each derived class shall call start().

        virtual ~Interface() {}

        // Reset counter to zero and start a new measurement.
        virtual void start() = 0;

        // Stop measurement.  After stop() all calls to value() return
        // the same result.
        virtual void stop() = 0;

        // Continue to measure after counter has been stopped.  Do not
        // reset the counter, instead accumulate.
        virtual void restart() = 0;

        // Answer the time passed between calls to start() and stop()
        // or restart() and stop() in seconds.
        virtual double value() const = 0;
    }; // class Interface


    class StandardWallClock : public Interface
    {
    public:
        StandardWallClock();

        void start();
        void stop();
        void restart();

        double value() const;

    private:
        clock_t value_;
        clock_t start_;
        clock_t stop_;
    }; // class StandardWallClock


#ifdef HAVE_CLOCK_GETTIME
    class RealTimeWallClock : public Interface
    {
    public:
        RealTimeWallClock();

        void restart();
        void start();
        void stop();

        double value() const;

    private:
        typedef std::uint64_t value_t;
        // typedef decltype(struct timespec . tv_sec) seconds_t;
        // typedef decltype(struct timespec . tv_nsec) nano_seconds_t;

        value_t value_;         // unit: nano seconds
        timespec start_;
        timespec stop_;
    }; // class RealTimeWallClock

#define WallClock RealTimeWallClock
#else
#define WallClock StandardWallClock
#endif // HAVE_CLOCK_GETTIME


#ifdef WIN32
    class ProcessorTime : public Interface
    {
    public:
        ProcessorTime();

        void start();
        void stop();
        void restart();

    protected:
        ULONGLONG user_value_;
        ULONGLONG system_value_;
        ULONGLONG start_user_value_;
        ULONGLONG start_system_value_;
        ULONGLONG start_idle_value_;
    }; // class ProcessorTime

#elif defined(HAVE_SYS_TIMES_H)

    class ProcessorTime : public Interface
    {
    public:
        ProcessorTime();

        void start();
        void stop();
        void restart();

    protected:
        clock_t user_value_;
        clock_t system_value_;
        tms start_;
        tms stop_;
    }; // class ProcessorTime

#else

    // Null class -- does nothing
    class ProcessorTime : public Interface
    {
    public:
        void start() {}
        void stop() {}
        void restart() {}
    }; // class ProcessorTime
#endif

    class UserTime : public ProcessorTime
    {
    public:
        double value() const;
    }; // class UserTime


    class SystemTime : public ProcessorTime
    {
    public:
        double value() const;
    }; // class SystemTime
} // namespace timer


#endif // TIMER_H_INCLUDED

// Local Variables:
// mode: c++
// End:
