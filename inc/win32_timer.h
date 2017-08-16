#ifndef __YS_WIN32_TIMER_HPP__
#define __YS_WIN32_TIMER_HPP__

#include <windows.h>
#include <intrin.h>
#include "maths.h"

namespace tools
{


class Timer final
{
public:
	using Duration_t = long long;
	using CycleCount_t = DWORD64;
	using Frequency_t = LARGE_INTEGER;

	Timer(char const *_name) :
		name_{ _name }
	{}

	inline void
	Add(Duration_t const &_time, CycleCount_t const &_cycles)
	{
		total_time_ += _time;
		if (_time < best_time_) best_time_ = _time;
		if (_time > worst_time_) worst_time_ = _time;
		total_cycles_ += _cycles;
		if (_cycles < best_cycles_) best_cycles_ = _cycles;
		if (_cycles > worst_cycles_) worst_cycles_ = _cycles;
		++call_count_;
	}

	inline void
	Add(Timer const &_other)
	{
		total_time_ += _other.total_time_;
		best_time_ = maths::Min(best_time_, _other.best_time_);
		total_cycles_ += _other.total_cycles_;
		best_cycles_ = maths::Min(best_cycles_, _other.best_cycles_);
		call_count_ += _other.call_count_;
	}

	static inline Frequency_t
	timer_frequency()
	{
		Frequency_t freq;
		QueryPerformanceFrequency(&freq);
		return freq;
	}
	static const Frequency_t		kTimerFrequency;


	inline char const *
	name() const { return name_; }
	inline unsigned long long
	call_count() const { return call_count_; }
	inline double
	total_time() const { return total_time_ / (double)kTimerFrequency.QuadPart; }
	inline double
	best_time() const { return best_time_ / (double)kTimerFrequency.QuadPart; }
	inline double
	worst_time() const { return worst_time_ / (double)kTimerFrequency.QuadPart; }
	inline double
	average_time() const { return total_time() / call_count_; }
	inline CycleCount_t
	total_cycles() const { return total_cycles_; }
	inline CycleCount_t
	best_cycles() const { return best_cycles_; }
	inline CycleCount_t
	worst_cycles() const { return worst_cycles_; }
	inline CycleCount_t
	average_cycles() const { return total_cycles() / call_count_; }
	inline Duration_t
	total_ticks() const { return total_time_; }
	inline Duration_t
	best_ticks() const { return best_time_; }
	inline Duration_t
	average_ticks() const { return total_time_ / call_count_; }

private:
	char const			*name_{};
	Duration_t			total_time_{ 0 };
	Duration_t			best_time_{ _I64_MAX };
	Duration_t			worst_time_{ 0 };
	CycleCount_t		total_cycles_{ 0 };
	CycleCount_t		best_cycles_{ _I64_MAX };
	CycleCount_t		worst_cycles_{ 0 };
	unsigned long long	call_count_{ 0 };
};


class TimeProbe final
{
public:
	using TimePoint_t = LARGE_INTEGER;

	TimeProbe() = delete;
	TimeProbe(Timer &_timer) :
		timer_{ _timer }
	{}

	~TimeProbe()
	{
		TimePoint_t				end{ now() };
		Timer::CycleCount_t		last_cycle{ __rdtsc() };
		timer_.Add(end.QuadPart - start_.QuadPart, last_cycle - first_cycle_);
	}

	TimeProbe &operator=(TimeProbe const &) = delete;

	static inline TimePoint_t
	now()
	{
		TimePoint_t	value;
		QueryPerformanceCounter(&value);
		return value;
	}

private:
	Timer					&timer_;
	TimePoint_t				start_{ now() };
	Timer::CycleCount_t		first_cycle_{ __rdtsc() };
};


} // namespace tools

#endif // __YS_WIN32_TIMER_HPP__
