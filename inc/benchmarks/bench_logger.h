#pragma once
#ifndef __YS_BENCH_LOGGER_HPP__
#define __YS_BENCH_LOGGER_HPP__

#include "globals.h"
#include "profiler.h"
#include "logger.h"
#include <thread>


#define BENCH_HS_MIN 1024
#define BENCH_HS_MAX 4096
#define BENCH_BS_MIN 32
#define BENCH_BS_MAX 4096
#define BENCH_ITER_COUNT 10000
#define BENCH_THREAD_COUNT 4


namespace benchmark {


template <size_t hs, size_t bs>
struct LogIntegrationBench
{
	void operator()(size_t _thread_index, tools::Logger<hs, bs> &_logger)
	{
		_logger.thread_index = _thread_index;
		for (size_t i = 0; i < BENCH_ITER_COUNT; ++i)
		{
			std::string message = "Thread " + std::to_string(_thread_index) + " iter " + std::to_string(i);
			_logger.Log(tools::kChannelGeneral, tools::kLevelInfo, message);
		}

		globals::profiler_aggregate.GrabTimers(globals::profiler);
	}
};

template <size_t hs, size_t bs>
void LoggerRoutine()
{
	const std::string test_suffix = std::to_string(hs) + "x" + std::to_string(bs);

	tools::Logger<hs, bs> yolo;
	yolo.BindPath(tools::kChannelGeneral, "general" + test_suffix + ".log");
	yolo.BindPath(tools::kChannelProfiling, "bench_results.log", false);

	const size_t thread_count = BENCH_THREAD_COUNT;
	yolo.AllowMultipleThreads(thread_count);
	LogIntegrationBench<hs, bs> callable{};
	std::thread threads[thread_count];
	{
		TIMED_SCOPE(LoggingBench);

		for (size_t i = 0; i < thread_count; ++i)
			threads[i] = std::thread(callable, i, std::ref(yolo));
		for (size_t i = 0; i < thread_count; ++i)
			threads[i].join();

		yolo.FlushAll();
	}

	yolo.Log(tools::kChannelProfiling, tools::kLevelInfo, test_suffix);
	globals::profiler_aggregate.GrabTimers(globals::profiler);
	{
		tools::Timer	loggingbench_time = globals::profiler_aggregate.GetTimer("LoggingBench");
		std::stringstream string;
		string << std::endl <<
			loggingbench_time.total_cycles() << " , " << loggingbench_time.total_time();
		yolo.Log(tools::kChannelProfiling, tools::kLevelInfo, string.str());
	}
	{
		tools::Timer	log_time = globals::profiler_aggregate.GetTimer("Logger_Log");
		std::stringstream string;
		string << std::endl <<
			log_time.total_cycles() << " , " << log_time.worst_cycles() << " , " <<
			log_time.average_cycles() << " , " << log_time.best_cycles() << std::endl <<
			log_time.total_time() << " , " << log_time.worst_time() << " , " <<
			log_time.average_time() << " , " << log_time.best_time() << std::endl;
		yolo.Log(tools::kChannelProfiling, tools::kLevelInfo, string.str());
	}

	yolo.FlushAll();
	globals::profiler_aggregate.Clear();
	globals::profiler.Clear();
}


template <size_t hs, size_t bs>
struct inner
{
	static void eval()
	{
		std::cout << "iteration " << hs << "x" << bs << " begin" << std::endl;
		LoggerRoutine<hs, bs>();
		std::cout << "iteration " << hs << "x" << bs << " done" << std::endl;
		inner<hs, bs * 2>::eval();
	}
};
template <size_t hs>
struct inner<hs, BENCH_BS_MAX>
{
	static void eval()
	{
		LoggerRoutine<hs, BENCH_BS_MAX>();
	}
};
template <size_t hs>
struct outer
{
	static void eval()
	{
		inner<hs, BENCH_BS_MIN>::eval();
		outer<hs * 2>::eval();
	}
};
template <>
struct outer<BENCH_HS_MAX>
{
	static void eval()
	{
		inner<BENCH_HS_MAX, BENCH_BS_MIN>::eval();
	}
};
void BenchRoutine()
{
	outer<BENCH_HS_MIN>::eval();
}


} // namespace benchmark


#endif // __YS_BENCH_LOGGER_HPP__
