#include <typeinfo>
#include <iostream>
#include <sstream>

#include "boost/filesystem.hpp"

#include "algorithms.h"
#include "primes.h"
#include "api/input_processor.h"
#include "benchmarks/bench_logger.h"
#include "core/logger.h"
#include "core/profiler.h"
#include "core/rng.h"

int main(int argc, char *argv[])
{
#if 0
	constexpr uint64_t count = 10;
	using CountSizedArray_t = std::array<uint64_t, count>;
	using DoubleSizedArray_t = std::array<uint64_t, count * 2>;
	//constexpr uint64_t prime0 = algo::prime<count>();
	//constexpr CountSizedArray_t first_primes = algo::get_primes<count>();
	//constexpr CountSizedArray_t next_primes = algo::get_primes<count, count>(first_primes);
	//constexpr DoubleSizedArray_t all_primes = algo::get_primes<count * 2>();
	//constexpr CountSizedArray_t primes = algo::prime_rec::get_primes<count>();
	//constexpr uint64_t prime = algo::compute_next_prime();
	//constexpr uint64_t prime3 = algo::compute_next_prime(2, 3, 5, 7, 11, 13);

	algo::func_struct<3> func{};
	algo::TestFunctionType copy{ func };
	uint64_t result = copy(3);
	using container = algo::FunctionContainer<algo::func_struct, 3>;
	using container_type = container::type;
	algo::TestFunctionType yolo{ container_type{} };
	result = container::value(3);
	using value_container = std::array<uint64_t, count>;
	using functor_container = std::array<algo::func_struct_container<0>::type, count>;
	//constexpr value_container values = algo::build_array<count, algo::UintContainer>();
	constexpr value_container values = algo::build_array_rec<count, algo::UintContainer>::value;
	constexpr functor_container functors = algo::build_array_rec<count, algo::func_struct_container>::value;
	//using test = algo::FunctionContainer<3, std::function<uint64_t(uint64_t)>, &algo::func>;
	//uint64_t result = test::value(3);

	return 0;
#endif

	uint32_t input{ 65421356u }, output{ 0u }, inverse{ 0u };
	uint64_t seed{ 345981108932u };
	core::RNG rng{ seed };
	{
		TIMED_SCOPE(RNG_Bench);
		for (uint64_t i = 0; i < 10000000u; ++i)
		{
			rng();
		}
	}

	globals::logger.BindPath(tools::kChannelGeneral, "general.log");
	globals::logger.BindPath(tools::kChannelProfiling, "profiling.log");
	globals::logger.BindPath(tools::kChannelParsing, "parsing.log");

	if (argc > 1 && false)
	{
		std::string const input_file{ argv[1] };
		std::string const absolute_path = boost::filesystem::absolute(input_file).string();
		if (boost::filesystem::exists(absolute_path))
		{
			api::ProcessInputFile(absolute_path);
		}
		else
		{
			std::cout << "No such file as " << absolute_path << std::endl;
		}
	}
	else
	{
		std::cout << "Please provide a path to an input file as first argument." << std::endl;
	}

	globals::profiler_aggregate.GrabTimers(globals::profiler);
	{
		tools::Profiler::TimerTable_t timers = globals::profiler_aggregate.timers();
		for (auto it = timers.begin(); it != timers.end(); ++it)
		{
			tools::Timer const &timer = it->second;
			char const *call_string = (timer.call_count() > 1) ? " calls. " : " call. ";
			std::stringstream string;
			string <<
				timer.name() << " : " << std::endl <<
				timer.call_count() << call_string << std::endl <<
				"real time : [ " <<
				timer.total_time() << " : " << timer.worst_time() << " + " <<
				timer.average_time() << " - " <<
				timer.best_time() << " ] " << std::endl <<
				"cycles :    [ " <<
				timer.total_cycles() << " : " << timer.worst_cycles() << " + " <<
				timer.average_cycles() << " - " <<
				timer.best_cycles() << " ]";
			globals::logger.Log(tools::kChannelProfiling, tools::kLevelInfo, string.str());
		}
	}

	globals::logger.FlushAll();

	system("pause");
	return 0;
}