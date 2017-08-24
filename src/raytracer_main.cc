#include "maths/vector.h"
#include "maths/point.h"
#include "maths/bounds.h"
#include "maths/matrix.h"
#include "maths/quaternion.h"
#include "maths/transform.h"
#include "raytracer/film.h"
#include "globals.h"
#include "maths/redecimal.h"
#include "raytracer/surface_interaction.h"
#include "raytracer/camera.h"
#include "raytracer/shapes/sphere.h"
#include "raytracer/triangle_mesh.h"
#include "raytracer/shapes/triangle.h"
#include "core/logger.h"
#include "core/profiler.h"
#include "api/input_processor.h"
#include "raytracer/bvh_accelerator.h"
#include "raytracer/primitive.h"
#include "core/memory_region.h"
#include "algorithms.h"
#include "primes.h"

#include <typeinfo>
#include <iostream>
#include <sstream>
#include "boost/filesystem.hpp"

#include "benchmarks/bench_logger.h"

int main(int argc, char *argv[])
{
	constexpr uint64_t index = 128;
	using IndexSizedArray_t = std::array<uint64_t, index>;
	constexpr IndexSizedArray_t primes = algo::primes::get_array<index>();

	return 0;

	globals::logger.BindPath(tools::kChannelGeneral, "general.log");
	globals::logger.BindPath(tools::kChannelProfiling, "profiling.log");
	globals::logger.BindPath(tools::kChannelParsing, "parsing.log");

	if (argc > 1)
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