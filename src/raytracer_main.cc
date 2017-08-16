#include "vector.h"
#include "point.h"
#include "bounds.h"
#include "matrix.h"
#include "quaternion.h"
#include "transform.h"
#include "film.h"
#include "globals.h"
#include "algorithms.h"
#include "redecimal.h"
#include "surface_interaction.h"
#include "camera.h"
#include "shapes/sphere.h"
#include "triangle_mesh.h"
#include "shapes/triangle.h"
#include "logger.h"
#include "profiler.h"
#include "input_processor.h"
#include "bvh_accelerator.h"
#include "primitive.h"
#include "memory_region.h"

#include <typeinfo>
#include <iostream>
#include <sstream>
#include "boost/filesystem.hpp"

#include "benchmarks/bench_logger.h"

int main(int argc, char *argv[])
{
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