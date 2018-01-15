#include <typeinfo>
#include <iostream>
#include <sstream>

#include "boost/filesystem.hpp"

#include "algorithms.h"
#include "primes.h"
#include "api/input_processor.h"
#include "api/translation_state.h"
#include "benchmarks/bench_logger.h"
#include "core/logger.h"
#include "core/profiler.h"
#include "core/rng.h"
#include "raytracer/samplers/halton_sampler.h"


void flush_profiler()
{
	globals::profiler_aggregate.GrabTimers(globals::profiler);
	{
		tools::Profiler::TimerTable_t timers = globals::profiler_aggregate.timers();
		for (auto it = timers.begin(); it != timers.end(); ++it)
		{
			tools::Timer const &timer = it->second;
			char const *call_string = (timer.call_count() > 1) ? " calls. " : " call. ";
			std::stringstream string;
			string << "\n" << 
				timer.name() << " : " << std::endl <<
				"	" << timer.call_count() << call_string << std::endl <<
				"	real time : [ " <<
				timer.total_time() << " : " << timer.worst_time() << " + " <<
				timer.average_time() << " - " <<
				timer.best_time() << " ] " << std::endl <<
				"	cycles :    [ " <<
				timer.total_cycles() << " : " << timer.worst_cycles() << " + " <<
				timer.average_cycles() << " - " <<
				timer.best_cycles() << " ]";
			globals::logger.Log(tools::kChannelProfiling, tools::kLevelInfo, string.str());
		}
	}
}


void flush_logger()
{
	globals::logger.FlushAll();
}


int main(int argc, char *argv[])
{
	// TODO: Make this depend on the file name
	globals::logger.BindPath(tools::kChannelGeneral, "general.log");
	globals::logger.BindPath(tools::kChannelProfiling, "profiling.log");
	globals::logger.BindPath(tools::kChannelParsing, "parsing.log");

	// TODO: Implement $ ./RayTracer.exe _JUNK/magikarp.txt
	std::string absolute_path{};
	bool interactive_mode = false;
	if (argc > 1)
	{
		for (int i = 1; i < argc; ++i)
		{
			std::string const arg{ argv[i] };
			if (arg == "--interactive")
			{
				interactive_mode = true;
			}
			else
			{
				absolute_path = boost::filesystem::absolute(arg).generic_string();
				if (!boost::filesystem::exists(absolute_path))
				{
					std::cout << "No such file as " << absolute_path << std::endl;
				}
			}
		}
	}
	else
	{
		std::cout << "Please provide a path to an input file as first argument." << std::endl;
	}

	{
		api::TranslationState translation_state{};
		if (interactive_mode)
		{
			std::string input_string{};
			bool running = true;
			while (running)
			{
				std::cout << "sup\'>> ";
				std::cin >> input_string;
				if (input_string == "render")
				{
					translation_state.ResetResourceCounters();
					if (!absolute_path.empty() && boost::filesystem::exists(absolute_path))
					{
						api::ProcessInputFile(absolute_path, translation_state);
					}
					else
					{
						std::cout << "No input file" << std::endl;
					}
					//
					if (translation_state.render_context().GoodForRender())
					{
						translation_state.render_context().RenderAndWrite(translation_state.output_path());
					}
					else
					{
						std::cout << "Error in input file" << std::endl;
					}
					flush_profiler();
					flush_logger();
				}
				else if (input_string == "exit")
				{
					running = false;
				}
				else
				{
					std::cout << "Unknown command" << std::endl;
				}
			}
		}
		else
		{
			if (!absolute_path.empty() && boost::filesystem::exists(absolute_path))
			{
				api::ProcessInputFile(absolute_path, translation_state);
			}
			else
			{
				std::cout << "No input file" << std::endl;
			}
			//
			if (translation_state.render_context().GoodForRender())
			{
				translation_state.render_context().RenderAndWrite(translation_state.output_path());
			}
			else
			{
				std::cout << "Error in input file" << std::endl;
			}
			flush_profiler();
			flush_logger();
		}
	}


	
	system("pause");
	return 0;
}
