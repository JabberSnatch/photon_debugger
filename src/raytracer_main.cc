#include <csignal>
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


void flush_profiler();
void flush_logger();

void render(std::string const &_path, api::TranslationState &_translation_state)
{
	_translation_state.ResetResourceCounters();
	if (!_path.empty() && boost::filesystem::exists(_path))
	{
		api::ProcessInputFile(_path, _translation_state);
	}
	else
	{
		std::cout << "No input file" << std::endl;
	}
	//
	if (_translation_state.render_context().GoodForRender())
	{
		_translation_state.render_context().RenderAndWrite(_translation_state.output_path());
	}
	else
	{
		std::cout << "Error in input file" << std::endl;
	}
	flush_profiler();
	flush_logger();
}


void signal_handler(int signal)
{
	std::cout << "unexpected error, shutting down.." << std::endl;
	globals::logger.~Logger();
}


int main(int argc, char *argv[])
{
	std::signal(SIGSEGV, signal_handler);
	std::signal(SIGTERM, signal_handler);
	std::signal(SIGABRT, signal_handler);

	globals::logger.BindPath(tools::kChannelGeneral, "general.log");
	globals::logger.BindPath(tools::kChannelProfiling, "profiling.log");
	globals::logger.BindPath(tools::kChannelParsing, "parsing.log");

	// TODO: make logging output dependant on output name
	// TODO: add planes (and quads)
	// TODO: add support for c4d files
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
					render(absolute_path, translation_state);
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
			render(absolute_path, translation_state);
		}
	}

	system("pause");

	std::signal(SIGSEGV, SIG_DFL);
	std::signal(SIGTERM, SIG_DFL);
	std::signal(SIGABRT, SIG_DFL);

	return 0;
}


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
