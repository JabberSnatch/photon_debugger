#ifndef __YS_GLOBALS_HPP__
#define __YS_GLOBALS_HPP__

#include <iostream>
#include "profiler.h"

namespace globals
{

thread_local extern tools::Profiler	profiler;

} // namespace globals

#define TIMED_SCOPE(name) TIMED_SCOPE_FULL(name, globals::profiler)

#ifdef _WIN32

#ifdef _DEBUG
#define YS_DEBUG
#define LOG(message) std::cout << message << std::endl
#else
#define LOG(message) {}
#endif

#endif // _WIN32

#endif // __YS_GLOBALS_HPP__
