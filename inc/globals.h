#pragma once
#ifndef __YS_GLOBALS_HPP__
#define __YS_GLOBALS_HPP__

#include <iostream>
//#include "profiler.h"
//#include "logger.h"

namespace tools {
class Profiler;
template <size_t hs, size_t bs> class Logger;
}

namespace globals
{

extern tools::Profiler				profiler_aggregate;
thread_local extern tools::Profiler	profiler;
extern tools::Logger<256, 32>		logger;

} // namespace globals

#define TIMED_SCOPE_FULL(name, profiler) tools::TimeProbe name##_probe {profiler.GetTimer(#name)}
#define TIMED_SCOPE(name) TIMED_SCOPE_FULL(name, globals::profiler)


#define LOG_FULL(channel, level, message, logger) logger.Log(channel, level, message)
#ifndef YS_NO_LOGS
#define LOG(channel, level, message) LOG_FULL(channel, level, message, globals::logger)
#define LOG_DEBUG(channel, message) LOG(channel, tools:kLevelDebug, message)
#define LOG_INFO(channel, message) LOG(channel, tools::kLevelInfo, message)
#define LOG_ERROR(channel, message) LOG(channel, tools::kLevelError, message)
#define LOG_WARNING(channel, message) LOG(channel, tools::kLevelWarning, message)
#else
#define LOG(channel, level, message) {}
#define LOG_DEBUG(channel, message) {}
#define LOG_INFO(channel, message) {}
#define LOG_ERROR(channel, message) {}
#define LOG_WARNING(channel, message) {}
#endif

#endif // __YS_GLOBALS_HPP__
