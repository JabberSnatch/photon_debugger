#include "globals.h"

#include "core/logger.h"
#include "core/profiler.h"

namespace globals
{

tools::Profiler											profiler_aggregate{};
thread_local tools::Profiler							profiler{};
tools::Logger<logger.kHistorySize, logger.kBufferSize>	logger{};

} // namespace globals

