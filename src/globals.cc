#include "globals.h"

#include "logger.h"
#include "profiler.h"

namespace globals
{

tools::Profiler											profiler_aggregate{};
thread_local tools::Profiler							profiler{};
tools::Logger<logger.kHistorySize, logger.kBufferSize>	logger{};

} // namespace globals

