#include "globals.h"

namespace globals
{

thread_local tools::Profiler							profiler{};
tools::Logger<logger.kHistorySize, logger.kBufferSize>	logger{};

} // namespace globals

