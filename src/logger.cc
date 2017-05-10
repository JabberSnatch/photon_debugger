#include "logger.h"

namespace tools {

LogEntry::Clock_t::time_point	LogEntry::start_time = LogEntry::Clock_t::now();

} // namespace tools
