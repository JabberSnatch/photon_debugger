#include "profiler.h"

int main()
{
	for (int j = 0; j < 10; ++j)
	{
		tools::TimeProbe	probe{ gProfiler.GetTimer("TEST") };
		for (int i = 0; i < 1000; ++i)
		{
		}
	}

	auto timers{ gProfiler.timers() };

	return 0;
}