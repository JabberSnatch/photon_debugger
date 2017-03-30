#include "profiler.h"
#include "vector.h"

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

	maths::Vec3f	vector3_A{};
	maths::Vec3f	vector3_B{ 1.f, 2.f, 3.f };
	maths::Vec3f	vector3_C{ 1.f, 0.f, 0.f };
	maths::Vec3f	vector3_D{ 0.f, 1.f, 0.f };

	maths::Vec3f	vector3_E{ maths::Cross(vector3_C, vector3_D) };
	float			dot_product{ maths::Dot(vector3_C, vector3_D) };

	return 0;
}