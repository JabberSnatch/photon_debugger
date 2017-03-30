#include "profiler.h"
#include "vector.h"
#include "film.h"

//static tools::Profiler gProfiler{};

int main()
{
	for (int j = 0; j < 10; ++j)
	{
		tools::TimeProbe	probe{ gProfiler.GetTimer("TEST") };
		for (int i = 0; i < 1000; ++i)
		{
		}
	}

	maths::Vec3f	vector3_A{};
	maths::Vec3f	vector3_B{ 1.f, 2.f, 3.f };
	maths::Vec3f	vector3_C{ 1.f, 0.f, 0.f };
	maths::Vec3f	vector3_D{ 0.f, 1.f, 0.f };

	maths::Vec3f	vector3_E{ maths::Cross(vector3_C, vector3_D) };
	maths::Vec3f	vector3_F{ vector3_B / 2.f };
	maths::Vec3f	vector3_G{ vector3_B * 2.f };
	maths::Vec3f	vector3_H{ vector3_B * vector3_C };
	maths::Vec3f	vector3_I{ vector3_B * vector3_D };
	maths::Vec3f	vector3_J{ vector3_B * vector3_E };
	float			dot_product{ maths::Dot(vector3_C, vector3_D) };

	raytracer::Film	test_film{ 100, 100 };

	for (int i = 0; i < 100; ++i)
	{
		for (int j = 0; j < 100; ++j)
		{
			test_film.SetPixel({ i / 50.f, j / 50.f, 0.f }, {i, j});
		}
	}

	test_film.WriteToFile("YOLO.png");

	auto timers{ gProfiler.timers() };

	return 0;
}