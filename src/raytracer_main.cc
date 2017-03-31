#include "profiler.h"
#include "vector.h"
#include "film.h"
#include "globals.h"

//static tools::Profiler gProfiler{};

int main()
{
	for (int j = 0; j < 10; ++j)
	{
		tools::TimeProbe	probe{ globals::profiler.GetTimer("TEST") };
		for (int i = 0; i < 1000; ++i)
		{
		}
	}

	maths::Vec3f	vector3_A{};
	maths::Vec3f	vector3_B{ 1.f, 2.f, 3.f };
	maths::Vec3f	vector3_C{ 1.f, 0.f, 0.f };
	maths::Vec3f	vector3_D{ 0.f, 1.f, 0.f };

	maths::Vec3f	vector3_E{ maths::vector::Cross(vector3_C, vector3_D) };
	maths::Vec3f	vector3_F{ vector3_B / 2.f };
	maths::Vec3f	vector3_G{ vector3_B * 2.f };
	maths::Vec3f	vector3_H{ vector3_B * vector3_C };
	maths::Vec3f	vector3_I{ vector3_B * vector3_D };
	maths::Vec3f	vector3_J{ vector3_B * vector3_E };
	float			dot_product{ maths::vector::Dot(vector3_C, vector3_D) };

	raytracer::Film	test_film{ 100, 100 };

	for (int32_t i = 0; i < 100; ++i)
	{
		for (int32_t j = 0; j < 100; ++j)
		{
			test_film.SetPixel({ i / 50.f, j / 50.f, 0.f }, {i, j});
		}
	}

	test_film.WriteToFile("YOLO.png");

	auto timers{ globals::profiler.timers() };

	maths::Vector<float, 4> A4{ 1.f, 2.f, 3.f, 4.f };
	maths::Vector<float, 4> B4{ 1.f, 2.f, 3.f, 4.f };
	maths::Vector<float, 3> A3{ 1.f, 2.f, 3.f };
	maths::Vector<float, 3> B3{ 1.f, 2.f, 3.f };

	A4 += B4;
	A3 += B3;

	//maths::Vector<float, 4> C4{ A4 + B4 };
	//maths::Vector<float, 3> C3{ A3 + B3 };

	return 0;
}