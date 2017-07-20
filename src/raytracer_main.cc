#include "vector.h"
#include "point.h"
#include "bounds.h"
#include "matrix.h"
#include "quaternion.h"
#include "transform.h"
#include "film.h"
#include "globals.h"
#include "algorithms.h"
#include "redecimal.h"
#include "surface_interaction.h"
#include "camera.h"
#include "shapes/sphere.h"
#include "triangle_mesh.h"
#include "shapes/triangle.h"
#include "logger.h"
#include "profiler.h"
#include "input_processor.h"
#include "bvh_accelerator.h"
#include "primitive.h"
#include "memory_region.h"

#include <typeinfo>
#include <iostream>
#include <sstream>

#include "benchmarks/bench_logger.h"

int main()
{
	api::ProcessInputFile("TestScene.txt");

	globals::logger.BindPath(tools::kChannelGeneral, "general.log");
	globals::logger.BindPath(tools::kChannelProfiling, "profiling.log");

	core::MemoryRegion<>	main_region;

	maths::Point3f		swizzle_test{ 1._d, 2._d, 3._d };
	maths::Point3f		swizzledA{ maths::Swizzle(swizzle_test, 2u, 0u, 1u) };

	raytracer::Film		*film_35mm = main_region.Alloc<raytracer::Film>();
	new (film_35mm) raytracer::Film{ 1000, 1000, 0.035_d };

	raytracer::Camera	*camera = main_region.Alloc<raytracer::Camera>();
	new (camera) raytracer::Camera{
		{0._d, -5._d, 2.5_d}, {0._d, 1._d, 2.5_d}, {0._d, 0._d, 1._d},
		60._d, film_35mm
	};
	//raytracer::Film		film_35mm{ 1000, 1000, 0.035_d };
	//raytracer::Camera	camera{
	//	film_35mm,
	//	{0._d, -5._d, 2.5_d}, {0._d, 1._d, 2.5_d}, {0._d, 0._d, 1._d},
	//	60._d
	//};
	maths::Transform	little_sphere_transform = maths::Translate({ 0._d, 1._d, .0_d }) * maths::Scale(0.02_d, 0.02_d, 0.02_d) * maths::RotateX(180._d);// *maths::RotateZ(90._d);
	raytracer::Sphere	little_sphere{ little_sphere_transform, false, .5_d, -.5_d, .2_d, 45.3_d };
	maths::Transform	big_sphere_transform = maths::Translate({ 0._d, 1._d, -500._d });
	raytracer::Sphere	big_sphere{ big_sphere_transform, false, 500._d, -500._d, 500._d, 360._d };

#ifdef TEST_TRIANGLE
	std::vector<int32_t>		indices{ 0, 1, 2 };
	std::vector<maths::Point3f>	vertices{ {0._d, 0._d, 1._d}, {-1._d, 0._d, 0._d}, {1._d, 0._d, 0._d} };
	std::vector<maths::Norm3f>	normals{ {0._d, -1._d, 0._d}, {.0_d, -1._d, 0._d}, {0._d, -1._d, 0._d} };
	raytracer::TriangleMesh		test_mesh(little_sphere_transform, 1, indices, vertices, &normals);
	raytracer::Triangle			test_triangle(little_sphere_transform, false, test_mesh, 0);

	maths::Ray					test_ray{ {0._d, 0._d, .5_d}, {0._d, 1._d, 0._d}, maths::infinity<maths::Decimal>, 0._d };
	maths::Decimal				t_hit;
	raytracer::SurfaceInteraction	hit_info;
	bool result = test_triangle.Intersect(test_ray, t_hit, hit_info);
	YS_ASSERT(result);
#endif // TEST_TRIANGLE

	core::MemoryRegion<>	shapes_region;
	core::MemoryRegion<>	primitives_region;
	std::vector<raytracer::Primitive*>	primitives;
	
#if 1
	raytracer::TriangleMesh	test_mesh =
		raytracer::ReadTriangleMeshFromFile("Leon/Leon.dae", little_sphere_transform);
	raytracer::Triangle				*test_triangles =
		shapes_region.Alloc<raytracer::Triangle>(test_mesh.triangle_count);
	raytracer::GeometryPrimitive	*test_primitives =
		primitives_region.Alloc<raytracer::GeometryPrimitive>(test_mesh.triangle_count);
	for (int32_t i = 0; i < test_mesh.triangle_count; ++i)
	{
		new (test_triangles + i) raytracer::Triangle(little_sphere_transform, false, test_mesh, i);
		new (test_primitives + i) raytracer::GeometryPrimitive(test_triangles[i]);
		primitives.push_back(&test_primitives[i]);
	}
#endif 

	raytracer::Primitive	*sphere = primitives_region.Alloc<raytracer::GeometryPrimitive>();
	new (sphere) raytracer::GeometryPrimitive(big_sphere);
	primitives.emplace_back(sphere);

	raytracer::BvhAccelerator	accelerator{ primitives, 1 };
	raytracer::Camera::Scene	scene = { &accelerator };

	camera->Expose(scene, 0._d);
	camera->WriteToFile("squall.png");

	maths::FloatBitsMapper start{ 34.23f };
	maths::FloatBitsMapper end{ maths::NextDecimalUp(start.value) };
	maths::FloatBitsMapper delta{ end.value - start.value };
	maths::FloatBitsMapper epsilon_machine{ 0x33800000u };	// 2^-24
	maths::FloatBitsMapper epsilon{ 0x34000000u };			// 2^-23
	maths::FloatBitsMapper verif{ start.value * epsilon.value };
	maths::FloatBitsMapper verif2{ start.value * (1.f + epsilon_machine.value) };
	maths::FloatBitsMapper verif3{ start.value * (1.f - epsilon_machine.value) };
	maths::FloatBitsMapper one_eps{ 1.f + epsilon.value };

	constexpr std::array<float, 5> arr_A{ algo::fill<5>::apply(0.4f) };

	maths::Vec3f								v3_default_ctor{};
	maths::Vec3f								v3_fill_ctor(0.4_d);
	constexpr maths::Vec3f						v3_initlist_ctor{ 1._d, 2._d, 3._d };
	maths::Decimal								v3_bracket{ v3_initlist_ctor[0] };
	maths::Vector<maths::Decimal, 2>			v2_fill_ctor(1._d);
	constexpr maths::Vector<maths::Decimal, 2>	v2_initlist_ctor{ 1._d, 2._d };
	maths::Vec3f								v3_extension_ctor(v2_initlist_ctor, 3._d);
	bool										v3_hasnans = v3_initlist_ctor.HasNaNs();
	bool										v3_equal = maths::Vec3f{ 1._d, 2._d, 3._d } == v3_initlist_ctor;
	maths::Vector<maths::Decimal, 5>			vn_fill_ctor(0.4_d);

	constexpr maths::Norm3f						n3_initlist_ctor{ 1._d, 2._d, 3._d };
	constexpr maths::Decimal					v3_n3_dot = maths::Dot(v3_initlist_ctor, n3_initlist_ctor);

	constexpr uint32_t	vn_size{ maths::Vector<float, 5>::size };

	maths::Mat4x4f								m44_default_ctor{};
	constexpr maths::Mat4x4f					m44_fill_ctor(0.4_d);

	maths::Point3f								p3_default_ctor{};
	maths::Point3f								p3_fill_ctor(0.4_d);
	constexpr maths::Point3f					p3_initlist_ctor{ 1._d, 2._d, 3._d };
	constexpr maths::Vec3f						v3_p3_conversion{ (maths::Vec3f)p3_initlist_ctor };

	maths::Quaternion							quat_angleaxis_ctor(v3_p3_conversion, 0.4_d);
	maths::Transform							transf_fromquat((maths::Transform)quat_angleaxis_ctor);

	constexpr maths::Decimal					decim_lowest{ maths::lowest_value<maths::Decimal> };
	maths::Bounds3f								b3_default_ctor{};
	constexpr maths::Bounds3f					b3_point_ctor{ p3_initlist_ctor };
	maths::Bounds3f								b3_points_ctor{ p3_initlist_ctor, p3_initlist_ctor };
	maths::Point3f								p3_b3brackets{ b3_point_ctor[1] };

	maths::REDecimal							redecim_value_ctor{ 0.4_d };

	maths::Decimal	muc = maths::NextDecimalDown(0._d);
	maths::REDecimal redecim_a{ 2._d }, redecim_b{ 4.3_d }, redecim_c{ -7.8_d }, redecim_t0, redecim_t1;
	if (maths::Quadratic(redecim_a, redecim_b, redecim_c, redecim_t0, redecim_t1))
		int i = 0;
	maths::Decimal	decim_t0, decim_t1;
	if (maths::Quadratic(2._d, 4.3_d, -7.8_d, decim_t0, decim_t1))
		int i = 0;

	maths::Decimal	a = 1.1_d;
	float			b = 1.1f;
	double			c = 1.1;
	bool			is_double = typeid(maths::Decimal) == typeid(double);
	bool			is_float = typeid(maths::Decimal) == typeid(float);

	maths::Vec3f	vector3_A{};
	maths::Vec3f	vector3_B{ 1._d, 2._d, 3._d };
	maths::Vec3f	vector3_C{ 1._d, 0._d, 0._d };
	maths::Vec3f	vector3_D{ 0._d, 1._d, 0._d };

	maths::Vec3f	vector3_E{ maths::Cross(vector3_C, vector3_D) };
	maths::Vec3f	vector3_F{ vector3_B / 2._d };
	maths::Vec3f	vector3_G{ vector3_B * 2._d };
	maths::Vec3f	vector3_H{ vector3_B * vector3_C };
	maths::Vec3f	vector3_I{ vector3_B * vector3_D };
	maths::Vec3f	vector3_J{ vector3_B * vector3_E };
	maths::Decimal	dot_product{ maths::Dot(vector3_C, vector3_D) };

	raytracer::Film	test_film{ 100, 100, 0.035_d };

	for (int32_t i = 0; i < 100; ++i)
	{
		for (int32_t j = 0; j < 100; ++j)
		{
			test_film.SetPixel({ i / 50._d, j / 50._d, 0._d }, { i, j });
		}
	}

	test_film.WriteToFile("YOLO.png");

	maths::Vec4f A4{ 1._d, 2._d, 3._d, 4._d };
	maths::Vec4f B4{ 1._d, 2._d, 3._d, 4._d };
	maths::Vec3f A3{ 1._d, 2._d, 3._d };
	maths::Vec3f B3{ 1._d, 2._d, 3._d };
	maths::Vec4f C4{ A4 + B4 };
	maths::Vec3f C3{ A3 + B3 };

	A4 += B4;
	A3 += B3;

	maths::Mat4x4f A4x4{};
	A4x4[0][3] = 3;
	A4x4[1][3] = 2;
	maths::Vec4f D4{ 1._d, 0._d, 0._d, 1._d };
	maths::Vec4f E4{ A4x4 * D4 };

	maths::Matrix<maths::Decimal, 4, 2> A4x2{};
	maths::Matrix<maths::Decimal, 2, 4> A2x4{};
	maths::Matrix<maths::Decimal, 3, 3> A3x3(0._d);

	maths::Matrix<maths::Decimal, 2, 2> A2x2{ A2x4 * A4x2 };
	maths::Matrix<maths::Decimal, 4, 4> B4x4{ A4x2 * A2x4 };

	A3x3[0][0] = 3;
	A3x3[0][2] = 2;
	A3x3[1][0] = 2;
	A3x3[1][2] = -2;
	A3x3[2][1] = 1;
	A3x3[2][2] = 1;

	maths::Matrix<maths::Decimal, 3, 3> A3x3_inv{ maths::Inverse(A3x3) };
	maths::Matrix<maths::Decimal, 3, 3> check{ A3x3 * A3x3_inv };

	A4x2[0][1] = 4._d;
	A4x2[2][0] = 7._d;

	maths::Mat4x4f C4x4{ A3x3 };
	maths::Mat3x3f B3x3{ C4x4 };
	maths::Matrix<maths::Decimal, 3, 5> A3x5{ A4x2 };
	maths::Matrix<maths::Decimal, 4, 2> B4x2{ A3x5 };

	maths::Decimal Adec = maths::Abs(-0.2_d);
	float Af = maths::Abs<float>(-0.2f);
	double Ad = maths::Abs<double>(-0.2);
	int Ai = maths::Abs<int>(1);

	float value = maths::zero<float>;

	maths::Quaternion				AQ{ {1._d, 2._d, 3._d}, 5._d };
	maths::Normalized(AQ);

	C4x4 = static_cast<maths::Mat4x4f>(B3x3);
	maths::Transform transformA{ maths::Mat4x4f{B3x3} };

	C4x4.SetColumn(2, maths::Vec4f{ 1._d, 2._d, 3._d, 4._d });

	maths::LookAt({ 0._d, 0._d, 0._d }, { 0._d, 0._d, -1._d }, { 0._d, 1._d, 0._d });

	maths::Point3f			Point3fA{ 1._d, 2._d, 3._d };
	maths::Vector<int, 3>	Vector3iA{ Point3fA };
	maths::Norm3f			Normal3fA{ 1._d, 2._d, 3._d };
	maths::Vec3f			Vector3fA{ Normal3fA };

	maths::Bounds3f		Bounds3fA{{-1._d, -1._d, -1._d}, {1._d, 1._d, 1._d}};
	maths::Point3f		Corner0{ Bounds3fA.Corner(0) };
	maths::Point3f		Corner4{ Bounds3fA.Corner(4) };
	maths::Point3f		Corner7{ Bounds3fA.Corner(7) };

	maths::Bounds4f					Bounds4fA{ {1._d, -1._d, 1._d, -1._d}, {-1._d, 1._d, -1._d, 1._d} };
	maths::Point<maths::Decimal, 4>	Corner4f0{ Bounds4fA.Corner(0) };
	maths::Point<maths::Decimal, 4>	Corner4f7{ Bounds4fA.Corner(7) };
	maths::Point<maths::Decimal, 4>	Corner4f15{ Bounds4fA.Corner(15) };

	//maths::Decimal yolo = 0._d;
	//maths::Decimal next = maths::NextDecimalUp(yolo);
	//maths::Decimal nextnext = maths::NextDecimalUp(yolo, 2);
	//maths::Decimal delta = next - yolo;
	//maths::Decimal nextdelta = nextnext - next;

	globals::profiler_aggregate.GrabTimers(globals::profiler);
	{
		tools::Profiler::TimerTable_t timers = globals::profiler_aggregate.timers();
		for (auto it = timers.begin(); it != timers.end(); ++it)
		{
			tools::Timer const &timer = it->second;
			char const *call_string = (timer.call_count() > 1) ? " calls. " : " call. ";
			std::stringstream string;
			string <<
				timer.name() << " : " << std::endl <<
				timer.call_count() << call_string << std::endl <<
				"real time : [ " <<
				timer.total_time() << " : " << timer.worst_time() << " + " <<
				timer.average_time() << " - " <<
				timer.best_time() << " ] " << std::endl <<
				"cycles :    [ " <<
				timer.total_cycles() << " : " << timer.worst_cycles() << " + " <<
				timer.average_cycles() << " - " <<
				timer.best_cycles() << " ]";
			globals::logger.Log(tools::kChannelProfiling, tools::kLevelInfo, string.str());
		}
	}

	globals::logger.FlushAll();

	system("pause");
	return 0;
}