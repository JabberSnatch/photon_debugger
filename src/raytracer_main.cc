#include "profiler.h"
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

//static tools::Profiler gProfiler{};

#include <typeinfo>
#include <iostream>

int main()
{
	raytracer::Film		film_35mm{ 480, 270, 0.035_d };
	raytracer::Camera	camera{
		film_35mm,
		{0._d, -5._d, .5_d}, {0._d, 1._d, .5_d}, {0._d, 0._d, 1._d},
		60._d
	};
	maths::Transform	little_sphere_transform = maths::Translate({ 0._d, 1._d, .5_d });
	raytracer::Sphere	little_sphere{ little_sphere_transform, false, .5_d, -5._d, 5._d, 360._d };
	maths::Transform	big_sphere_transform = maths::Translate({ 0._d, 1._d, -500._d });
	raytracer::Sphere	big_sphere{ big_sphere_transform, false, 500._d, -500._d, 500._d, 360._d };

	raytracer::Shape const &sphereA = little_sphere;
	raytracer::Shape const &sphereB = big_sphere;
	raytracer::Camera::Scene	scene = { std::ref(sphereA), std::ref(sphereB) };

	camera.Expose(scene, 0._d);
	camera.film.WriteToFile("big_little_sphere.png");


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


	tools::Profiler::TimerTable_t timers = globals::profiler.timers();
	for (auto it = timers.begin(); it != timers.end(); ++it)
	{
		tools::Timer const &timer = it->second;
		char const *call_string = (timer.call_count() > 1) ? " calls. " : " call. ";
		std::cout <<
			timer.name() << " : " <<
			timer.call_count() << call_string << std::endl <<
			timer.total_time() << " seconds, " <<
			timer.best_time() << " seconds at least, " <<
			timer.total_cycles() << " cycles, " <<
			timer.best_cycles() << " cycles as least" <<
		std::endl;
	}

	system("pause");
	return 0;
}