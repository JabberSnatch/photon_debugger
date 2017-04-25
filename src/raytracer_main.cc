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

//static tools::Profiler gProfiler{};

#include <typeinfo>

int main()
{
	maths::Decimal	a = 1.1_d;
	float			b = 1.1f;
	double			c = 1.1;
	bool			is_double = typeid(maths::Decimal) == typeid(double);
	bool			is_float = typeid(maths::Decimal) == typeid(float);

	for (int j = 0; j < 10; ++j)
	{
		tools::TimeProbe	probe{ globals::profiler.GetTimer("TEST") };
		for (int i = 0; i < 1000; ++i)
		{
		}
	}

	maths::Vec3f	vector3_A{};
	maths::Vec3f	vector3_B{ 1._d, 2._d, 3._d };
	maths::Vec3f	vector3_C{ 1._d, 0._d, 0._d };
	maths::Vec3f	vector3_D{ 0._d, 1._d, 0._d };

	maths::Vec3f	vector3_E{ maths::vector::Cross(vector3_C, vector3_D) };
	maths::Vec3f	vector3_F{ vector3_B / 2._d };
	maths::Vec3f	vector3_G{ vector3_B * 2._d };
	maths::Vec3f	vector3_H{ vector3_B * vector3_C };
	maths::Vec3f	vector3_I{ vector3_B * vector3_D };
	maths::Vec3f	vector3_J{ vector3_B * vector3_E };
	maths::Decimal	dot_product{ maths::vector::Dot(vector3_C, vector3_D) };

	raytracer::Film	test_film{ 100, 100 };

	for (int32_t i = 0; i < 100; ++i)
	{
		for (int32_t j = 0; j < 100; ++j)
		{
			test_film.SetPixel({ i / 50._d, j / 50._d, 0._d }, { i, j });
		}
	}

	test_film.WriteToFile("YOLO.png");

	auto timers{ globals::profiler.timers() };

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

	maths::Matrix<maths::Decimal, 3, 3> A3x3_inv{ maths::matrix::Inverse(A3x3) };
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

	maths::Quaternion<float>	AQ{};
	maths::quaternion::Normalized(AQ);

	C4x4 = static_cast<maths::Mat4x4f>(B3x3);
	maths::Transform transformA{ maths::Mat4x4f{B3x3} };

	C4x4.SetColumn(2, maths::Vec4f{ 1._d, 2._d, 3._d, 4._d });

	maths::transform::LookAt({ 0._d, 0._d, 0._d }, { 0._d, 0._d, -1._d }, { 0._d, 1._d, 0._d });

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

	maths::Decimal yolo = 0._d;
	maths::Decimal next = maths::NextDecimalUp(yolo);
	maths::Decimal nextnext = maths::NextDecimalUp(yolo, 2);
	maths::Decimal delta = next - yolo;
	maths::Decimal nextdelta = nextnext - next;

	int decimal_shift;
#ifdef YS_DECIMAL_IS_DOUBLE
	decimal_shift = 53;
#else
	decimal_shift = 24;
#endif
	/*
	maths::Decimal nextexp = maths::BitsToDecimal((maths::DecimalBits(1) << decimal_shift) - 1);
	maths::DecimalBits nextexpbits = maths::DecimalToBits(nextexp);
	maths::Decimal nextexpnext = maths::NextDecimalUp(nextexp);
	maths::DecimalBits nextexpnextbits = maths::DecimalToBits(nextexpnext);
	maths::Decimal nextexpdelta = nextexpnext - nextexp;
	maths::Decimal nextexpdeltanext = maths::NextDecimalUp(nextexpnext) - nextexpnext;
	
	maths::DecimalBits fullexpbits = (((~((maths::DecimalBits(1) << decimal_shift) - 1)) << 1) >> 1);
	maths::Decimal fullexp = maths::BitsToDecimal(fullexpbits);
	maths::Decimal fullexpnext = maths::NextDecimalUp(fullexp);
	maths::Decimal fullexppenult = maths::NextDecimalDown(fullexp);
	maths::DecimalBits fullexpnextbits = maths::DecimalToBits(fullexpnext);
	maths::Decimal fullexpdelta = fullexpnext - fullexp;

	maths::Decimal fullexpnextepsilon = fullexpnext * (1._d + maths::machine_epsilon);
	maths::Decimal fullexpnextmepsilon = fullexpnext * (1._d - maths::machine_epsilon);
	maths::Decimal fullexpnextepsilondelta = (fullexpnextepsilon - fullexpnextmepsilon);

	maths::DecimalBits epsilonbits = maths::DecimalToBits(fullexpnextepsilon);
	maths::DecimalBits mepsilonbits = maths::DecimalToBits(fullexpnextmepsilon);
	maths::DecimalBits epsilondeltabits = maths::DecimalToBits(fullexpnextepsilondelta);

	auto foo = maths::GetDecimalBits(3._d).bits;
	auto baz = maths::GetDecimalBits(1_db).value;
	*/

	system("pause");
	return 0;
}