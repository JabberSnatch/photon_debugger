#ifndef __YS_QUATERNION_HPP__
#define __YS_QUATERNION_HPP__

#include "maths.h"
#include "matrix.h"
#include "transform.h"


namespace maths
{


struct Quaternion
{
	Quaternion() :
		v(0._d), w{ 0._d }
	{}
	Quaternion(Vector3<Decimal> _v, Decimal _w) :
		v(_v), w(_w)
	{}

	explicit Quaternion(Transform const &_v);
	explicit operator Transform() const;

	Vector3<Decimal>	v;
	Decimal				w;
};


bool operator==(Quaternion const &_lhs, Quaternion const &_rhs);
bool operator!=(Quaternion const &_lhs, Quaternion const &_rhs);

Quaternion &operator+=(Quaternion &_lhs, Quaternion const &_rhs);
Quaternion &operator-=(Quaternion &_lhs, Quaternion const &_rhs);
Quaternion &operator*=(Quaternion &_lhs, Quaternion const &_rhs);

Quaternion operator+(Quaternion const &_lhs, Quaternion const &_rhs);
Quaternion operator-(Quaternion const &_lhs, Quaternion const &_rhs);
Quaternion operator*(Quaternion const &_lhs, Quaternion const &_rhs);

Quaternion &operator*=(Quaternion &_lhs, Decimal _rhs);
Quaternion &operator/=(Quaternion &_lhs, Decimal _rhs);

Quaternion operator*(Quaternion const &_lhs, Decimal _rhs);
Quaternion operator*(Decimal _rhs, Quaternion const &_lhs);
Quaternion operator/(Quaternion const &_lhs, Decimal _rhs);


namespace quaternion
{

Quaternion Conjugate(Quaternion const &_v);
Quaternion Inverse(Quaternion const &_v);
Decimal Norm(Quaternion const &_v);
Decimal Dot(Quaternion const &_lhs, Quaternion const &_rhs);
Quaternion Normalized(Quaternion const &_v);
Quaternion Slerp(Quaternion const &_a, Quaternion const &_b, Decimal _t);

} // namespace quaternion
} // namespace maths


#endif // __YS_QUATERNION_HPP__
