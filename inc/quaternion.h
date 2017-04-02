#ifndef __YS_QUATERNION_HPP__
#define __YS_QUATERNION_HPP__

#include "matrix.h"


namespace maths
{


template <typename T>
struct Quaternion
{
	constexpr Quaternion() :
		v(zero<T>), w{ one<T> }
	{}
	constexpr Quaternion(Vector3<T> _v, T _w) :
		v(_v), w(_w)
	{}

	Vector3<T>		v;
	T				w;
};

using Quatf = Quaternion<float>;


template <typename T>
constexpr bool operator==(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs);
template <typename T>
constexpr bool operator!=(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs);

template <typename T>
constexpr Quaternion<T> &operator+=(Quaternion<T> &_lhs, Quaternion<T> const &_rhs);
template <typename T>
constexpr Quaternion<T> &operator-=(Quaternion<T> &_lhs, Quaternion<T> const &_rhs);
template <typename T>
constexpr Quaternion<T> &operator*=(Quaternion<T> &_lhs, Quaternion<T> const &_rhs);

template <typename T>
constexpr Quaternion<T> operator+(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs);
template <typename T>
constexpr Quaternion<T> operator-(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs);
template <typename T>
constexpr Quaternion<T> operator*(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs);

template <typename T>
constexpr Quaternion<T> &operator*=(Quaternion<T> &_lhs, T _rhs);
template <typename T>
constexpr Quaternion<T> &operator/=(Quaternion<T> &_lhs, T _rhs);

template <typename T>
constexpr Quaternion<T> operator*(Quaternion<T> const &_lhs, T _rhs);
template <typename T>
constexpr Quaternion<T> operator*(T _rhs, Quaternion<T> const &_lhs);
template <typename T>
constexpr Quaternion<T> operator/(Quaternion<T> const &_lhs, T _rhs);


namespace quaternion
{

template <typename T>
constexpr T Dot(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs);
template <typename T>
constexpr Quaternion<T> Normalized(Quaternion<T> const &_v);
template <typename T>
Quaternion<T> Slerp(Quaternion<T> const &_a, Quaternion<T> const &_b, float _t);

} // namespace quaternion
} // namespace maths


#include "quaternion.inl"

#endif // __YS_QUATERNION_HPP__
