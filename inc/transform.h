#ifndef __YS_TRANSFORM_HPP__
#define __YS_TRANSFORM_HPP__

#include "matrix.h"
#include "raytracer.h"
#include "maths.h"
#include "ray.h"

namespace maths
{


class Transform final
{
public:
	enum OpDirection { kForward = 0, kInverse };

	explicit Transform() :
		m_{}, mInv_{}
	{}
	Transform(Mat4x4f _m) :
		m_{ _m }, mInv_{ Inverse(_m) }
	{}
	Transform(Mat4x4f _m, Mat4x4f _mInv) :
		m_{ _m }, mInv_{ _mInv }
	{}

	bool IsIdentity() const;
	bool SwapsHandedness() const;

	template <typename T> Point<T, 3> operator()(Point<T, 3> const &_v, OpDirection _dir = kForward) const;
	template <typename T> Vector<T, 3> operator()(Vector<T, 3> const &_v, OpDirection _dir = kForward) const;
	template <typename T> Vector<T, 4> operator()(Vector<T, 4> const &_v, OpDirection _dir = kForward) const;
	template <typename T> Normal<T, 3> operator()(Normal<T, 3> const &_v, OpDirection _dir = kForward) const;
	template <typename T> Bounds<T, 3> operator()(Bounds<T, 3> const &_v, OpDirection _dir = kForward) const;
	inline Ray operator()(Ray const &_v, OpDirection _dir = kForward) const;
	raytracer::SurfaceInteraction operator()(raytracer::SurfaceInteraction const &_v, OpDirection _dir = kForward) const;

	Mat4x4f const &m() const { return m_; }
	Mat4x4f const &mInv() const { return mInv_; }

private:
	Mat4x4f		m_, mInv_;
};


inline bool operator==(Transform const &_lhs, Transform const &_rhs);
inline bool operator!=(Transform const &_lhs, Transform const &_rhs);

inline Transform operator*(Transform const &_lhs, Transform const &_rhs);


inline Transform Inverse(Transform const &_v);
inline Transform Transpose(Transform const &_v);

inline Transform Translate(Vec3f const &_v);
inline Transform Scale(Decimal _x, Decimal _y, Decimal _z);

Transform RotateX(Decimal _theta);
Transform RotateY(Decimal _theta);
Transform RotateZ(Decimal _theta);
Transform Rotate(Decimal _theta, Vec3f const &_axis);

Transform LookAt(Vec3f const &_position, Vec3f const &_target, Vec3f const &_up);


class AnimatedTransform final
{
public:
	template <typename T> constexpr Point<T, 3> operator()(Decimal _t, Point<T, 3> const &_v) const;
	template <typename T> constexpr Vector<T, 3> operator()(Decimal _t, Vector<T, 3> const &_v) const;
	template <typename T> constexpr Vector<T, 4> operator()(Decimal _t, Vector<T, 4> const &_v) const;
	template <typename T> constexpr Normal<T, 3> operator()(Decimal _t, Normal<T, 3> const &_v) const;
	template <typename T> constexpr Bounds<T, 3> operator()(Decimal _t, Bounds<T, 3> const &_v) const;
	inline Ray operator()(Ray const &_v) const;

private:
	struct DerivativeTerm
	{
	};
};
} // namespace maths



namespace maths
{

template <typename T>
Point<T, 3>
Transform::operator()(Point<T, 3> const &_v, OpDirection _dir) const
{
	Mat4x4f	const &m = (_dir == kForward) ? m_ : mInv_;

	T xp = m[0][0] * _v.x + m[0][1] * _v.y + m[0][2] * _v.z + m[0][3];
	T yp = m[1][0] * _v.x + m[1][1] * _v.y + m[1][2] * _v.z + m[1][3];
	T zp = m[2][0] * _v.x + m[2][1] * _v.y + m[2][2] * _v.z + m[2][3];
	T wp = m[3][0] * _v.x + m[3][1] * _v.y + m[3][2] * _v.z + m[3][3];
	if (wp == one<T>)
		return Point<T, 3>{xp, yp, zp};
	else
		return Point<T, 3>{xp, yp, zp} / wp;
}
template <typename T>
Vector<T, 3>
Transform::operator()(Vector<T, 3> const &_v, OpDirection _dir) const
{
	Mat4x4f	const &m = (_dir == kForward) ? m_ : mInv_;

	return Vector<T, 3>{
		m[0][0] * _v.x + m[0][1] * _v.y + m[0][2] * _v.z,
		m[1][0] * _v.x + m[1][1] * _v.y + m[1][2] * _v.z,
		m[2][0] * _v.x + m[2][1] * _v.y + m[2][2] * _v.z
	};
}
template <typename T>
Vector<T, 4>
Transform::operator()(Vector<T, 4> const &_v, OpDirection _dir) const
{
	Mat4x4f	const &m = (_dir == kForward) ? m_ : mInv_;

	Vector<T, 4> result{ m * _v };
	if (result.w != one<T> && result.w != zero<T>)
		result /= result.w;
	return result;
}
template <typename T>
Normal<T, 3>
Transform::operator()(Normal<T, 3> const &_v, OpDirection _dir) const
{
	Mat4x4f	const &m = (_dir == kForward) ? m_ : mInv_;

	return Normal<T, 3>{
		m[0][0] * _v.x + m[1][0] * _v.y + m[2][0] * _v.z,
		m[0][1] * _v.x + m[1][1] * _v.y + m[2][1] * _v.z,
		m[0][2] * _v.x + m[1][2] * _v.y + m[2][2] * _v.z
	};
}
template <typename T>
Bounds<T, 3>
Transform::operator()(Bounds<T, 3> const &_v, OpDirection _dir) const
{
	Transform const &M{ *this };
	Bounds<T, 3> result{ M(_v.min, _dir) };
	result = Union(result, M(Point<T,3>{ _v.min.x, _v.min.y, _v.max.z }, _dir));
	result = Union(result, M(Point<T, 3>{ _v.min.x, _v.max.y, _v.min.z }, _dir));
	result = Union(result, M(Point<T, 3>{ _v.min.x, _v.max.y, _v.max.z }, _dir));
	result = Union(result, M(Point<T, 3>{ _v.max.x, _v.min.y, _v.min.z }, _dir));
	result = Union(result, M(Point<T, 3>{ _v.max.x, _v.min.y, _v.max.z }, _dir));
	result = Union(result, M(Point<T, 3>{ _v.max.x, _v.max.y, _v.min.z }, _dir));
	result = Union(result, M(_v.max, _dir));
	return result;
}
inline Ray
Transform::operator()(Ray const &_v, OpDirection _dir) const
{
	return Ray{ (*this)(_v.origin, _dir), (*this)(_v.direction, _dir), _v.tMax, _v.time };
}



inline bool
operator==(Transform const &_lhs, Transform const &_rhs)
{
	return _lhs.m() == _rhs.m() && _lhs.mInv() == _rhs.mInv();
}
inline bool
operator!=(Transform const &_lhs, Transform const &_rhs)
{
	return _lhs.m() != _rhs.m() || _lhs.mInv() != _rhs.mInv();
}

inline Transform
operator*(Transform const &_lhs, Transform const &_rhs)
{
	return Transform(_lhs.m() * _rhs.m(), _rhs.mInv() * _lhs.mInv());
}



inline Transform
Inverse(Transform const &_v)
{
	return Transform{ _v.mInv(), _v.m() };
}
inline Transform
Transpose(Transform const &_v)
{
	return Transform{ Transpose(_v.m()), Transpose(_v.mInv()) };
}

inline Transform
Translate(Vec3f const &_v)
{
	return Transform{
		{	1._d, 0._d, 0._d, _v.x,
			0._d, 1._d, 0._d, _v.y,
			0._d, 0._d, 1._d, _v.z,
			0._d, 0._d, 0._d, 1._d },
		{	1._d, 0._d, 0._d, -_v.x,
			0._d, 1._d, 0._d, -_v.y,
			0._d, 0._d, 1._d, -_v.z,
			0._d, 0._d, 0._d, 1._d }
	};
}
inline Transform
Scale(Decimal _x, Decimal _y, Decimal _z)
{
	return Transform{
		{	_x, 0._d, 0._d, 0._d,
			0._d, _y, 0._d, 0._d,
			0._d, 0._d, _z, 0._d,
			0._d, 0._d, 0._d, 1._d },
		{	1._d / _x, 0._d, 0._d, 0._d,
			0._d, 1._d / _y, 0._d, 0._d,
			0._d, 0._d, 1._d / _z, 0._d,
			0._d, 0._d, 0._d, 1._d }
	};
}

} // namespace maths


#endif // __YS_TRANSFORM_HPP__
