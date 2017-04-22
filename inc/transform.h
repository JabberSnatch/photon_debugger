#ifndef __YS_TRANSFORM_HPP__
#define __YS_TRANSFORM_HPP__

#include "matrix.h"
#include "maths.h"
#include "ray.h"

namespace maths
{


class Transform final
{
public:
	Transform() = delete;
	Transform(Mat4x4f _m) :
		m_{ _m }, mInv_{ matrix::Inverse(_m) }
	{}
	constexpr Transform(Mat4x4f _m, Mat4x4f _mInv) :
		m_{ _m }, mInv_{ _mInv }
	{}

	bool IsIdentity() const;

	template <typename T> constexpr Point<T, 3> operator()(Point<T, 3> const &_v) const;
	template <typename T> constexpr Vector<T, 3> operator()(Vector<T, 3> const &_v) const;
	template <typename T> constexpr Vector<T, 4> operator()(Vector<T, 4> const &_v) const;
	template <typename T> constexpr Normal<T, 3> operator()(Normal<T, 3> const &_v) const;
	template <typename T> constexpr Bounds<T, 3> operator()(Bounds<T, 3> const &_v) const;
	inline Ray operator()(Ray const &_v) const;

	Mat4x4f const &m() const { return m_; }
	Mat4x4f const &mInv() const { return mInv_; }

private:
	Mat4x4f		m_, mInv_;
};


inline bool operator==(Transform const &_lhs, Transform const &_rhs);
inline bool operator!=(Transform const &_lhs, Transform const &_rhs);

inline Transform operator*(Transform const &_lhs, Transform const &_rhs);


namespace transform
{

inline Transform Inverse(Transform const &_v);
inline Transform Transpose(Transform const &_v);

inline Transform Translate(Vec3f const &_v);
inline Transform Scale(Decimal _x, Decimal _y, Decimal _z);

Transform RotateX(float _theta);
Transform RotateY(float _theta);
Transform RotateZ(float _theta);
Transform Rotate(float _theta, Vec3f const &_axis);

Transform LookAt(Vec3f const &_position, Vec3f const &_target, Vec3f const &_up);

} // namespace transform


class AnimatedTransform final
{
public:
	template <typename T> constexpr Point<T, 3> operator()(float _t, Point<T, 3> const &_v) const;
	template <typename T> constexpr Vector<T, 3> operator()(float _t, Vector<T, 3> const &_v) const;
	template <typename T> constexpr Vector<T, 4> operator()(float _t, Vector<T, 4> const &_v) const;
	template <typename T> constexpr Normal<T, 3> operator()(float _t, Normal<T, 3> const &_v) const;
	template <typename T> constexpr Bounds<T, 3> operator()(float _t, Bounds<T, 3> const &_v) const;
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
constexpr Point<T, 3>
Transform::operator()(Point<T, 3> const &_v) const
{
	T xp = m_[0][0] * _v.x + m_[0][1] * _v.y + m_[0][2] * _v.z + m_[0][3];
	T yp = m_[1][0] * _v.x + m_[1][1] * _v.y + m_[1][2] * _v.z + m_[1][3];
	T zp = m_[2][0] * _v.x + m_[2][1] * _v.y + m_[2][2] * _v.z + m_[2][3];
	T wp = m_[3][0] * _v.x + m_[3][1] * _v.y + m_[3][2] * _v.z + m_[3][3];
	if (wp == one<T>)
		return Point<T, 3>{xp, yp, zp};
	else
		return Point<T, 3>{xp, yp, zp} / wp;
}
template <typename T>
constexpr Vector<T, 3>
Transform::operator()(Vector<T, 3> const &_v) const
{
	return Vector<T, 3>{
		m_[0][0] * _v.x + m_[0][1] * _v.y + m_[0][2] * _v.z,
		m_[1][0] * _v.x + m_[1][1] * _v.y + m_[1][2] * _v.z,
		m_[2][0] * _v.x + m_[2][1] * _v.y + m_[2][2] * _v.z
	};
}
template <typename T>
constexpr Vector<T, 4>
Transform::operator()(Vector<T, 4> const &_v) const
{
	Vector<T, 4> result{ m_ * _v };
	if (result.w != one<T> && result.w != zero<T>)
		result /= result.w;
	return result;
}
template <typename T>
constexpr Normal<T, 3>
Transform::operator()(Normal<T, 3> const &_v) const
{
	return Normal<T, 3>{
		m_[0][0] * _v.x + m_[1][0] * _v.y + m_[2][0] * _v.z,
		m_[0][1] * _v.x + m_[1][1] * _v.y + m_[2][1] * _v.z,
		m_[0][2] * _v.x + m_[1][2] * _v.y + m_[2][2] * _v.z
	};
}
template <typename T>
constexpr Bounds<T, 3>
Transform::operator()(Bounds<T, 3> const &_v) const
{
	Transform const &M{ *this };
	Bounds<T, 3> result{ M(_v.min) };
	result = bounds::Union(result, M({ _v.min.x, _v.min.y, _v.max.z }));
	result = bounds::Union(result, M({ _v.min.x, _v.max.y, _v.min.z }));
	result = bounds::Union(result, M({ _v.min.x, _v.max.y, _v.max.z }));
	result = bounds::Union(result, M({ _v.max.x, _v.min.y, _v.min.z }));
	result = bounds::Union(result, M({ _v.max.x, _v.min.y, _v.max.z }));
	result = bounds::Union(result, M({ _v.max.x, _v.max.y, _v.min.z }));
	result = bounds::Union(result, M(_v.max));
	return result;
}
inline Ray
Transform::operator()(Ray const &_v) const
{
	return Ray{ (*this)(_v.origin), (*this)(_v.direction), _v.tMax, _v.time };
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


namespace transform
{

inline Transform
Inverse(Transform const &_v)
{
	return Transform{ _v.mInv(), _v.m() };
}
inline Transform
Transpose(Transform const &_v)
{
	return Transform{ matrix::Transpose(_v.m()), matrix::Transpose(_v.mInv()) };
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

} // namespace transform
} // namespace maths


#endif // __YS_TRANSFORM_HPP__
