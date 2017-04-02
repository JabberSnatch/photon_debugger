#ifndef __YS_TRANSFORM_HPP__
#define __YS_TRANSFORM_HPP__

#include "matrix.h"

namespace maths
{


class Transform
{
public:
	Transform() = delete;
	Transform(Mat4x4f _m) :
		m_{ _m }, mInv_{ matrix::Inverse(_m) }
	{}
	Transform(Mat4x4f _m, Mat4x4f _mInv) :
		m_{ _m }, mInv_{ _mInv }
	{}

	bool IsIdentity() const;

	Mat4x4f const &m() const { return m_; }
	Mat4x4f const &mInv() const { return mInv_; }

private:
	Mat4x4f		m_, mInv_;
};


inline bool operator==(Transform const &_lhs, Transform const &_rhs)
{ return _lhs.m() == _rhs.m() && _lhs.mInv() == _rhs.mInv(); }
inline bool operator!=(Transform const &_lhs, Transform const &_rhs)
{ return _lhs.m() != _rhs.m() || _lhs.mInv() != _rhs.mInv(); }


namespace transform
{

inline Transform Inverse(Transform const &_v) 
{ return Transform{ _v.mInv(), _v.m() }; }
inline Transform Transpose(Transform const &_v) 
{ return Transform{ matrix::Transpose(_v.m()), matrix::Transpose(_v.mInv()) }; }

inline Transform Translate(Vec3f const &_v) 
{
	return Transform{
		{1.f, 0.f, 0.f, _v.x, 0.f, 1.f, 0.f, _v.y, 0.f, 0.f, 1.f, _v.z, 0.f, 0.f, 0.f, 1.f},
		{1.f, 0.f, 0.f, -_v.x, 0.f, 1.f, 0.f, -_v.y, 0.f, 0.f, 1.f, -_v.z, 0.f, 0.f, 0.f, 1.f}
	};
}
inline Transform Scale(float _x, float _y, float _z)
{
	return Transform{
		{_x, 0.f, 0.f, 0.f, 0.f, _y, 0.f, 0.f, 0.f, 0.f, _z, 0.f, 0.f, 0.f, 0.f, 1.f},
		{1.f / _x, 0.f, 0.f, 0.f, 1.f / _y, 0.f, 0.f, 0.f, 0.f, 1.f / _z, 0.f, 0.f, 0.f, 0.f, 1.f}
	};
}

Transform RotateX(float _theta);
Transform RotateY(float _theta);
Transform RotateZ(float _theta);
Transform Rotate(float _theta, Vec3f const &_axis);

Transform LookAt(Vec3f const &_position, Vec3f const &_target, Vec3f const &_up);

} // namespace transform
} // namespace maths


#endif // __YS_TRANSFORM_HPP__
