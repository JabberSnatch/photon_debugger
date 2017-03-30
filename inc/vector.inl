#ifndef __YS_VECTOR_INL__
#define __YS_VECTOR_INL__

#include "maths.h"


namespace maths
{

// ============================================================
// Vec3f
// ============================================================

constexpr Vec3f& 
Vec3f::operator+=(Vec3f const &_rhs)
{
	x += _rhs.x; y += _rhs.y; z += _rhs.z;
	return *this;
}

constexpr Vec3f&
Vec3f::operator-=(Vec3f const &_rhs)
{
	x -= _rhs.x; y -= _rhs.y; z -= _rhs.z;
	return *this;
}

constexpr Vec3f&
Vec3f::operator*=(Vec3f const &_rhs)
{
	x *= _rhs.x; y *= _rhs.y; z *= _rhs.z;
	return *this;
}

constexpr Vec3f&
Vec3f::operator/=(Vec3f const &_rhs)
{
	x /= _rhs.x; y /= _rhs.y; z /= _rhs.z;
	return *this;
}

constexpr Vec3f&
Vec3f::operator*=(float _rhs)
{
	x *= _rhs; y *= _rhs; z *= _rhs;
	return *this;
}

constexpr Vec3f&
Vec3f::operator/=(float _rhs)
{
	x /= _rhs; y /= _rhs; z /= _rhs;
	return *this;
}

constexpr Vec3f
Vec3f::normalized() const
{
	return *this / length();
}

constexpr float
Vec3f::sqr_length() const
{
	return x * x + y * y + z * z;
}

inline float
Vec3f::length() const
{
	return sqrt(sqr_length());
}


constexpr Vec3f
Clamp(Vec3f const &_v, float _min, float _max)
{
	return Vec3f{ clampf(_v.x, _min, _max),
				  clampf(_v.y, _min, _max),
				  clampf(_v.z, _min, _max) };
}

constexpr float
Dot(Vec3f const &_lhs, Vec3f const &_rhs)
{
	return _lhs.x * _rhs.x + _lhs.y * _rhs.y + _lhs.z * _rhs.z;
}

constexpr Vec3f
Cross(Vec3f const &_lhs, Vec3f const &_rhs)
{
	return	Vec3f{ _lhs.y * _rhs.z - _lhs.z * _rhs.y, 
				   _lhs.z * _rhs.x - _lhs.x * _rhs.z, 
				   _lhs.x * _rhs.y - _lhs.y * _rhs.x };
}

constexpr Vec3f
Reflect(Vec3f const &_v, Vec3f const &_n)
{
	return _v - 2.f * Dot(_v, _n) * _n;
}

constexpr Vec3f
operator+(Vec3f const &_lhs, Vec3f const &_rhs)
{
	return Vec3f{ _lhs.x + _rhs.x, _lhs.y + _rhs.y, _lhs.z + _rhs.z };
}

constexpr Vec3f
operator-(Vec3f const &_lhs, Vec3f const &_rhs)
{
	return Vec3f{ _lhs.x - _rhs.x, _lhs.y - _rhs.y, _lhs.z - _rhs.z };
}

constexpr Vec3f
operator*(Vec3f const &_lhs, Vec3f const &_rhs)
{
	return Vec3f{ _lhs.x * _rhs.x, _lhs.y * _rhs.y, _lhs.z * _rhs.z };
}

constexpr Vec3f
operator/(Vec3f const &_lhs, Vec3f const &_rhs)
{
	return Vec3f{ _lhs.x / _rhs.x, _lhs.y / _rhs.y, _lhs.z / _rhs.z };
}

constexpr Vec3f
operator*(float _s, Vec3f const &_v)
{
	return Vec3f{ _v.x * _s, _v.y * _s, _v.z * _s };
}

constexpr Vec3f
operator*(Vec3f const &_v, float _s)
{
	return _s * _v;
}

constexpr Vec3f
operator/(Vec3f const &_v, float _s)
{
	return Vec3f{ _v.x / _s, _v.y / _s, _v.z / _s };
}


// ============================================================
// Vec2i32
// ============================================================

constexpr Vec2i32&
Vec2i32::operator+=(Vec2i32 const &_rhs)
{
	x += _rhs.x; y += _rhs.y;
	return *this;
}

constexpr Vec2i32&
Vec2i32::operator-=(Vec2i32 const &_rhs)
{
	x -= _rhs.x; y -= _rhs.y;
	return *this;
}

constexpr Vec2i32&
Vec2i32::operator*=(Vec2i32 const &_rhs)
{
	x *= _rhs.x; y *= _rhs.y;
	return *this;
}

constexpr Vec2i32&
Vec2i32::operator/=(Vec2i32 const &_rhs)
{
	x /= _rhs.x; y /= _rhs.y;
	return *this;
}

constexpr Vec2i32&
Vec2i32::operator*=(int32_t _rhs)
{
	x *= _rhs; y *= _rhs;
	return *this;
}

constexpr Vec2i32&
Vec2i32::operator/=(int32_t _rhs)
{
	x /= _rhs; y /= _rhs;
	return *this;
}


constexpr Vec2i32
operator+(Vec2i32 const &_lhs, Vec2i32 const &_rhs)
{
	return Vec2i32{ _lhs.x + _rhs.x, _lhs.y + _rhs.y };
}

constexpr Vec2i32
operator-(Vec2i32 const &_lhs, Vec2i32 const &_rhs)
{
	return Vec2i32{ _lhs.x - _rhs.x, _lhs.y - _rhs.y };
}

constexpr Vec2i32
operator*(Vec2i32 const &_lhs, Vec2i32 const &_rhs)
{
	return Vec2i32{ _lhs.x * _rhs.x, _lhs.y * _rhs.y };
}

constexpr Vec2i32
operator/(Vec2i32 const &_lhs, Vec2i32 const &_rhs)
{
	return Vec2i32{ _lhs.x / _rhs.x, _lhs.y / _rhs.y };
}

constexpr Vec2i32
operator*(int32_t _s, Vec2i32 const &_v)
{
	return Vec2i32{ _v.x * _s, _v.y * _s };
}

constexpr Vec2i32
operator*(Vec2i32 const &_v, int32_t _s)
{
	return _s * _v;
}

constexpr Vec2i32
operator/(Vec2i32 const &_v, int32_t _s)
{
	return Vec2i32{ _v.x / _s, _v.y / _s };
}

} // namespace maths


#endif // __YS_VECTOR_INL__
