#ifndef __YS_VECTOR_INL__
#define __YS_VECTOR_INL__


namespace maths
{

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

} // namespace maths


#endif // __YS_VECTOR_INL__
