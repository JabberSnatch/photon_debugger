#ifndef __YS_VECTOR_INL__
#define __YS_VECTOR_INL__

#include "maths.h"


namespace maths
{


// ============================================================
// Vector<typename T, int n> operations
// ============================================================

template <typename T, uint32_t n>
constexpr bool
operator==(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		if (_lhs[i] != _rhs[i]) return false;
	return true;
}
template <typename T, uint32_t n>
constexpr bool
operator!=(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		if (_lhs[i] != _rhs[i]) return true;
	return false;
}

template <typename T, uint32_t n> 
constexpr Vector<T, n>
&operator+=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] += _rhs[i];
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator-=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] -= _rhs[i];
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator*=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] *= _rhs[i];
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator/=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] /= _rhs[i];
	return _lhs;
}

template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator+=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] += _rhs;
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator-=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] -= _rhs;
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator*=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] *= _rhs;
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator/=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] /= _rhs;
	return _lhs;
}

template <typename T, uint32_t n>
constexpr Vector<T, n>
operator+(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] + _rhs[i];
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator-(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] - _rhs[i];
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator*(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs[i];
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator/(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] / _rhs[i];
	return result;
}

template <typename T, uint32_t n>
constexpr Vector<T, n>
operator+(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] + _rhs;
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator-(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] - _rhs;
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator*(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs;
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator/(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] / _rhs;
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator*(T _lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs;
	return result;
}

template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator+(Vector<T, n> const &_op)
{
	return _op;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator-(Vector<T, n> const &_op)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = -_op[i];
	return result;
}


namespace vector
{

template <typename T, uint32_t n>
constexpr Vector<T, n>
Clamp(Vector<T, n> const &_v, T _min, T _max)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Clamp<T>(_v[i], _min, _max);
	return result;
}

template <typename T, uint32_t n>
constexpr T
SqrLength(Vector<T, n> const &_v)
{
	T result{ zero<T> };
	for (uint32_t i = 0; i < n; ++i)
		result += _v[i] * _v[i];
	return result;
}

template <typename T, uint32_t n>
constexpr T
FoldProduct(Vector<T, n> const &_v)
{
	T result{ zero<T> };
	for (uint32_t i = 0; i < n; ++i)
		result *= _v;
	return result;
}
template <typename T, uint32_t n>
constexpr T
FoldSum(Vector<T, n> const &_v)
{
	T result{ zero<T> };
	for (uint32_t i = 0; i < n; ++i)
		result += _v;
	return result;
}


// ============================================================
// Vector<typename T, 3> operations
// ============================================================

template <typename T>
constexpr T
Dot(Vector<T, 3> const &_lhs, Vector<T, 3> const & _rhs)
{
	return _lhs.x * _rhs.x + _lhs.y * _rhs.y + _lhs.z * _rhs.z;
}
template <typename T>
constexpr Vector<T, 3>
Cross(Vector<T, 3> const &_lhs, Vector<T, 3> const &_rhs)
{
	return Vector<T, 3>{ _lhs.y * _rhs.z - _lhs.z * _rhs.y,
						 _lhs.z * _rhs.x - _lhs.x * _rhs.z,
						 _lhs.x * _rhs.y - _lhs.y * _rhs.x };
}
template <typename T>
constexpr Vector<T, 3>
Reflect(Vector<T, 3> const &_v, Vector<T, 3> const &_n)
{
	return _v + 2 * Dot(_v, _n) * _n;
}


// ============================================================
// Vector<float, uint32_t n> operations
// ============================================================

template <uint32_t n>
inline float
Length(Vector<float, n> const &_v)
{
	return std::sqrtf(SqrLength(_v));
}
template <uint32_t n>
constexpr Vector<float, n>
Normalized(Vector<float, n> const &_v)
{
	return _v / Length(_v);
}

} // namespace vector
} // namespace maths


#endif // __YS_VECTOR_INL__
