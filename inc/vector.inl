#ifndef __YS_VECTOR_INL__
#define __YS_VECTOR_INL__

#include "maths.h"
#include "ys_assert.h"


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
	YS_ASSERT(!_lhs.HasNaNs());
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator-=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] -= _rhs[i];
	YS_ASSERT(!_lhs.HasNaNs());
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator*=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] *= _rhs[i];
	YS_ASSERT(!_lhs.HasNaNs());
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator/=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] /= _rhs[i];
	YS_ASSERT(!_lhs.HasNaNs());
	return _lhs;
}

template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator+=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] += _rhs;
	YS_ASSERT(!_lhs.HasNaNs());
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator-=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] -= _rhs;
	YS_ASSERT(!_lhs.HasNaNs());
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator*=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] *= _rhs;
	YS_ASSERT(!_lhs.HasNaNs());
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
&operator/=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] /= _rhs;
	YS_ASSERT(!_lhs.HasNaNs());
	return _lhs;
}

template <typename T, uint32_t n>
constexpr Vector<T, n>
operator+(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] + _rhs[i];
	YS_ASSERT(!result.HasNaNs());
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator-(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] - _rhs[i];
	YS_ASSERT(!result.HasNaNs());
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator*(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs[i];
	YS_ASSERT(!result.HasNaNs());
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator/(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] / _rhs[i];
	YS_ASSERT(!result.HasNaNs());
	return result;
}

template <typename T, uint32_t n>
constexpr Vector<T, n>
operator+(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] + _rhs;
	YS_ASSERT(!result.HasNaNs());
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator-(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] - _rhs;
	YS_ASSERT(!result.HasNaNs());
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator*(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs;
	YS_ASSERT(!result.HasNaNs());
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator/(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] / _rhs;
	YS_ASSERT(!result.HasNaNs());
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
operator*(T _lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs;
	YS_ASSERT(!result.HasNaNs());
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
constexpr bool
HasNaNs(Vector<T, n> const &_v)
{
	for (uint32_t i = 0; i < n; ++i)
		if (std::isnan(_v[i])) return true;
	return false;
}

template <typename T, uint32_t n>
constexpr Vector<T, n>
Min(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Min(_lhs[i], _rhs[i]);
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
Max(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Max(_lhs[i], _rhs[i]);
	return result;
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
Clamp(Vector<T, n> const &_v, T _min, T _max)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Clamp<T>(_v[i], _min, _max);
	YS_ASSERT(!result.HasNaNs());
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
inline T
Length(Vector<T, n> const &_v)
{
	return std::sqrt(SqrLength(_v));
}
template <typename T, uint32_t n>
constexpr Vector<T, n>
Normalized(Vector<T, n> const &_v)
{
	Vector<T, n> result{ _v / Length(_v) };
	YS_ASSERT(!result.HasNaNs());
	return result;
}

template <typename T, uint32_t n>
constexpr T
FoldProduct(Vector<T, n> const &_v)
{
	T result{ one<T> };
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

template <typename T, uint32_t n>
constexpr uint32_t
MinimumDimension(Vector<T, n> const &_v)
{
	uint32_t	result = 0u;
	T			min = _v[0];

	for (uint32_t i = 1; i < n; ++i)
	{
		if (_v[i] < min)
		{
			result = i;
			min = _v[i];
		}
	}

	return result;
}
template <typename T, uint32_t n>
constexpr uint32_t
MaximumDimension(Vector<T, n> const &_v)
{
	uint32_t	result = 0u;
	T			max = _v[0];
	
	for (uint32_t i = 1; i < n; ++i)
	{
		if (_v[i] > max)
		{
			result = i;
			max = _v[i];
		}
	}

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
	Vector<T, 3> result{ _lhs.y * _rhs.z - _lhs.z * _rhs.y,
						 _lhs.z * _rhs.x - _lhs.x * _rhs.z,
						 _lhs.x * _rhs.y - _lhs.y * _rhs.x };
	YS_ASSERT(!result.HasNaNs());
	return result;
}
template <typename T>
constexpr Vector<T, 3>
Reflect(Vector<T, 3> const &_v, Vector<T, 3> const &_n)
{
	Vector<T, 3> result = _v + 2 * Dot(_v, _n) * _n;
	YS_ASSERT(!result.HasNaNs());
	return result;
}


} // namespace vector


// ============================================================
// Normal<typename T, uint32_t n> operations
// ============================================================

template <typename T, uint32_t n>
constexpr bool
operator==(Normal<T, n> const &_lhs, Normal<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		if (_lhs[i] != _rhs[i]) return false;
	return true;
}
template <typename T, uint32_t n>
constexpr bool
operator!=(Normal<T, n> const &_lhs, Normal<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		if (_lhs[i] != _rhs[i]) return true;
	return false;
}

template <typename T, uint32_t n>
constexpr Normal<T, n>
&operator+(Normal<T, n> const &_op)
{
	return _op;
}
template <typename T, uint32_t n>
constexpr Normal<T, n>
operator-(Normal<T, n> const &_op)
{
	Normal<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = -_op[i];
	YS_ASSERT(!result.HasNaNs());
	return result;
}

template <typename T, uint32_t n>
constexpr Normal<T, n>
operator*(Normal<T, n> const &_lhs, T _rhs)
{
	Normal<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs;
	YS_ASSERT(!result.HasNaNs());
	return result;
}
template <typename T, uint32_t n>
constexpr Normal<T, n>
operator*(T _lhs, Normal<T, n> const &_rhs)
{
	return _rhs * _lhs;
}
template <typename T, uint32_t n>
constexpr Normal<T, n>
operator/(Normal<T, n> const &_lhs, T _rhs)
{
	Normal<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] / _rhs;
	YS_ASSERT(!result.HasNaNs());
	return result;
}

template <typename T, uint32_t n>
constexpr Normal<T, n>
&operator*=(Normal<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] *= _rhs;
	YS_ASSERT(!_lhs.HasNaNs());
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Normal<T, n>
&operator/=(Normal<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] /= _rhs;
	YS_ASSERT(!_lhs.HasNaNs());
	return _lhs;
}


namespace normal
{

template <typename T, uint32_t n>
constexpr bool
HasNaNs(Normal<T, n> const &_v)
{
	for (uint32_t i = 0; i < n; ++i)
		if (std::isnan(_v[i])) return true;
	return false;
}

template <typename T, uint32_t n>
constexpr T
SqrLength(Normal<T, n> const &_v)
{
	T result{ zero<T> };
	for (uint32_t i = 0; i < n; ++i)
		result += _v[i] * _v[i];
	return result;
}
template <typename T, uint32_t n>
inline T
Length(Normal<T, n> const &_v)
{
	return std::sqrt(SqrLength(_v));
}
template <typename T, uint32_t n>
constexpr Normal<T, n>
Normalized(Normal<T, n> const &_v)
{
	Normal<T, n> result{ _v / Length(_v) };
	YS_ASSERT(!result.HasNaNs());
	return result;
}

template <typename T, uint32_t n>
constexpr T
Dot(Normal<T, n> const &_lhs, Normal<T, n> const &_rhs)
{
	T result{ zero<T> };
	for (uint32_t i = 0; i < n; ++i)
		result += _lhs[i] * _rhs[i];
	return result;
}
template <typename T, uint32_t n>
constexpr T
Dot(Normal<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	T result{ zero<T> };
	for (uint32_t i = 0; i < n; ++i)
		result += _lhs[i] * _rhs[i];
	return result;
}
template <typename T, uint32_t n>
constexpr T
Dot(Vector<T, n> const &_lhs, Normal<T, n> const &_rhs)
{
	T result{ zero<T> };
	for (uint32_t i = 0; i < n; ++i)
		result += _lhs[i] * _rhs[i];
	return result;
}

template <typename T, uint32_t n>
constexpr Normal<T, n>
FaceForward(Normal<T, n> const &_value, Normal<T, n> const &_direction)
{
	return (Dot(_value, _direction) > zero<T>) ? _value : -_value;
}
template <typename T, uint32_t n>
constexpr Normal<T, n>
FaceForward(Normal<T, n> const &_value, Vector<T, n> const &_direction)
{
	return (Dot(_value, _direction) > zero<T>) ? _value : -_value;
}

} // namespace normal
} // namespace maths


#endif // __YS_VECTOR_INL__
