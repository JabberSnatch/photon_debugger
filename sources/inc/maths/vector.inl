#ifndef __YS_VECTOR_INL__
#define __YS_VECTOR_INL__

#include "maths/maths.h"


namespace maths
{

// ============================================================
// Vector<typename T, int n> operations
// ============================================================


template <typename T, uint32_t n>
bool
operator==(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		if (_lhs[i] != _rhs[i]) return false;
	return true;
}
template <typename T, uint32_t n>
bool
operator!=(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		if (_lhs[i] != _rhs[i]) return true;
	return false;
}

template <typename T, uint32_t n> 
Vector<T, n>
&operator+=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] += _rhs[i];
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}
template <typename T, uint32_t n>
Vector<T, n>
&operator-=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] -= _rhs[i];
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}
template <typename T, uint32_t n>
Vector<T, n>
&operator*=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] *= _rhs[i];
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}
template <typename T, uint32_t n>
Vector<T, n>
&operator/=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] /= _rhs[i];
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}

template <typename T, uint32_t n>
Vector<T, n>
&operator+=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] += _rhs;
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}
template <typename T, uint32_t n>
Vector<T, n>
&operator-=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] -= _rhs;
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}
template <typename T, uint32_t n>
Vector<T, n>
&operator*=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] *= _rhs;
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}
template <typename T, uint32_t n>
Vector<T, n>
&operator/=(Vector<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] /= _rhs;
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}

template <typename T, uint32_t n>
Vector<T, n>
operator+(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] + _rhs[i];
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
operator-(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] - _rhs[i];
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
operator*(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs[i];
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
operator/(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] / _rhs[i];
	YS_ASSERT(!HasNaNs(result));
	return result;
}

template <typename T, uint32_t n>
Vector<T, n>
operator+(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] + _rhs;
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
operator-(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] - _rhs;
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
operator*(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs;
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
operator/(Vector<T, n> const &_lhs, T _rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] / _rhs;
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
operator*(T _lhs, Vector<T, n> const &_rhs)
{
	return _rhs * _lhs;
}

template <typename T, uint32_t n>
Vector<T, n>
&operator+(Vector<T, n> const &_op)
{
	return _op;
}
template <typename T, uint32_t n>
Vector<T, n>
operator-(Vector<T, n> const &_op)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = -_op[i];
	return result;
}



template <typename T, uint32_t n>
struct HasNaNs_Impl_
{
	static bool value(Vector<T, n> const &_v)
	{
		return false;
	}
};
template <uint32_t n>
struct HasNaNs_Impl_<float, n>
{
	static bool value(Vector<float, n> const &_v)
	{
		for (uint32_t i = 0; i < n; ++i)
			if (std::isnan(_v[i])) return true;
		return false;
	}
};
template <uint32_t n>
struct HasNaNs_Impl_<double, n>
{
	static bool value(Vector<double, n> const &_v)
	{
		for (uint32_t i = 0; i < n; ++i)
			if (std::isnan(_v[i])) return true;
		return false;
	}
};
template <typename T, uint32_t n>
bool
HasNaNs(Vector<T, n> const &_v)
{
	return HasNaNs_Impl_<T, n>::value(_v);
}

template <typename T, uint32_t n>
Vector<T, n>
Min(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Min(_lhs[i], _rhs[i]);
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
Max(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Max(_lhs[i], _rhs[i]);
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
Clamp(Vector<T, n> const &_v, T _min, T _max)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Clamp<T>(_v[i], _min, _max);
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
Abs(Vector<T, n> const &_v)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Abs<T>(_v[i]);
	YS_ASSERT(!HasNaNs(result));
	return result;
}

template <typename T, uint32_t n>
T
SqrLength(Vector<T, n> const &_v)
{
	T result{ zero<T> };
	for (uint32_t i = 0; i < n; ++i)
		result += _v[i] * _v[i];
	return result;
}
template <typename T, uint32_t n>
T
Length(Vector<T, n> const &_v)
{
	return std::sqrt(SqrLength(_v));
}
template <typename T, uint32_t n>
Vector<T, n>
Normalized(Vector<T, n> const &_v)
{
	Vector<T, n> result{ _v / Length(_v) };
	YS_ASSERT(!HasNaNs(result));
	return result;
}

template <typename T, uint32_t n>
T
FoldProduct(Vector<T, n> const &_v)
{
	T result{ one<T> };
	for (uint32_t i = 0; i < n; ++i)
		result *= _v[i];
	return result;
}
template <typename T, uint32_t n>
T
FoldSum(Vector<T, n> const &_v)
{
	T result{ zero<T> };
	for (uint32_t i = 0; i < n; ++i)
		result += _v[i];
	return result;
}

template <typename T, uint32_t n>
uint32_t
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
uint32_t
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
template <typename T, uint32_t n>
T
MaximumComponent(Vector<T, n> const &_v)
{
	T  result = zero<T>;
	for (uint32_t i = 0; i < n; ++i)
	{
		result = maths::Max(result, _v[i]);
	}
	return result;
}

template <typename T, uint32_t n, typename... Indices>
Vector<T, n>
Swizzle(Vector<T, n> const &_v, Indices... _indices)
{
	static_assert(!(sizeof...(_indices) < (size_t)n), "Not enough swizzle parameters.");
	static_assert(!(sizeof...(_indices) > (size_t)n), "Too many swizzle parameters.");

	std::vector<uint32_t> const	swizzle{ _indices... };
#ifdef YS_DEBUG
	std::set<uint32_t> const	duplicate_check{ _indices... };
	YS_ASSERT(duplicate_check.size() == n);
	for (uint32_t i = 0; i < n; ++i)
		YS_ASSERT(swizzle[i] < n);
#endif
	Vector<T, n> result;
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _v[swizzle[i]];
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
Vector<T, 3>
Cross(Vector<T, 3> const &_lhs, Vector<T, 3> const &_rhs)
{
	Vector<T, 3> result{ _lhs.y * _rhs.z - _lhs.z * _rhs.y,
						 _lhs.z * _rhs.x - _lhs.x * _rhs.z,
						 _lhs.x * _rhs.y - _lhs.y * _rhs.x };
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T>
Vector<T, 3>
Reflect(Vector<T, 3> const &_v, Vector<T, 3> const &_n)
{
	Vector<T, 3> result = 2 * Dot(_v, _n) * _n - _v;
	YS_ASSERT(!HasNaNs(result));
	return result;
}

template <typename T>
void
OrthonormalBasis(Vector<T, 3> const &_v0, Vector<T, 3> &_v1, Vector<T, 3> &_v2)
{
	if (Abs(_v0.x) > Abs(_v0.y))
		_v1 = maths::Vector<T, 3>{ -_v0.z, 0._d, _v0.x } / std::sqrt(_v0.x * _v0.x + _v0.z * _v0.z);
	else
		_v1 = maths::Vector<T, 3>{ 0._d, _v0.z, -_v0.y } / std::sqrt(_v0.y * _v0.y + _v0.z * _v0.z);
	_v2 = Cross(_v0, _v1);
}


// ============================================================
// Normal<typename T, uint32_t n> operations
// ============================================================

template <typename T, uint32_t n>
bool
operator==(Normal<T, n> const &_lhs, Normal<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		if (_lhs[i] != _rhs[i]) return false;
	return true;
}
template <typename T, uint32_t n>
bool
operator!=(Normal<T, n> const &_lhs, Normal<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		if (_lhs[i] != _rhs[i]) return true;
	return false;
}

template <typename T, uint32_t n>
Normal<T, n>
&operator+(Normal<T, n> const &_op)
{
	return _op;
}
template <typename T, uint32_t n>
Normal<T, n>
operator-(Normal<T, n> const &_op)
{
	Normal<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = -_op[i];
	YS_ASSERT(!HasNaNs(result));
	return result;
}

template <typename T, uint32_t n>
Normal<T, n>
operator+(Normal<T, n> const &_lhs, Normal<T, n> const &_rhs)
{
	Normal<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] + _rhs[i];
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T, uint32_t n>
Normal<T, n>&
operator+=(Normal<T, n> &_lhs, Normal<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] += _rhs[i];
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}
template <typename T, uint32_t n>
Normal<T, n>
operator-(Normal<T, n> const &_lhs, Normal<T, n> const &_rhs)
{
	Normal<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] - _rhs[i];
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T, uint32_t n>
Normal<T, n>&
operator-=(Normal<T, n> &_lhs, Normal<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] -= _rhs[i];
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}

template <typename T, uint32_t n>
Normal<T, n>
operator*(Normal<T, n> const &_lhs, T _rhs)
{
	Normal<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs;
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T, uint32_t n>
Normal<T, n>
operator*(T _lhs, Normal<T, n> const &_rhs)
{
	return _rhs * _lhs;
}
template <typename T, uint32_t n>
Normal<T, n>
operator/(Normal<T, n> const &_lhs, T _rhs)
{
	Normal<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] / _rhs;
	YS_ASSERT(!HasNaNs(result));
	return result;
}

template <typename T, uint32_t n>
Normal<T, n>
&operator*=(Normal<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] *= _rhs;
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}
template <typename T, uint32_t n>
Normal<T, n>
&operator/=(Normal<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] /= _rhs;
	YS_ASSERT(!HasNaNs(_lhs));
	return _lhs;
}


template <typename T, uint32_t n>
bool
HasNaNs(Normal<T, n> const &_v)
{
	for (uint32_t i = 0; i < n; ++i)
		if (std::isnan(_v[i])) return true;
	return false;
}

template <typename T, uint32_t n>
T
SqrLength(Normal<T, n> const &_v)
{
	T result{ zero<T> };
	for (uint32_t i = 0; i < n; ++i)
		result += _v[i] * _v[i];
	return result;
}
template <typename T, uint32_t n>
T
Length(Normal<T, n> const &_v)
{
	return std::sqrt(SqrLength(_v));
}
template <typename T, uint32_t n>
Normal<T, n>
Normalized(Normal<T, n> const &_v)
{
	Normal<T, n> result{ _v / Length(_v) };
	YS_ASSERT(!HasNaNs(result));
	return result;
}

template <typename T, uint32_t n>
Normal<T, n>
FaceForward(Normal<T, n> const &_value, Normal<T, n> const &_direction)
{
	return (Dot(_value, _direction) > zero<T>) ? _value : -_value;
}
template <typename T, uint32_t n>
Normal<T, n>
FaceForward(Normal<T, n> const &_value, Vector<T, n> const &_direction)
{
	return (Dot(_value, _direction) > zero<T>) ? _value : -_value;
}

template <typename T, uint32_t n, typename... Indices>
Normal<T, n>
Swizzle(Normal<T, n> const &_v, Indices... _indices)
{
	static_assert(!(sizeof...(_indices) < (size_t)n), "Not enough swizzle parameters.");
	static_assert(!(sizeof...(_indices) > (size_t)n), "Too many swizzle parameters.");

	std::vector<uint32_t> const	swizzle{ _indices... };
#ifdef YS_DEBUG
	std::set<uint32_t> const	duplicate_check{ _indices... };
	YS_ASSERT(duplicate_check.size() == n);
	for (uint32_t i = 0; i < n; ++i)
		YS_ASSERT(swizzle[i] < n);
#endif
	Normal<T, n> result;
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _v[swizzle[i]];
	return result;
}

template <typename T>
constexpr T Dot(Normal<T, 3> const &_lhs, Normal<T, 3> const &_rhs)
{
	return _lhs.x * _rhs.x + _lhs.y * _rhs.y + _lhs.z * _rhs.z;
}
template <typename T>
constexpr T Dot(Normal<T, 3> const &_lhs, Vector<T, 3> const &_rhs)
{
	return _lhs.x * _rhs.x + _lhs.y * _rhs.y + _lhs.z * _rhs.z;
}
template <typename T>
constexpr T
Dot(Vector<T, 3> const &_lhs, Normal<T, 3> const &_rhs)
{
	return _lhs.x * _rhs.x + _lhs.y * _rhs.y + _lhs.z * _rhs.z;
}

template <typename T>
Vector<T, 3>
Cross(Normal<T, 3> const &_lhs, Vector<T, 3> const &_rhs)
{
	Vector<T, 3> result{ _lhs.y * _rhs.z - _lhs.z * _rhs.y,
						 _lhs.z * _rhs.x - _lhs.x * _rhs.z, 
						 _lhs.x * _rhs.y - _lhs.y * _rhs.x };
	YS_ASSERT(!HasNaNs(result));
	return result;
}
template <typename T>
Vector<T, 3>
Cross(Vector<T, 3> const &_lhs, Normal<T, 3> const &_rhs)
{
	Vector<T, 3> result{ _lhs.y * _rhs.z - _lhs.z * _rhs.y,
						 _lhs.z * _rhs.x - _lhs.x * _rhs.z,
						 _lhs.x * _rhs.y - _lhs.y * _rhs.x };
	YS_ASSERT(!HasNaNs(result));
	return result;
}

} // namespace maths


#endif // __YS_VECTOR_INL__
