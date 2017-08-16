#ifndef __YS_POINT_INL__
#define __YS_POINT_INL__

#include "globals.h"
#include "logger.h"


namespace maths
{

template <typename T, uint32_t n>
Point<T, n>
operator+(Point<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] + _rhs[i];
	return result;
}
template <typename T, uint32_t n>
Point<T, n>
&operator+=(Point<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] += _rhs[i];
	return _lhs;
}
template <typename T, uint32_t n>
Point<T, n>
operator-(Point<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] - _rhs[i];
	return result;
}
template <typename T, uint32_t n>
Point<T, n>
&operator-=(Point<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] -= _rhs[i];
	return _lhs;
}

template <typename T, uint32_t n>
Vector<T, n>
operator-(Point<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] - _rhs[i];
	return result;
}

template <typename T, uint32_t n>
Point<T, n>
operator*(Point<T, n> const &_lhs, T _rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs;
	return result;
}
template <typename T, uint32_t n>
Point<T, n>
operator*(T _lhs, Point<T, n> const &_rhs)
{
	return _rhs * _lhs;
}
template <typename T, uint32_t n>
Point<T, n>
&operator*=(Point<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] *= _rhs;
	return _lhs;
}
template <typename T, uint32_t n>
Point<T, n>
operator/(Point<T, n> const &_lhs, T _rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] / _rhs;
	return result;
}
template <typename T, uint32_t n>
Point<T, n>
&operator/=(Point<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] /= _rhs;
	return _lhs;
}

template <typename T, uint32_t n>
bool
HasNaNs(Point<T, n> const &_v)
{
	for (uint32_t i = 0; i < n; ++i)
		if (std::isnan(_v[i])) return true;
	return false;
}

template <typename T, uint32_t n>
Point<T, n>
Min(Point<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Min(_lhs[i], _rhs[i]);
	return result;
}
template <typename T, uint32_t n>
Point<T, n>
Max(Point<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Max(_lhs[i], _rhs[i]);
	return result;
}
template <typename T, uint32_t n>
Point<T, n>
Clamp(Point<T, n> const &_v, T _min, T _max)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Clamp(_min[i], _max[i]);
	return result;
}
template <typename T, uint32_t n>
Point<T, n>
Blend<Point<T, n>>::Do(std::initializer_list<std::pair<Point<T, n> const &, Decimal>> _args)
{
	Decimal		summed_weights = 0._d;
	Point<T, n>	blended_point(0._d);
	for (auto &&pair : _args)
	{
		summed_weights += pair.second;
		for (uint32_t i = 0; i < n; ++i)
			blended_point[i] += pair.first[i] * pair.second;
	}
	
	YS_ASSERT(summed_weights == 1._d);
	if (summed_weights != 1._d)
		LOG_WARNING(tools::kChannelGeneral,
					"Point::Blend has been passed a list of weight with a total of " + std::to_string(summed_weights));

	return blended_point;
}

template <typename T, uint32_t n>
Point<T, n>
Floor(Point<T, n> const &_v)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = std::floor(_v[i]);
	return result;
}
template <typename T, uint32_t n>
Point<T, n>
Ceil(Point<T, n> const &_v)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = std::ceil(_v[i]);
	return result;
}
template <typename T, uint32_t n>
Point<T, n>
Abs(Point<T, n> const &_v)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = Abs(_v[i]);
	return result;
}

template <typename T, uint32_t n>
Decimal
Distance(Point<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	return Length(_lhs - _rhs);
}
template <typename T, uint32_t n>
Decimal
SqrDistance(Point<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	return SqrLength(_lhs - _rhs);
}

template <typename T, uint32_t n, typename... Indices>
Point<T, n>
Swizzle(Point<T, n> const &_v, Indices... _indices)
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
	Point<T, n> result;
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _v[swizzle[i]];
	return result;
}

} // namespace maths

#endif // __YS_POINT_INL__
