#ifndef __YS_POINT_INL__
#define __YS_POINT_INL__


namespace maths
{

template <typename T, uint32_t n>
constexpr Point<T, n>
operator+(Point<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] + _rhs[i];
	return result;
}
template <typename T, uint32_t n>
constexpr Point<T, n>
&operator+=(Point<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] += _rhs[i];
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Point<T, n>
operator-(Point<T, n> const &_lhs, Vector<T, n> const &_rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] - _rhs[i];
	return result;
}
template <typename T, uint32_t n>
constexpr Point<T, n>
&operator-=(Point<T, n> &_lhs, Vector<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] -= _rhs[i];
	return _lhs;
}

template <typename T, uint32_t n>
constexpr Vector<T, n>
operator-(Point<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	Vector<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] - _rhs[i];
	return result;
}

template <typename T, uint32_t n>
constexpr Point<T, n>
operator*(Point<T, n> const &_lhs, T _rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] * _rhs;
	return result;
}
template <typename T, uint32_t n>
constexpr Point<T, n>
operator*(T _lhs, Point<T, n> const &_rhs)
{
	return _rhs * _lhs;
}
template <typename T, uint32_t n>
constexpr Point<T, n>
&operator*=(Point<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] *= _rhs;
	return _lhs;
}
template <typename T, uint32_t n>
constexpr Point<T, n>
operator/(Point<T, n> const &_lhs, T _rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = _lhs[i] / _rhs;
	return result;
}
template <typename T, uint32_t n>
constexpr Point<T, n>
&operator/=(Point<T, n> &_lhs, T _rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		_lhs[i] /= _rhs;
	return _lhs;
}

namespace point
{

template <typename T, uint32_t n>
constexpr bool
HasNaNs(Point<T, n> const &_v)
{
	for (uint32_t i = 0; i < n; ++i)
		if (std::isnan(_v[i])) return true;
	return false;
}

template <typename T, uint32_t n>
constexpr Point<T, n>
Min(Point<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Min(_lhs[i], _rhs[i]);
	return result;
}
template <typename T, uint32_t n>
constexpr Point<T, n>
Max(Point<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Max(_lhs[i], _rhs[i]);
	return result;
}
template <typename T, uint32_t n>
constexpr Point<T, n>
Clamp(Point<T, n> const &_v, T _min, T _max)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Clamp(_lhs[i], _rhs[i]);
	return result;
}

template <typename T, uint32_t n>
constexpr Point<T, n>
Floor(Point<T, n> const &_v)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = std::floor(_v[i]);
	return result;
}
template <typename T, uint32_t n>
constexpr Point<T, n>
Ceil(Point<T, n> const &_v)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = std::ceil(_v[i]);
	return result;
}
template <typename T, uint32_t n>
constexpr Point<T, n>
Abs(Point<T, n> const &_v)
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = Abs(_v[i]);
	return result;
}

template <typename T, uint32_t n>
constexpr Decimal
Distance(Point<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	return vector::Length(_lhs - _rhs);
}
template <typename T, uint32_t n>
constexpr Decimal
SqrDistance(Point<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	return vector::SqrLength(_lhs - _rhs);
}

} // namespace point
} // namespace maths

#endif // __YS_POINT_INL__
