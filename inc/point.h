#ifndef __YS_POINT_HPP__
#define __YS_POINT_HPP__

#include <array>
#include "algorithms.h"
#include "vector.h"


namespace maths
{


template <typename T, uint32_t n>
struct Point final
{
	constexpr Point() :
		Point(zero<T>)
	{}
	explicit constexpr Point(T _value) :
		e{ algo::fill<n>::apply(_value) }
	{}
	constexpr Point(std::initializer_list<T> _args) {
		std::copy(_args.begin(), _args.end(), e.begin());
	}

	std::array<T, n>	e;

	constexpr bool HasNaNs() {
		for (uint32_t i = 0; i < n; ++i)
			if (std::isnan(e[i])) return true;
		return false;
	}

	constexpr T operator[](uint32_t _i) const { return e[_i]; };
	constexpr T& operator[](uint32_t _i) { return e[_i]; };
};


template <typename T> using Point3 = Point<T, 3>;
template <typename T> using Point2 = Point<T, 2>;
using Point3f = Point3<float>;
using Point2f = Point2<float>;


// ============================================================
// Point<typename T, int n> operations
// ============================================================


template <typename T, uint32_t n>
constexpr Point<T, n> operator+(Point<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Point<T, n> &operator+=(Point<T, n> &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Point<T, n> operator-(Point<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Point<T, n> &operator-=(Point<T, n> &_lhs, Vector<T, n> const &_rhs);

template <typename T, uint32_t n>
constexpr Vector<T, n> operator-(Point<T, n> const &_lhs, Point<T, n> const &_rhs);

template <typename T, uint32_t n>
constexpr Point<T, n> operator*(Point<T, n> const &_lhs, T _rhs);
template <typename T, uint32_t n>
constexpr Point<T, n> operator*(T _lhs, Point<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Point<T, n> &operator*=(Point<T, n> &_lhs, T _rhs);
template <typename T, uint32_t n>
constexpr Point<T, n> operator/(Point<T, n> const &_lhs, T _rhs);
template <typename T, uint32_t n>
constexpr Point<T, n> &operator/=(Point<T, n> &_lhs, T _rhs);

namespace point
{

template <typename T, uint32_t n>
constexpr bool HasNaNs(Point<T, n> const &_v);

template <typename T, uint32_t n>
constexpr Point<T, n> Min(Point<T, n> const &_lhs, Point<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Point<T, n> Max(Point<T, n> const &_lhs, Point<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Point<T, n> Clamp(Point<T, n> const &_v, T _min, T _max);

template <typename T, uint32_t n>
constexpr Point<T, n> Floor(Point<T, n> const &_v);
template <typename T, uint32_t n>
constexpr Point<T, n> Ceil(Point<T, n> const &_v);
template <typename T, uint32_t n>
constexpr Point<T, n> Abs(Point<T, n> const &_v);

template <typename T, uint32_t n>
constexpr float Distance(Point<T, n> const &_lhs, Point<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr float SqrDistance(Point<T, n> const &_lhs, Point<T, n> const &_rhs);

} // namespace point
} // namespace maths

#include "point.inl"


#endif // __YS_POINT_HPP__