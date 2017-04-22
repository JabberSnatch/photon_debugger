#ifndef __YS_BOUNDS_HPP__
#define __YS_BOUNDS_HPP__

#include "point.h"

namespace maths
{


template <typename T, uint32_t n>
struct Bounds final
{
	constexpr Bounds() :
		min(LowestValue<T>), max(HighestValue<T>)
	{}
	explicit constexpr Bounds(Point<T, n> const &_v) :
		min{ _v }, max{ _v }
	{}
	constexpr Bounds(Point<T, n> const &_a, Point<T, n> const &_b) :
		min{ point::Min(_a, _b) }, max{ point::Max(_a, _b) }
	{}

	constexpr Point<T, n> const &operator[](uint32_t _i) const;
	constexpr Point<T, n> &operator[] (uint32_t _i);

	constexpr Point<T, n> Corner(uint32_t _index) const;
	constexpr Vector<T, n> Diagonal() const;
	constexpr T SurfaceArea() const;
	constexpr T Volume() const;
	constexpr uint32_t MaximumExtent() const;

	constexpr Point<T, n> Lerp(Point<float, n> const &_t) const;
	constexpr Vector<T, n> Offset(Point<T, n> const &_p) const;

	constexpr void BoundingSphere(Point<T, n> &_center, T &_radius) const;

	Point<T, n> min, max;
};


template <typename T, uint32_t n>
constexpr bool operator==(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr bool operator!=(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs);


namespace bounds
{

template <typename T, uint32_t n>
constexpr Bounds<T, n> Union(Bounds<T, n> const &_lhs, Point<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Bounds<T, n> Union(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Bounds<T, n> Intersect(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr bool Overlap(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr bool Inside(Bounds<T, n> const &_lhs, Point<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr bool InsideExclusive(Bounds<T, n> const &_lhs, Point<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Bounds<T, n> Expand(Bounds<T, n> const &_v, T _delta);

} // namespace bounds


using Bounds2f = Bounds<maths::Decimal, 2>;
using Bounds3f = Bounds<maths::Decimal, 3>;
using Bounds4f = Bounds<maths::Decimal, 4>;


} // namespace maths

#include "bounds.inl"


#endif // __YS_BOUNDS_HPP__
