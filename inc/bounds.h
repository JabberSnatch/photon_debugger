#ifndef __YS_BOUNDS_HPP__
#define __YS_BOUNDS_HPP__

#include "point.h"

namespace maths
{


template <typename T, uint32_t n>
struct Bounds final
{
	Bounds() :
		min(highest_value<T>), max(-highest_value<T>)
	{}
	explicit constexpr Bounds(Point<T, n> const &_v) :
		min{ _v }, max{ _v }
	{}
	Bounds(Point<T, n> const &_a, Point<T, n> const &_b) :
		min{ Min(_a, _b) }, max{ Max(_a, _b) }
	{}

	// NOTE: These would be contexpr without assertion.
	Point<T, n> const &operator[](uint32_t _i) const;
	Point<T, n> &operator[] (uint32_t _i);

	Point<T, n>		Corner(uint32_t _index) const;
	Vector<T, n>	Diagonal() const;
	T				SurfaceArea() const;
	T				Volume() const;
	uint32_t		MaximumExtent() const;

	Point<T, n>		Lerp(Point<float, n> const &_t) const;
	Vector<T, n>	Offset(Point<T, n> const &_p) const;

	void			BoundingSphere(Point<T, n> &_center, T &_radius) const;

	Point<T, n> min, max;
};


template <typename T, uint32_t n>
bool operator==(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs);
template <typename T, uint32_t n>
bool operator!=(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs);


template <typename T, uint32_t n>
Bounds<T, n> Union(Bounds<T, n> const &_lhs, Point<T, n> const &_rhs);
template <typename T, uint32_t n>
Bounds<T, n> Union(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs);
template <typename T, uint32_t n>
Bounds<T, n> Intersect(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs);
template <typename T, uint32_t n>
bool Overlap(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs);
template <typename T, uint32_t n>
bool Inside(Bounds<T, n> const &_lhs, Point<T, n> const &_rhs);
template <typename T, uint32_t n>
bool InsideExclusive(Bounds<T, n> const &_lhs, Point<T, n> const &_rhs);
template <typename T, uint32_t n>
Bounds<T, n> Expand(Bounds<T, n> const &_v, T _delta);


using Bounds2f = Bounds<maths::Decimal, 2>;
using Bounds3f = Bounds<maths::Decimal, 3>;
using Bounds4f = Bounds<maths::Decimal, 4>;


} // namespace maths

#include "bounds.inl"


#endif // __YS_BOUNDS_HPP__
