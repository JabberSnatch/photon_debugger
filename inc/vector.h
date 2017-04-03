#ifndef __YS_VECTOR_HPP__
#define __YS_VECTOR_HPP__

#include <array>
#include "algorithms.h"


namespace maths
{


/// This class represents column vectors.
// NOTE: Nice post about vector libraries http://reedbeta.com/blog/on-vector-math-libraries/

template <typename T, uint32_t n>
struct Vector final
{
	constexpr Vector() :
		Vector(zero<T>)
	{}
	explicit constexpr Vector(T _value) :
		e{ algo::fill<n>::apply(_value) }
	{}
	constexpr Vector(std::initializer_list<T> _args) {
		std::copy(_args.begin(), _args.end(), e.begin());
	}
	constexpr Vector(Vector<T, n - 1> const &_v, T _value)
	{ std::copy(_v.e.begin(), _v.e.end(), e.begin()); e[n - 1] = _value; }

	std::array<T, n>	e;

	constexpr bool HasNaNs() {
		for (uint32_t i = 0; i < n; ++i)
			if (std::isnan(e[i])) return true;
		return false;
	}

	constexpr T operator[](uint32_t _i) const { return e[_i]; };
	constexpr T& operator[](uint32_t _i) { return e[_i]; };
};

template <typename T>
struct Vector<T, 0> final
{
	constexpr Vector(){}
	explicit constexpr Vector(T _value){}
	constexpr Vector(std::initializer_list<T> _args){}

	std::array<T, 0>	e;

	constexpr bool HasNaNs() { return false; }

	constexpr T operator[](uint32_t _i) const { return zero<T>; };
	constexpr T& operator[](uint32_t _i) { return zero<T>; };
};

template <typename T>
struct Vector<T, 4> final
{
	constexpr Vector() :
		Vector(zero<T>)
	{}
	explicit constexpr Vector(T _value) :
		e{ algo::fill<4>::apply(_value) }
	{}
	constexpr Vector(T _e0, T _e1, T _e2, T _e3) :
		x{ _e0 }, y{ _e1 }, z{ _e2 }, w{ _e3 }
	{}
	constexpr Vector(Vector<T, 3> const &_v, T _value)
	{ std::copy(_v.e.begin(), _v.e.end(), e.begin()); e[3] = _value; }

	union
	{
		std::array<T, 4>	e;
		struct { T	x, y, z, w; };
		struct { T	r, g, b, a; };
	};

	constexpr bool HasNaNs() { return std::isnan(x) || std::isnan(y) || std::isnan(z) || std::isnan(w); }

	constexpr T operator[](uint32_t _i) const { return e[_i]; };
	constexpr T& operator[](uint32_t _i) { return e[_i]; };
};

template <typename T>
struct Vector<T, 3> final
{
	constexpr Vector() :
		Vector(zero<T>)
	{}
	explicit constexpr Vector(T _value) :
		e{ algo::fill<3>::apply(_value) }
	{}
	constexpr Vector(T _e0, T _e1, T _e2) :
		x{ _e0 }, y{ _e1 }, z{ _e2 }
	{}
	constexpr Vector(Vector<T, 2> const &_v, T _value)
	{ std::copy(_v.e.begin(), _v.e.end(), e.begin()); e[2] = _value; }

	union
	{
		std::array<T, 3>	e;
		struct { T	x, y, z; };
		struct { T	r, g, b; };
	};

	constexpr bool HasNaNs() { return std::isnan(x) || std::isnan(y) || std::isnan(z); }

	constexpr T operator[](uint32_t _i) const { return e[_i]; };
	constexpr T& operator[](uint32_t _i) { return e[_i]; };
};

template <typename T>
struct Vector<T, 2> final
{
	constexpr Vector() :
		Vector(zero<T>)
	{}
	explicit constexpr Vector(T _value) :
		e{ algo::fill<2>::apply(_value) }
	{}
	constexpr Vector(T _e0, T _e1) :
		x{ _e0 }, y{ _e1 }
	{}
	constexpr Vector(Vector<T, 1> const &_v, T _value)
	{ std::copy(_v.e.begin(), _v.e.end(), e.begin()); e[1] = _value; }

	union
	{
		std::array<T, 2>	e;
		struct { T	x, y; };
		struct { T	u, v; };
		struct { T	w, h; };
	};

	constexpr bool HasNaNs() { return std::isnan(x) || std::isnan(y); }

	constexpr T operator[](uint32_t _i) const { return e[_i]; };
	constexpr T& operator[](uint32_t _i) { return e[_i]; };
};


// ============================================================
// Vector<typename T, int n> operations
// ============================================================


template <typename T, uint32_t n>
constexpr bool operator==(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr bool operator!=(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);

template <typename T, uint32_t n>
constexpr Vector<T, n> &operator+=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> &operator-=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> &operator*=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> &operator/=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs);

template <typename T, uint32_t n>
constexpr Vector<T, n> &operator+=(Vector<T, n> &_lhs, T _rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> &operator-=(Vector<T, n> &_lhs, T _rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> &operator*=(Vector<T, n> &_lhs, T _rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> &operator/=(Vector<T, n> &_lhs, T _rhs);

template <typename T, uint32_t n>
constexpr Vector<T, n> operator+(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> operator-(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> operator*(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> operator/(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);

template <typename T, uint32_t n>
constexpr Vector<T, n> operator+(Vector<T, n> const &_lhs, T _rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> operator-(Vector<T, n> const &_lhs, T _rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> operator*(Vector<T, n> const &_lhs, T _rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> operator/(Vector<T, n> const &_lhs, T _rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> operator*(T _lhs, Vector<T, n> const &_rhs);

template <typename T, uint32_t n>
constexpr Vector<T, n> &operator+(Vector<T, n> const &_op);
template <typename T, uint32_t n>
constexpr Vector<T, n> operator-(Vector<T, n> const &_op);


namespace vector
{

template <typename T, uint32_t n>
constexpr bool HasNaNs(Vector<T, n> const &_v);

template <typename T, uint32_t n>
constexpr Vector<T, n> Min(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> Max(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
constexpr Vector<T, n> Clamp(Vector<T, n> const &_v, T _min, T _max);

template <typename T, uint32_t n>
constexpr T SqrLength(Vector<T, n> const &_v);
template <typename T, uint32_t n>
inline T Length(Vector<T, n> const &_v);
template <typename T, uint32_t n>
constexpr Vector<T, n> Normalized(Vector<T, n> const &_v);

template <typename T, uint32_t n>
constexpr T FoldProduct(Vector<T, n> const &_v);
template <typename T, uint32_t n>
constexpr T FoldSum(Vector<T, n> const &_v);


// ============================================================
// Vector<typename T, 3> operations
// ============================================================

template <typename T>
constexpr T Dot(Vector<T, 3> const &_lhs, Vector<T, 3> const & _rhs);
template <typename T>
constexpr Vector<T, 3> Cross(Vector<T, 3> const &_lhs, Vector<T, 3> const &_rhs);
template <typename T>
constexpr Vector<T, 3> Reflect(Vector<T, 3> const &_v, Vector<T, 3> const &_n);

} // namespace vector


template <typename T> using Vector2 = Vector<T, 2>;
template <typename T> using Vector3 = Vector<T, 3>;
template <typename T> using Vector4 = Vector<T, 4>;

using Vec3f = Vector3<float>;
using Vec4f = Vector4<float>;
using Vec2i32 = Vector2<int32_t>;

} // maths


#include "vector.inl"

#endif // __YS_VECTOR_HPP__
