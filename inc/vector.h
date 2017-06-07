#ifndef __YS_VECTOR_HPP__
#define __YS_VECTOR_HPP__

#include <array>
#include <vector>
#include "maths.h"
#include "algorithms.h"

#ifdef YS_DEBUG
#include <set>
#endif


namespace maths
{

template <typename T, uint32_t n> struct Normal;


/// This class represents column vectors.
// NOTE: Nice post about vector libraries http://reedbeta.com/blog/on-vector-math-libraries/

template <typename T, uint32_t n>
struct Vector final
{
	Vector() :
		Vector(zero<T>)
	{}
	constexpr explicit Vector(T _value) :
		e{ std::move(algo::fill<n>::apply(_value)) }
	{}
	explicit Vector(Normal<T, n> const &_v) {
		std::copy(_v.e.begin(), _v.e.end(), e.begin());
	}
	Vector(std::initializer_list<T> _args) {
		// So close to being a constant expression. Can't wait for
		// constexpr reference std::array::operator[]
		//size_t args_size = _args.end() - _args.begin();
		//for (size_t i = 0; i < args_size; ++i)
		//	e[i] = *(_args.begin() + args_size);
		std::copy(_args.begin(), _args.end(), e.begin());
	}
	Vector(Vector<T, n - 1> const &_v, T _value)
	{ std::copy(_v.e.begin(), _v.e.end(), e.begin()); e[n - 1] = _value; }
	
	template<typename U>
	explicit Vector(Vector<U, n> const &_v)
	{
		for (uint32_t i = 0; i < size; ++i)
			e[i] = (T)_v[i];
	}

	using						value_type = T;
	static constexpr uint32_t	size{ n };

	std::array<T, n>	e;

	bool HasNaNs() const {
		for (uint32_t i = 0; i < n; ++i)
			if (std::isnan(e[i])) return true;
		return false;
	}

	T operator[](uint32_t _i) const { return e[_i]; };
	T& operator[](uint32_t _i) { return e[_i]; };
};

template <typename T>
struct Vector<T, 0> final
{
	Vector(){}
	explicit Vector(T _value){}
	Vector(std::initializer_list<T> _args){}

	using						value_type = T;
	static constexpr uint32_t	size{ 0 };

	std::array<T, 0>	e;

	bool HasNaNs() const { return false; }

	T operator[](uint32_t _i) const { return zero<T>; };
	T& operator[](uint32_t _i) { return zero<T>; };
};

template <typename T>
struct Vector<T, 4> final
{
	Vector() :
		Vector(zero<T>)
	{}
	constexpr explicit Vector(T _value) :
		e{ std::move(algo::fill<4>::apply(_value)) }
	{}
	explicit Vector(Normal<T, 4> const &_v) :
		Vector{ _v[0], _v[1], _v[2], _v[3] }
	{}
	constexpr Vector(T _e0, T _e1, T _e2, T _e3) :
		x{ _e0 }, y{ _e1 }, z{ _e2 }, w{ _e3 }
	{}
	Vector(Vector<T, 3> const &_v, T _value) :
		Vector{ _v[0], _v[1], _v[2], _value }
	{}
	template<typename U>
	explicit Vector(Vector<U, 4> const &_v)
	{
		for (uint32_t i = 0; i < size; ++i)
			e[i] = (T)_v[i];
	}

	using						value_type = T;
	static constexpr uint32_t	size{ 4 };

	union
	{
		std::array<T, 4>	e;
		struct { T	x, y, z, w; };
		struct { T	r, g, b, a; };
	};

	bool HasNaNs() const { return std::isnan(x) || std::isnan(y) || std::isnan(z) || std::isnan(w); }

	T operator[](uint32_t _i) const { return e[_i]; };
	T& operator[](uint32_t _i) { return e[_i]; };
};

template <typename T>
struct Vector<T, 3> final
{
	Vector() :
		Vector(zero<T>)
	{}
	constexpr explicit Vector(T _value) :
		e{ std::move(algo::fill<3>::apply(_value)) }
	{}
	explicit Vector(Normal<T, 3> const &_v) :
		Vector{_v.x, _v.y, _v.z}
	{}
	constexpr Vector(T _e0, T _e1, T _e2) :
		x{ _e0 }, y{ _e1 }, z{ _e2 }
	{}
	Vector(Vector<T, 2> const &_v, T _value) :
		Vector { _v.e[0], _v.e[1], _value }
	{}
	template<typename U>
	explicit Vector(Vector<U, 3> const &_v)
	{
		for (uint32_t i = 0; i < size; ++i)
			e[i] = (T)_v[i];
	}

	using						value_type = T;
	static constexpr uint32_t	size{ 3 };

	union
	{
		std::array<T, 3>	e;
		struct { T	x, y, z; };
		struct { T	r, g, b; };
	};

	bool HasNaNs() const { return std::isnan(x) || std::isnan(y) || std::isnan(z); }

	T operator[](uint32_t _i) const { return e[_i]; };
	T& operator[](uint32_t _i) { return e[_i]; };
};

template <typename T>
struct Vector<T, 2> final
{
	Vector() :
		Vector(zero<T>)
	{}
	constexpr explicit Vector(T _value) :
		e{ std::move(algo::fill<2>::apply(_value)) }
	{}
	explicit Vector(Normal<T, 2> const &_v) :
		Vector{ _v[0], _v[1] }
	{}
	constexpr Vector(T _e0, T _e1) :
		x{ _e0 }, y{ _e1 }
	{}
	Vector(Vector<T, 1> const &_v, T _value) :
		Vector{ _v[0], value }
	{}
	template<typename U>
	explicit Vector(Vector<U, 2> const &_v)
	{
		for (uint32_t i = 0; i < size; ++i)
			e[i] = (T)_v[i];
	}

	using						value_type = T;
	static constexpr uint32_t	size{ 2 };

	union
	{
		std::array<T, 2>	e;
		struct { T	x, y; };
		struct { T	u, v; };
		struct { T	w, h; };
	};

	bool HasNaNs() const { return std::isnan(x) || std::isnan(y); }

	T operator[](uint32_t _i) const { return e[_i]; };
	T& operator[](uint32_t _i) { return e[_i]; };
};


// ============================================================
// Vector<typename T, int n> operations
// ============================================================


template <typename T, uint32_t n>
bool operator==(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
bool operator!=(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);

template <typename T, uint32_t n>
Vector<T, n> &operator+=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
Vector<T, n> &operator-=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
Vector<T, n> &operator*=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
Vector<T, n> &operator/=(Vector<T, n> &_lhs, Vector<T, n> const &_rhs);

template <typename T, uint32_t n>
Vector<T, n> &operator+=(Vector<T, n> &_lhs, T _rhs);
template <typename T, uint32_t n>
Vector<T, n> &operator-=(Vector<T, n> &_lhs, T _rhs);
template <typename T, uint32_t n>
Vector<T, n> &operator*=(Vector<T, n> &_lhs, T _rhs);
template <typename T, uint32_t n>
Vector<T, n> &operator/=(Vector<T, n> &_lhs, T _rhs);

template <typename T, uint32_t n>
Vector<T, n> operator+(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
Vector<T, n> operator-(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
Vector<T, n> operator*(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
Vector<T, n> operator/(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);

template <typename T, uint32_t n>
Vector<T, n> operator+(Vector<T, n> const &_lhs, T _rhs);
template <typename T, uint32_t n>
Vector<T, n> operator-(Vector<T, n> const &_lhs, T _rhs);
template <typename T, uint32_t n>
Vector<T, n> operator*(Vector<T, n> const &_lhs, T _rhs);
template <typename T, uint32_t n>
Vector<T, n> operator/(Vector<T, n> const &_lhs, T _rhs);
template <typename T, uint32_t n>
Vector<T, n> operator*(T _lhs, Vector<T, n> const &_rhs);

template <typename T, uint32_t n>
Vector<T, n> &operator+(Vector<T, n> const &_op);
template <typename T, uint32_t n>
Vector<T, n> operator-(Vector<T, n> const &_op);


template <typename T, uint32_t n>
bool HasNaNs(Vector<T, n> const &_v);

template <typename T, uint32_t n>
Vector<T, n> Min(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
Vector<T, n> Max(Vector<T, n> const &_lhs, Vector<T, n> const &_rhs);
template <typename T, uint32_t n>
Vector<T, n> Clamp(Vector<T, n> const &_v, T _min, T _max);
template <typename T, uint32_t n>
Vector<T, n> Abs(Vector<T, n> const &_v);

template <typename T, uint32_t n>
T SqrLength(Vector<T, n> const &_v);
template <typename T, uint32_t n>
T Length(Vector<T, n> const &_v);
template <typename T, uint32_t n>
Vector<T, n> Normalized(Vector<T, n> const &_v);

template <typename T, uint32_t n>
T FoldProduct(Vector<T, n> const &_v);
template <typename T, uint32_t n>
T FoldSum(Vector<T, n> const &_v);

template <typename T, uint32_t n>
uint32_t MinimumDimension(Vector<T, n> const &_v);
template <typename T, uint32_t n>
uint32_t MaximumDimension(Vector<T, n> const &_v);

template <typename T, uint32_t n, typename... Indices>
Vector<T, n> Swizzle(Vector<T, n> const &_v, Indices... _indices);


// ============================================================
// Vector<typename T, 3> operations
// ============================================================

template <typename T>
constexpr T Dot(Vector<T, 3> const &_lhs, Vector<T, 3> const & _rhs);
template <typename T>
Vector<T, 3> Cross(Vector<T, 3> const &_lhs, Vector<T, 3> const &_rhs);
template <typename T>
Vector<T, 3> Reflect(Vector<T, 3> const &_v, Vector<T, 3> const &_n);

template <typename T>
void OrthonormalBasis(Vector<T, 3> const &_v0, Vector<T, 3> &_v1, Vector<T, 3> &_v2);


template <typename T, uint32_t n>
struct Normal final
{
	Normal() :
		Normal(zero<T>)
	{}
	explicit Normal(T _value) :
		e{ std::move(algo::fill<n>::apply(_value)) }
	{}
	explicit Normal(Vector<T, n> const &_v) {
		std::copy(_v.e.begin(), _v.e.end(), e.begin());
	}
	Normal(std::initializer_list<T> _args) {
		std::copy(_args.begin(), _args.end(), e.begin());
	}

	using						value_type = T;
	static constexpr uint32_t	size{ n };

	std::array<T, n>	e;

	bool HasNaNs() const {
		for (uint32_t i = 0; i < n; ++i)
			if (std::isnan(e[i])) return true;
		return false;
	}

	T operator[](uint32_t _i) const { return e[_i]; };
	T& operator[](uint32_t _i) { return e[_i]; };
};

template <typename T>
struct Normal<T, 3> final
{
	Normal() :
		Normal(zero<T>)
	{}
	explicit Normal(T _value) :
		e{ std::move(algo::fill<3>::apply(_value)) }
	{}
	explicit Normal(Vector<T, 3> const &_v) :
		Normal(_v[0], _v[1], _v[2])
	{}
	constexpr Normal(T _e0, T _e1, T _e2) :
		x{ _e0 }, y{ _e1 }, z{ _e2 }
	{}

	using						value_type = T;
	static constexpr uint32_t	size{ 3 };

	union
	{
		std::array<T, 3>	e;
		struct { T x, y, z; };
	};

	bool HasNaNs() const {
		for (uint32_t i = 0; i < 3; ++i)
			if (std::isnan(e[i])) return true;
		return false;
	}

	T operator[](uint32_t _i) const { return e[_i]; };
	T& operator[](uint32_t _i) { return e[_i]; };
};


// ============================================================
// Normal<typename T, uint32_t n> operations
// ============================================================

template <typename T, uint32_t n>
bool operator==(Normal<T, n> const &_lhs, Normal<T, n> const &_rhs);
template <typename T, uint32_t n>
bool operator!=(Normal<T, n> const &_lhs, Normal<T, n> const &_rhs);

template <typename T, uint32_t n>
Normal<T, n> &operator+(Normal<T, n> const &_op);
template <typename T, uint32_t n>
Normal<T, n> operator-(Normal<T, n> const &_op);

template <typename T, uint32_t n>
Normal<T, n> operator+(Normal<T, n> const &_lhs, Normal<T, n> const &_rhs);
template <typename T, uint32_t n>
Normal<T, n> &operator+=(Normal<T, n> &_lhs, Normal<T, n> const &_rhs);
template <typename T, uint32_t n>
Normal<T, n> operator-(Normal<T, n> const &_lhs, Normal<T, n> const &_rhs);
template <typename T, uint32_t n>
Normal<T, n> &operator-=(Normal<T, n> &_lhs, Normal<T, n> const &_rhs);

template <typename T, uint32_t n>
Normal<T, n> operator*(Normal<T, n> const &_lhs, T _rhs);
template <typename T, uint32_t n>
Normal<T, n> operator*(T _lhs, Normal<T, n> const &_rhs);
template <typename T, uint32_t n>
Normal<T, n> operator/(Normal<T, n> const &_lhs, T _rhs);

template <typename T, uint32_t n>
Normal<T, n> &operator*=(Normal<T, n> &_lhs, T _rhs);
template <typename T, uint32_t n>
Normal<T, n> &operator/=(Normal<T, n> &_lhs, T _rhs);


template <typename T, uint32_t n>
bool HasNaNs(Normal<T, n> const &_v);

template <typename T, uint32_t n>
T SqrLength(Normal<T, n> const &_v);
template <typename T, uint32_t n>
T Length(Normal<T, n> const &_v);
template <typename T, uint32_t n>
Normal<T, n> Normalized(Normal<T, n> const &_v);

template <typename T, uint32_t n>
Normal<T, n> FaceForward(Normal<T, n> const &_value, Normal<T, n> const &_direction);
template <typename T, uint32_t n>
Normal<T, n> FaceForward(Normal<T, n> const &_value, Vector<T, n> const &_direction);

template <typename T, uint32_t n, typename... Indices>
Normal<T, n> Swizzle(Normal<T, n> const &_v, Indices... _indices);

template <typename T>
constexpr T Dot(Normal<T, 3> const &_lhs, Normal<T, 3> const &_rhs);
template <typename T>
constexpr T Dot(Normal<T, 3> const &_lhs, Vector<T, 3> const &_rhs);
template <typename T>
constexpr T Dot(Vector<T, 3> const &_lhs, Normal<T, 3> const &_rhs);

template <typename T>
Vector<T, 3> Cross(Normal<T, 3> const &_lhs, Vector<T, 3> const &_rhs);
template <typename T>
Vector<T, 3> Cross(Vector<T, 3> const &_lhs, Normal<T, 3> const &_rhs);

template <typename T> using Vector2 = Vector<T, 2>;
template <typename T> using Vector3 = Vector<T, 3>;
template <typename T> using Vector4 = Vector<T, 4>;

using Vec3f = Vector3<Decimal>;
using Vec4f = Vector4<Decimal>;
using Vec2i32 = Vector2<int32_t>;

template <typename T> using Normal3 = Normal<T, 3>;
using Norm3f = Normal3<Decimal>;

template <typename T, uint32_t n> struct Zero<Vector<T, n>> { static constexpr Vector<T, n> value = Vector<T, n>(zero<T>); };
template <typename T, uint32_t n> struct Zero<Normal<T, n>> { static constexpr Normal<T, n> value = Normal<T, n>(zero<T>); };

} // maths


#include "vector.inl"

#endif // __YS_VECTOR_HPP__
