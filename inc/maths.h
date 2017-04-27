#pragma once
#ifndef __YS_MATHS_HPP__
#define __YS_MATHS_HPP__

#include <cstdint>
#include <limits>
#include <algorithm>

#define YS_DECIMAL_IS_DOUBLE

#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif


namespace maths
{

template <typename T, uint32_t n> struct Vector;
template <typename T, uint32_t n> struct Normal;
template <typename T, uint32_t n> struct Point;
template <typename T, uint32_t n> struct Bounds;
template <typename T, uint32_t R, uint32_t C> struct Matrix;

struct Ray;

class Transform;
struct Quaternion;

struct REDecimal;

} // namespace maths


namespace maths
{

template <typename T> struct Zero { static constexpr T value = T::zero; };
template <> struct Zero<float> { static constexpr float value = 0.f; };
template <> struct Zero<double> { static constexpr double value = 0.; };
template <> struct Zero<int32_t> { static constexpr int32_t value = 0; };
template <> struct Zero<int64_t> { static constexpr int64_t value = 0; };
template <> struct Zero<uint32_t> { static constexpr uint32_t value = 0u; };
template <> struct Zero<uint64_t> { static constexpr uint64_t value = 0u; };
template <typename T> constexpr T zero = Zero<T>::value;

template <typename T> struct One { static constexpr T value = T::one; };
template <> struct One<float> { static constexpr float value = 1.f; };
template <> struct One<double> { static constexpr double value = 1.; };
template <> struct One<int32_t> { static constexpr int32_t value = 1; };
template <> struct One<int64_t> { static constexpr int64_t value = 1; };
template <> struct One<uint32_t> { static constexpr uint32_t value = 1u; };
template <> struct One<uint64_t> { static constexpr uint64_t value = 1u; };
template <typename T> constexpr T one = One<T>::value;

} // namespace maths


// Floating point utilities
namespace maths
{

// NOTE: IEEE754 float is 1 sign, 8 exponent, 23 significand (32bits)
//		 (sign=0, exponent=1, significand=0 is at 1 << 24)
//		 IEEE754 double is 1 sign, 11 exponent, 52 significand (64bits)
//		 (sign=0, exponent=1, significand=0 is at 1 << 53)
//		 IEEE754 quadruple is 1 sign, 15 exponent, 113 significand (128bits) (TMYK)
template <typename T> struct FloatMeta
{
	using BitfieldType = uint64_t;
	static constexpr auto sign_mask = 0Ui64;
	static constexpr auto exponent_mask = 0Ui64;
	static constexpr auto significand_mask = 0Ui64;
};
template <> struct FloatMeta<float>
{
	using BitfieldType = uint32_t;
	static constexpr auto sign_mask = 0x80000000U;
	static constexpr auto exponent_mask = 0x7F800000U;
	static constexpr auto significand_mask = 0x7fffffU;
};
template <> struct FloatMeta<double>
{
	using BitfieldType = uint64_t;
	static constexpr auto sign_mask = 0x8000000000000000Ui64;
	static constexpr auto exponent_mask = 0x7ff0000000000000Ui64;
	static constexpr auto significand_mask = 0xfffffffffffffUi64;
};


#ifdef YS_DECIMAL_IS_DOUBLE
using Decimal = double;
using DecimalBits = uint64_t;
#else
using Decimal = float;
using DecimalBits = uint32_t;
#endif

static_assert(sizeof(Decimal) == sizeof(DecimalBits), "Decimal and DecimalBits are of different sizes.");
using DecimalMeta = FloatMeta<Decimal>;


template <typename V, typename B>
union ValueBitsMapper
{
	using ValueType = V;
	using BitsType = B;
	ValueBitsMapper(ValueType _value) : value{ _value } {}
	ValueBitsMapper(BitsType _bits) : bits{ _bits } {}
	ValueType	value;
	BitsType	bits;
};

using DecimalBitsMapper = ValueBitsMapper<Decimal, DecimalBits>;
using FloatBitsMapper = ValueBitsMapper<float, uint32_t>;
using DoubleBitsMapper = ValueBitsMapper<double, uint64_t>;

double	NextDecimalUp(double _v, uint64_t _delta = 1);
double	NextDecimalDown(double _v, uint64_t _delta = 1);
float	NextDecimalUp(float _v, uint32_t _delta = 1);
float	NextDecimalDown(float _v, uint32_t _delta = 1);

} // namespace maths
constexpr maths::Decimal operator "" _d(long double _v) { return maths::Decimal(_v); }
constexpr maths::DecimalBits operator "" _db(unsigned long long _v) { return maths::DecimalBits(_v); }
//inline maths::Decimal operator "" _d(char const *_v) { return maths::Decimal(std::atof(_v)); }


namespace maths
{

bool	Quadratic(Decimal _a, Decimal _b, Decimal _c, Decimal &_t0, Decimal &_t1);

template <typename T> static constexpr T lowest_value = std::numeric_limits<T>::lowest();
template <typename T> static constexpr T highest_value = std::numeric_limits<T>::max();
template <typename T> static constexpr T infinity = std::numeric_limits<T>::infinity();
template <typename T> static constexpr T almost_one = one<T> - std::numeric_limits<T>::epsilon();

template <typename T> static constexpr T pi = T( 3.14159235658979323846 );

template <typename T> constexpr T Radians(T _degrees) { return (pi<T> / T( 180 )) * _degrees; }
template <typename T> constexpr T Degrees(T _radians) { return (T( 180 ) / pi<T>) * _radians; }

// NOTE: Maybe a Scalar<T> class could make these functions a little more specific.
template <typename T> constexpr T Min(T _lhs, T _rhs) { return (_lhs < _rhs) ? _lhs : _rhs; }
template <typename T> constexpr T Max(T _lhs, T _rhs) { return (_lhs > _rhs) ? _lhs : _rhs; }
template <typename T> constexpr T Clamp(T _v, T _min, T _max) { return Min(Max(_v, _min), _max); }
template <typename T> constexpr T SafeClamp(T _v, T _a, T _b) { return Min(Max(_v, Min(_a, _b)), Max(_a, _b)); }

template <typename T> constexpr T Abs(T _v) { return (_v > zero<T>) ? _v : -_v; }

template <typename T> constexpr T Lerp(T _a, T _b, float _t) { return _a*(1.f - _t) + _b*_t; }
template <typename T> constexpr T Lerp(T _a, T _b, double _t) { return _a*(1. - _t) + _b*_t; }

} // namespace maths


// ============================================================================
//				 /!\	BEYOND THIS POINT LIES DARKNESS   /!\
// ============================================================================

namespace crappy_legacy {
namespace maths {

template <typename T> bool	IsInf(T _v) { return false; }
template <typename T> bool	IsNaN(T _v) { return false; }
template <> bool			IsInf<float>(float _v);
template <> bool			IsNaN<float>(float _v);
template <> bool			IsInf<double>(double _v);
template <> bool			IsNaN<double>(double _v);

float	InvSqrt(float _v);
double	InvSqrt(double _v);
float	Sqrt(float _v);
double	Sqrt(double _v);


} // namespace maths
} // namespace crappy_legacy


#endif // __YS_MATHS_HPP__
