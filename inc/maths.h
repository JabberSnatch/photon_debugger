#pragma once
#ifndef __YS_MATHS_HPP__
#define __YS_MATHS_HPP__


//#define YS_DECIMAL_IS_DOUBLE

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
template <typename T> struct Quaternion;

struct Ray;

class Transform;

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


} // namespace maths
constexpr maths::Decimal operator "" _d(long double _v) { return maths::Decimal(_v); }
constexpr maths::DecimalBits operator "" _db(unsigned long long _v) { return maths::DecimalBits(_v); }
//inline maths::Decimal operator "" _d(char const *_v) { return maths::Decimal(std::atof(_v)); }

namespace maths
{

static constexpr Decimal machine_epsilon = std::numeric_limits<Decimal>::epsilon() * 0.5_d;
// NOTE: Higham (2002, section 3.1) bounding term for (1+-Epsilon)^n for n*Epsilon < 1
constexpr Decimal gamma(uint32_t _n) { return (_n * machine_epsilon) / (1._d - _n * machine_epsilon); }

template <typename V, typename B>
union ValueBitsMapper
{
	using ValueType = V;
	using BitsType = B;
	constexpr ValueBitsMapper(ValueType _value) : value{ _value } {}
	constexpr ValueBitsMapper(BitsType _bits) : bits{ _bits } {}
	ValueType	value;
	BitsType	bits;
};

using DecimalBitsMapper = ValueBitsMapper<Decimal, DecimalBits>;
using FloatBitsMapper = ValueBitsMapper<float, uint32_t>;
using DoubleBitsMapper = ValueBitsMapper<double, uint64_t>;


template <typename T> constexpr bool IsInf(T _v) { return false; }
template <typename T> constexpr bool IsNaN(T _v) { return false; }
template <> constexpr bool IsInf<float>(float _v);
template <> constexpr bool IsNaN<float>(float _v);
template <> constexpr bool IsInf<double>(double _v);
template <> constexpr bool IsNaN<double>(double _v);

constexpr double NextDecimalUp(double _v, uint64_t _delta = 1);
constexpr double NextDecimalDown(double _v, uint64_t _delta = 1);
constexpr float NextDecimalUp(float _v, uint32_t _delta = 1);
constexpr float NextDecimalDown(float _v, uint32_t _delta = 1);

} // namespace maths


namespace maths
{

template <typename T> static constexpr T lowest_value = std::numeric_limits<T>::lowest();
template <typename T> static constexpr T highest_value = std::numeric_limits<T>::max();
template <typename T> static constexpr T infinity = std::numeric_limits<T>::infinity();
template <typename T> static constexpr T almost_one = one<T> - std::numeric_limits<T>::epsilon();

template <typename T> static constexpr T pi = T( 3.14159235658979323846 );

template <typename T>
constexpr T Radians(T _degrees) { return (pi<T> / T( 180 )) * _degrees; }
template <typename T>
constexpr T Degrees(T _radians) { return (T( 180 ) / pi<T>) * _radians; }

// NOTE: Maybe a Scalar<T> class could make these functions a little more specific.
template <typename T>
constexpr T Min(T _lhs, T _rhs) { return (_lhs < _rhs) ? _lhs : _rhs; }
template <typename T>
constexpr T Max(T _lhs, T _rhs) { return (_lhs > _rhs) ? _lhs : _rhs; }
template <typename T>
constexpr T Clamp(T _v, T _min, T _max) { return Min(Max(_v, _min), _max); }
template <typename T>
constexpr T SafeClamp(T _v, T _a, T _b) { return Min(Max(_v, Min(_a, _b)), Max(_a, _b)); }

template <typename T>
constexpr T Abs(T _v) { return (_v > zero<T>) ? _v : -_v; }

template <typename T>
constexpr T Lerp(T _a, T _b, float _t) { return _a*(1.f - _t) + _b*_t; }
template <typename T>
constexpr T Lerp(T _a, T _b, double _t) { return _a*(1. - _t) + _b*_t; }

constexpr float InvSqrt(float _v)
{
	// Lomont (2003) for the constant
	float half_v = 0.5f * _v;
	ValueBitsMapper<float, int32_t>	mapper(_v);
	mapper.bits = 0x5f375a86 - (mapper.bits >> 1);
	FloatBitsMapper x0{ mapper.value }, x1{ 0.f }, xold{ 0.f };
	x1.value = x0.value * (1.5f - half_v * x0.value * x0.value);
	xold.bits = x1.bits + 1;
	while ((x1.bits - x0.bits) > 1u && xold.bits != x1.bits)
	{
		xold = x0;
		x0 = x1;
		x1.value = x0.value * (1.5f - half_v * x0.value * x0.value);
	}
	if (x1.bits == x0.bits)
		return x1.value;
	else
		return x1.value * 0.5f + x0.value * 0.5f;
}

constexpr double InvSqrt(double _v)
{
	double half_v = 0.5 * _v;
	ValueBitsMapper<double, int64_t>	mapper(_v);
	mapper.bits = 0x5fe6ec85e7de30da - (mapper.bits >> 1);
	DoubleBitsMapper x0{ mapper.value }, x1{ 0. }, xold{ 0. };
	x1.value = x0.value * (1.5 - half_v * x0.value * x0.value);
	xold.bits = x1.bits + 1;
	while ((x1.bits - x0.bits) > 1u && xold.bits != x1.bits)
	{
		xold = x0;
		x0 = x1;
		x1.value = x0.value * (1.5 - half_v * x0.value * x0.value);
	}
	if (x1.bits == x0.bits)
		return x1.value;
	else
		return x1.value * 0.5 + x0.value * 0.5;
}

constexpr float Sqrt(float _v)
{
	if (_v <= 0.f || maths::IsInf(_v) || maths::IsNaN(_v)) return _v;
	return _v * InvSqrt(_v);
}

constexpr double Sqrt(double _v)
{
	if (_v <= 0.0 || maths::IsInf(_v) || maths::IsNaN(_v)) return _v;
	return _v * InvSqrt(_v);
}


bool Quadratic(Decimal _a, Decimal _b, Decimal _c, Decimal &_t0, Decimal &_t1)
{
	_t0 = _t1 = 0.0;
	double delta = double(_b) * double(_b) - 4.0 * double(_a) * double(_c);
	if (delta < 0.0) return false;
	double sqrt_delta = std::sqrt(delta);
	double q = (_b < 0._d) ?
		-.5_d * (_b - sqrt_delta) :
		-.5_d * (_b + sqrt_delta);
	_t0 = maths::Decimal(q / _a);
	_t1 = maths::Decimal(_c / q);
	if (_t1 > _t0) std::swap(_t0, _t1);
	return true;
}

} // namespace maths


namespace maths
{

template <> constexpr bool IsInf<float>(float _v)
{
	return !((FloatBitsMapper(_v).bits & ~FloatMeta<float>::sign_mask) ^ FloatMeta<float>::exponent_mask);
}
template <> constexpr bool IsInf<double>(double _v)
{
	return !((DoubleBitsMapper(_v).bits & ~FloatMeta<double>::sign_mask) ^ FloatMeta<double>::exponent_mask);
}

template <> constexpr bool IsNaN<float>(float _v)
{
	FloatBitsMapper::BitsType bits = FloatBitsMapper(_v).bits & ~FloatMeta<float>::sign_mask;
	bool exponent_check = !(bits & FloatMeta<float>::exponent_mask ^ FloatMeta<float>::exponent_mask);
	bool significand_check = (bits & FloatMeta<float>::significand_mask);
	return exponent_check && significand_check;
}
template <> constexpr bool IsNaN<double>(double _v)
{
	DoubleBitsMapper::BitsType bits = DoubleBitsMapper(_v).bits & ~FloatMeta<double>::sign_mask;
	bool exponent_check = !(bits & FloatMeta<double>::exponent_mask ^ FloatMeta<double>::exponent_mask);
	bool significand_check = (bits & FloatMeta<double>::significand_mask);
	return exponent_check && significand_check;
}


// NOTE: Even though this function allows for different-than-one deltas, it doesn't
//		 handle them well when close to zero.
constexpr double NextDecimalUp(double _v, uint64_t _delta)
{
	if (IsInf(_v) && _v > 0.)
		return _v;
	if (_v == -0.)
		_v = 0.;

	DoubleBitsMapper	mapper(_v);
	if (_v >= 0.f)
		mapper.bits += _delta;
	else
		mapper.bits -= _delta;
	return mapper.value;
}
constexpr double NextDecimalDown(double _v, uint64_t _delta)
{
	if (IsInf(_v) && _v < 0.)
		return _v;
	if (_v == 0.)
		_v = -0.;

	DoubleBitsMapper	mapper(_v);
	if (_v >= 0.f)
		mapper.bits -= _delta;
	else
		mapper.bits += _delta;
	return mapper.value;
}

constexpr float NextDecimalUp(float _v, uint32_t _delta)
{
	if (IsInf(_v) && _v > 0.f)
		return _v;
	if (_v == -0.f)
		_v = 0.f;

	FloatBitsMapper	mapper(_v);
	if (_v >= 0.f)
		mapper.bits += _delta;
	else
		mapper.bits -= _delta;
	return mapper.value;
}
constexpr float NextDecimalDown(float _v, uint32_t _delta)
{
	if (IsInf(_v) && _v < 0.f)
		return _v;
	if (_v == 0.f)
		_v = -0.f;

	FloatBitsMapper	mapper(_v);
	if (_v >= 0.f)
		mapper.bits -= _delta;
	else
		mapper.bits += _delta;
	return mapper.value;
}


} // namespace maths


#endif // __YS_MATHS_HPP__
