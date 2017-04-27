#include "maths.h"

namespace maths {

// NOTE: Even though this function allows for different-than-one deltas, it doesn't
//		 handle them well when close to zero.
double NextDecimalUp(double _v, uint64_t _delta)
{
	if (std::isinf(_v) && _v > 0.)
		return _v;
	if (_v == -0.)
		_v = 0.;
	DoubleBitsMapper	mapper(_v);
	if (_v >= 0.)		mapper.bits += _delta;
	else				mapper.bits -= _delta;
	return mapper.value;
}
double NextDecimalDown(double _v, uint64_t _delta)
{
	if (std::isinf(_v) && _v < 0.)
		return _v;
	if (_v == 0.)
		_v = -0.;
	DoubleBitsMapper	mapper(_v);
	if (_v > 0.)		mapper.bits -= _delta;
	else				mapper.bits += _delta;
	return mapper.value;
}

float NextDecimalUp(float _v, uint32_t _delta)
{
	if (std::isinf(_v) && _v > 0.f)
		return _v;
	if (_v == -0.f)
		_v = 0.f;
	FloatBitsMapper		mapper(_v);
	if (_v >= 0.f)		mapper.bits += _delta;
	else				mapper.bits -= _delta;
	return mapper.value;
}
float NextDecimalDown(float _v, uint32_t _delta)
{
	if (std::isinf(_v) && _v < 0.f)
		return _v;
	if (_v == 0.f)
		_v = -0.f;
	FloatBitsMapper		mapper(_v);
	if (_v > 0.f)		mapper.bits -= _delta;
	else				mapper.bits += _delta;
	return mapper.value;
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
	if (_t1 < _t0) std::swap(_t0, _t1);
	return true;
}

} // namespace maths



// ============================================================================
//				 /!\	BEYOND THIS POINT LIES DARKNESS   /!\
// ============================================================================

namespace crappy_legacy {
namespace maths {

template <> bool IsInf<float>(float _v)
{
	return !((::maths::FloatBitsMapper(_v).bits & ~::maths::FloatMeta<float>::sign_mask) ^ ::maths::FloatMeta<float>::exponent_mask);
}
template <> bool IsInf<double>(double _v)
{
	return !((::maths::DoubleBitsMapper(_v).bits & ~::maths::FloatMeta<double>::sign_mask) ^ ::maths::FloatMeta<double>::exponent_mask);
}

template <> bool IsNaN<float>(float _v)
{
	::maths::FloatBitsMapper::BitsType bits = ::maths::FloatBitsMapper(_v).bits & ~::maths::FloatMeta<float>::sign_mask;
	bool exponent_check = !(bits & ::maths::FloatMeta<float>::exponent_mask ^ ::maths::FloatMeta<float>::exponent_mask);
	bool significand_check = (bits & ::maths::FloatMeta<float>::significand_mask);
	return exponent_check && significand_check;
}
template <> bool IsNaN<double>(double _v)
{
	::maths::DoubleBitsMapper::BitsType bits = ::maths::DoubleBitsMapper(_v).bits & ~::maths::FloatMeta<double>::sign_mask;
	bool exponent_check = !(bits & ::maths::FloatMeta<double>::exponent_mask ^ ::maths::FloatMeta<double>::exponent_mask);
	bool significand_check = (bits & ::maths::FloatMeta<double>::significand_mask);
	return exponent_check && significand_check;
}

float InvSqrt(float _v)
{
	// Lomont (2003) for the constant
	float half_v = 0.5f * _v;
	::maths::ValueBitsMapper<float, int32_t>	mapper(_v);
	mapper.bits = 0x5f375a86 - (mapper.bits >> 1);
	::maths::FloatBitsMapper x0{ mapper.value }, x1{ 0.f }, xold{ 0.f };
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

double InvSqrt(double _v)
{
	double half_v = 0.5 * _v;
	::maths::ValueBitsMapper<double, int64_t>	mapper(_v);
	mapper.bits = 0x5fe6ec85e7de30da - (mapper.bits >> 1);
	::maths::DoubleBitsMapper x0{ mapper.value }, x1{ 0. }, xold{ 0. };
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

float Sqrt(float _v)
{
	if (_v <= 0.f || std::isinf(_v) || std::isnan(_v)) return _v;
	return _v * InvSqrt(_v);
}

double Sqrt(double _v)
{
	if (_v <= 0.0 || std::isinf(_v) || std::isnan(_v)) return _v;
	return _v * InvSqrt(_v);
}

} // namespace maths
} // namespace crappy_legacy

