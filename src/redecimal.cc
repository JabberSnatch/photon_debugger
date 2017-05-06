#include "redecimal.h"


namespace maths
{

REDecimal
REDecimal::operator+(REDecimal const &_rhs) const
{
	REDecimal result;
#if YS_REDECIMAL_HAS_PRECISE
	result.precise = _rhs.precise + precise;
#endif
	result.value = _rhs.value + value;
	result.low_bound = NextDecimalDown(LowerBound() + _rhs.LowerBound());
	result.high_bound = NextDecimalUp(UpperBound() + _rhs.UpperBound());
	result.Check();
	return result;
}
REDecimal
REDecimal::operator-(REDecimal const &_rhs) const
{
	REDecimal result;
#if YS_REDECIMAL_HAS_PRECISE
	result.precise = precise - _rhs.precise;
#endif
	result.value = value - _rhs.value;
	result.low_bound = NextDecimalDown(LowerBound() - _rhs.UpperBound());
	result.high_bound = NextDecimalUp(UpperBound() - _rhs.LowerBound());
	result.Check();
	return result;
}
REDecimal
REDecimal::operator*(REDecimal const &_rhs) const
{
	REDecimal result;
#if YS_REDECIMAL_HAS_PRECISE
	result.precise = _rhs.precise * precise;
#endif
	result.value = _rhs.value * value;
	Decimal prod[4] = {
		LowerBound() * _rhs.LowerBound(), UpperBound() * _rhs.LowerBound(),
		LowerBound() * _rhs.UpperBound(), UpperBound() * _rhs.UpperBound()
	};
	result.low_bound = maths::Min(maths::Min(prod[0], prod[1]), maths::Min(prod[2], prod[3]));
	result.high_bound = maths::Max(maths::Max(prod[0], prod[1]), maths::Max(prod[2], prod[3]));
	result.Check();
	return result;
}
REDecimal
REDecimal::operator/(REDecimal const &_rhs) const
{
	REDecimal result;
#if YS_REDECIMAL_HAS_PRECISE
	result.precise = precise / _rhs.precise;
#endif
	result.value = value / _rhs.value;
	if (_rhs.low_bound < 0._d && _rhs.high_bound > 0._d)
	{
		result.low_bound = -maths::infinity<Decimal>;
		result.high_bound = maths::infinity<Decimal>;
	}
	else
	{
		Decimal div[4] = {
			LowerBound() / _rhs.LowerBound(), UpperBound() / _rhs.LowerBound(),
			LowerBound() / _rhs.UpperBound(), UpperBound() / _rhs.UpperBound()
		};
		result.low_bound = maths::Min(maths::Min(div[0], div[1]), maths::Min(div[2], div[3]));
		result.high_bound = maths::Max(maths::Max(div[0], div[1]), maths::Max(div[2], div[3]));
	}
	result.Check();
	return result;
}

REDecimal&
REDecimal::operator+=(REDecimal const &_rhs)
{
#if YS_REDECIMAL_HAS_PRECISE
	precise += _rhs.precise;
#endif
	value += _rhs.value;
	low_bound = NextDecimalDown(LowerBound() + _rhs.LowerBound());
	high_bound = NextDecimalUp(UpperBound() + _rhs.UpperBound());
	Check();
	return *this;
}
REDecimal&
REDecimal::operator-=(REDecimal const &_rhs)
{
#if YS_REDECIMAL_HAS_PRECISE
	precise -= _rhs.precise;
#endif
	value -= _rhs.value;
	low_bound = NextDecimalDown(_rhs.LowerBound() - UpperBound());
	high_bound = NextDecimalUp(_rhs.UpperBound() - LowerBound());
	Check();
	return *this;
}
REDecimal&
REDecimal::operator*=(REDecimal const &_rhs)
{
#if YS_REDECIMAL_HAS_PRECISE
	precise *= _rhs.precise;
#endif
	value *= _rhs.value;
	Decimal prod[4] = {
		LowerBound() * _rhs.LowerBound(), UpperBound() * _rhs.LowerBound(),
		LowerBound() * _rhs.UpperBound(), UpperBound() * _rhs.UpperBound()
	};
	low_bound = NextDecimalDown(
		maths::Min(maths::Min(prod[0], prod[1]), maths::Min(prod[2], prod[3])));
	high_bound = NextDecimalUp(
		maths::Max(maths::Max(prod[0], prod[1]), maths::Max(prod[2], prod[3])));
	Check();
	return *this;
}
REDecimal&
REDecimal::operator/=(REDecimal const &_rhs)
{
#if YS_REDECIMAL_HAS_PRECISE
	precise /= _rhs.precise;
#endif
	value /= _rhs.value;
	if (_rhs.low_bound < 0._d && _rhs.high_bound > 0._d)
	{
		low_bound = -maths::infinity<Decimal>;
		high_bound = maths::infinity<Decimal>;
	}
	else
	{
		Decimal div[4] = {
			LowerBound() / _rhs.LowerBound(), UpperBound() / _rhs.LowerBound(),
			LowerBound() / _rhs.UpperBound(), UpperBound() / _rhs.UpperBound()
		};
		low_bound = NextDecimalDown(
			maths::Min(maths::Min(div[0], div[1]), maths::Min(div[2], div[3])));
		high_bound = NextDecimalUp(
			maths::Max(maths::Max(div[0], div[1]), maths::Max(div[2], div[3])));
	}
	Check();
	return *this;
}

REDecimal&
REDecimal::operator+()
{
	return *this;
}
REDecimal
REDecimal::operator-() const
{
	REDecimal result;
#if YS_REDECIMAL_HAS_PRECISE
	result.precise = -precise;
#endif
	result.value = -value;
	result.low_bound = -low_bound;
	result.high_bound = -high_bound;
	result.Check();
	return result;
}

bool
REDecimal::operator < (Decimal _rhs) const
{
	return low_bound < _rhs;
}
bool
REDecimal::operator > (Decimal _rhs) const
{
	return high_bound > _rhs;
}

void
REDecimal::Check() const
{
	YS_ASSERT(!std::isnan(low_bound));
	YS_ASSERT(!std::isnan(high_bound));
	if (!std::isinf(low_bound) && !std::isinf(high_bound))
		YS_ASSERT(low_bound <= high_bound);

#if YS_REDECIMAL_HAS_PRECISE
	YS_ASSERT(!std::isnan(precise));
	if (!std::isinf(precise))
		YS_ASSERT(low_bound <= precise && high_bound >= precise);
#endif
}


REDecimal
Sqrt(REDecimal const &_v)
{
	REDecimal result;
	result.value = std::sqrt(_v.value);
	result.low_bound = NextDecimalDown(result.value - std::sqrt(_v.value - _v.low_bound));
	result.high_bound = NextDecimalUp(result.value + std::sqrt(_v.high_bound - _v.value));
#if YS_REDECIMAL_HAS_PRECISE
	result.precise = std::sqrt(_v.precise);
#endif
	result.Check();
	return result;
}

bool
Quadratic(REDecimal const &_a, REDecimal const &_b, REDecimal const &_c,
		  REDecimal &_t0, REDecimal &_t1)
{
	_t0 = _t1 = REDecimal(0._d);

	// NOTE: Not really sure which one to pick right now..
	//		 I guess REDecimal Quadratic makes the most sense when error is measured across the
	//		 whole process, but it tends to be a bit over-conservative.
	//		 Do we really want to discard results just to be safe ?
//#define QUADRATIC_RUNNING_ERROR_ON_DISCRIMINANT

#ifdef QUADRATIC_RUNNING_ERROR_ON_DISCRIMINANT
	REDecimal delta = _b * _b - REDecimal(4._d) * _a * _c;
	if (delta.value < 0._d) return false;
#if YS_REDECIMAL_HAS_PRECISE
	if (delta.precise < 0.) return false;
#endif
	REDecimal sqrt_delta = Sqrt(delta);
#else
	double delta = _b.value * _b.value - 4. * _a.value * _c.value;
	if (delta < 0.) return false;
	double double_sqrt_delta = std::sqrt(delta);
	REDecimal sqrt_delta((maths::Decimal)double_sqrt_delta);
#endif

	REDecimal q = (_b.value < 0._d) ?
		REDecimal(-.5_d) * (_b - sqrt_delta) :
		REDecimal(-.5_d) * (_b + sqrt_delta);
	_t0 = q / _a;
	_t1 = _c / q;
	if (_t1.value < _t0.value) std::swap(_t0, _t1);
	return true;
}

} // namespace maths

