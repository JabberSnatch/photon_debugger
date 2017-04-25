#pragma once
#ifndef __YS_REDECIMAL_HPP__
#define __YS_REDECIMAL_HPP__

#include "ys_assert.h"
#include "maths.h"

#ifdef _DEBUG
#ifndef YS_DECIMAL_IS_DOUBLE
#define YS_REDECIMAL_HAS_PRECISE 1
#endif
#endif

namespace maths
{

// Short for Running Error Decimal
struct REDecimal
{
	constexpr REDecimal() :
		REDecimal(0._d)
	{}
	constexpr REDecimal(Decimal _value) :
		value{ _value }, low_bound{ NextDecimalDown(_value) }, high_bound{ NextDecimalUp(_value) }
#if YS_REDECIMAL_HAS_PRECISE
		, precise{ _value }
#endif
	{ Check(); }
	constexpr REDecimal(Decimal _value, Decimal _error) :
		value{ _value },
		low_bound{ NextDecimalDown(_value - _error) },
		high_bound{ NextDecimalUp(_value + _error) }
#if YS_REDECIMAL_HAS_PRECISE
		, precise{ _value }
#endif
	{ Check(); }

	explicit operator Decimal() const { return value; }

	Decimal		value = 0._d;
	Decimal		low_bound;
	Decimal		high_bound;

	constexpr Decimal AbsoluteError() const { return high_bound - low_bound; }
	constexpr Decimal UpperBound() const { return high_bound; }
	constexpr Decimal LowerBound() const { return low_bound; }

	constexpr REDecimal operator+(REDecimal const &_rhs) const;
	constexpr REDecimal operator-(REDecimal const &_rhs) const;
	constexpr REDecimal operator*(REDecimal const &_rhs) const;
	constexpr REDecimal operator/(REDecimal const &_rhs) const;

	constexpr REDecimal &operator+=(REDecimal const &_rhs);
	constexpr REDecimal &operator-=(REDecimal const &_rhs);
	constexpr REDecimal &operator*=(REDecimal const &_rhs);
	constexpr REDecimal &operator/=(REDecimal const &_rhs);

	constexpr REDecimal &operator+();
	constexpr REDecimal operator-() const;

	constexpr void Check() const;

#if YS_REDECIMAL_HAS_PRECISE
	double		precise;
#endif
};


constexpr bool Quadratic(REDecimal const &_a, REDecimal const &_b, REDecimal const &_c,
						 REDecimal &_t0, REDecimal &_t1);

} // namespace maths

#include "redecimal.inl"

#endif // __YS_REDECIMAL_HPP__
