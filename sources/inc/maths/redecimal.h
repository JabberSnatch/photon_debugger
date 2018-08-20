#pragma once
#ifndef __YS_REDECIMAL_HPP__
#define __YS_REDECIMAL_HPP__

#include "maths/maths.h"
#include "common_macros.h"



namespace maths
{


// Short for Running Error Decimal
struct REDecimal
{
	constexpr REDecimal() :
		value{ 0._d }, low_bound{ 0._d }, high_bound{ 0._d }
#ifdef YS_REDECIMAL_HAS_PRECISE
		, precise{ 0._d }
#endif
	{}
	REDecimal(Decimal _value) :
		value{ _value }, low_bound{ NextDecimalDown(_value) }, high_bound{ NextDecimalUp(_value) }
#ifdef YS_REDECIMAL_HAS_PRECISE
		, precise{ _value }
#endif
	{ Check(); }
	REDecimal(Decimal _value, Decimal _error) :
		value{ _value },
		low_bound{ NextDecimalDown(_value - _error) },
		high_bound{ NextDecimalUp(_value + _error) }
#ifdef YS_REDECIMAL_HAS_PRECISE
		, precise{ _value }
#endif
	{ Check(); }
	REDecimal(Decimal _value, Decimal _low_bound, Decimal _high_bound) :
		value{ _value },
		low_bound{ _low_bound },
		high_bound{ _high_bound }
#ifdef YS_REDECIMAL_HAS_PRECISE
		, precise{ _value }
#endif
	{ Check(); }

	explicit operator Decimal() const { return value; }

	Decimal		value = 0._d;
	Decimal		low_bound;
	Decimal		high_bound;

	uint64_t	round_count = 0u;

	Decimal AbsoluteError() const { return high_bound - low_bound; }
	Decimal UpperBound() const { return high_bound; }
	Decimal LowerBound() const { return low_bound; }

	inline bool operator==(REDecimal const &_rhs) const { return value == _rhs.value; }

	REDecimal operator+(REDecimal const &_rhs) const;
	REDecimal operator-(REDecimal const &_rhs) const;
	REDecimal operator*(REDecimal const &_rhs) const;
	REDecimal operator/(REDecimal const &_rhs) const;

	REDecimal &operator+=(REDecimal const &_rhs);
	REDecimal &operator-=(REDecimal const &_rhs);
	REDecimal &operator*=(REDecimal const &_rhs);
	REDecimal &operator/=(REDecimal const &_rhs);

	REDecimal &operator+();
	REDecimal operator-() const;

	// Return true only if the whole error interval is below _rhs
	bool	operator < (Decimal _rhs) const;
	// Return true only if the whole error interval is above _rhs
	bool	operator > (Decimal _rhs) const;
	// Return true only if _rhs is within the error interval
	bool    operator == (Decimal _rhs) const;

	bool	operator < (REDecimal _rhs) const;
	bool	operator > (REDecimal _rhs) const;

	void Check() const;

#ifdef YS_REDECIMAL_HAS_PRECISE
	double		precise;
#endif
};


template <> struct Zero<REDecimal> { static constexpr REDecimal value{}; };

REDecimal Sqrt(REDecimal const &_v);
bool Quadratic(REDecimal const &_a, REDecimal const &_b, REDecimal const &_c,
			   REDecimal &_t0, REDecimal &_t1);

} // namespace maths

//#include "maths/redecimal.inl"

#endif // __YS_REDECIMAL_HPP__
