#pragma once
#ifndef __YS_REDECIMAL_HPP__
#define __YS_REDECIMAL_HPP__

#include "maths/maths.h"
#include "common_macros.h"



namespace maths
{


static constexpr Decimal machine_epsilon = std::numeric_limits<Decimal>::epsilon() * 0.5_d;
// NOTE: Higham (2002, section 3.1) bounding term for (1+-Epsilon)^n for n*Epsilon < 1
constexpr Decimal gamma(uint32_t _n) { return (_n * machine_epsilon) / (1._d - _n * machine_epsilon); }


// Short for Running Error Decimal
struct REDecimal
{
	REDecimal() :
		REDecimal(0._d)
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

	void Check() const;

#ifdef YS_REDECIMAL_HAS_PRECISE
	double		precise;
#endif
};


REDecimal Sqrt(REDecimal const &_v);
bool Quadratic(REDecimal const &_a, REDecimal const &_b, REDecimal const &_c,
			   REDecimal &_t0, REDecimal &_t1);

} // namespace maths

//#include "maths/redecimal.inl"

#endif // __YS_REDECIMAL_HPP__
