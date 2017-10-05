#include "raytracer/samplers/hammersley_sampler.h"

namespace raytracer {


maths::REDecimal
HammersleySampler::vdCInverse(uint64_t const _base, uint64_t const _value)
{
	maths::REDecimal const inv_base{ 1._d / static_cast<maths::Decimal>(_base) };
	maths::REDecimal reversed{ 0._d };
	maths::REDecimal inv_base_n = 1._d;
	uint64_t working_value = _value;
	while (working_value)
	{
		uint64_t const next = working_value / _base;
		uint64_t const digit = working_value - next * _base;
		inv_base_n *= inv_base;
		reversed += static_cast<maths::REDecimal>(digit) * inv_base_n;
		working_value = next;
	}
	return reversed;
}
maths::REDecimal
HammersleySampler::RadicalInverse(uint64_t const _base, uint64_t const _value)
{
	maths::REDecimal const inv_base{ 1._d / static_cast<maths::Decimal>(_base) };
	uint64_t reversed_digits = 0;
	maths::REDecimal inv_base_n = 1._d;
	uint64_t working_value = _value;
	while (working_value)
	{
		uint64_t const next = working_value / _base;
		uint64_t const digit = working_value - next * _base;
		reversed_digits = reversed_digits * _base + digit;
		inv_base_n *= inv_base;
		working_value = next;
	}
	return static_cast<maths::REDecimal>(reversed_digits) * inv_base_n;
}

} // namespace raytracers