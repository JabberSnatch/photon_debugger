#include "raytracer/samplers/hammersley_sampler.h"


#include <vector>

#include "boost/numeric/conversion/cast.hpp"

#include "core/profiler.h"



namespace raytracer {



maths::Decimal
HammersleySampler::VDCInverse(uint64_t const _base, uint64_t const _value)
{
	maths::Decimal const inv_base{ 1._d / static_cast<maths::Decimal>(_base) };
	uint64_t reversed_digits = 0;
	maths::Decimal inv_base_n = 1._d;
	uint64_t working_value = _value;
	while (working_value)
	{
		uint64_t const next = working_value / _base;
		uint64_t const digit = working_value - next * _base;
		reversed_digits = reversed_digits * _base + digit;
		inv_base_n *= inv_base;
		working_value = next;
	}
	maths::Decimal result = static_cast<maths::Decimal>(reversed_digits) * inv_base_n;
	YS_ASSERT(result <= maths::almost_one<maths::Decimal>);
	return result;
}


HammersleySampler::PrimesVector_t &
HammersleySampler::primes_vector()
{
	static std::vector<uint64_t> result;
	return result;
}


void
HammersleySampler::ExtendPrimesSequence(PrimesVector_t &_primes_vector, uint64_t const _count)
{
	TIMED_SCOPE(HammersleySampler_ExtendPrimesSequence);
	if (_count > 0 && _primes_vector.empty())
	{
		_primes_vector.emplace_back(2u);
	}
	//
	uint64_t const primes_count = boost::numeric_cast<uint64_t>(_primes_vector.size());
	for (uint64_t prime_index = primes_count; prime_index < _count; ++prime_index)
	{
		uint64_t candidate_prime = _primes_vector.back() + 1u;
		PrimesVector_t::const_iterator pcit = _primes_vector.cbegin();
		uint64_t candidate_divisor = 0u;
		while ((pcit != _primes_vector.cend()) && (candidate_divisor * 2u < candidate_prime))
		{
			candidate_divisor = *pcit;
			if (candidate_prime % candidate_divisor == 0u)
			{
				++candidate_prime;
				pcit = _primes_vector.cbegin();
			}
			else
			{
				++pcit;
			}
		}
		_primes_vector.emplace_back(candidate_prime);
	}
}

maths::Decimal
HammersleySampler::RadicalInverse(uint64_t const _prime_index, uint64_t const _value)
{
	TIMED_SCOPE(HammersleySampler_RadicalInverse);
	static std::vector<uint64_t> &primes = primes_vector();
	maths::Decimal result = maths::infinity<maths::Decimal>;
	YS_ASSERT(primes.size() >= _prime_index);
	uint64_t const prime_base = primes[_prime_index];
	result = VDCInverse(prime_base, _value);
	YS_ASSERT(result != maths::infinity<maths::Decimal>);
	return result;
}


HammersleySampler::HammersleySampler(uint64_t const _seed,
									 uint64_t const _samples_per_pixel,
									 uint64_t const _dimensions_per_sample,
									 maths::Vec2u const &_tile_resolution) :
	Sampler(_seed, _samples_per_pixel, _dimensions_per_sample),
	tile_resolution_{ _tile_resolution }
{
	std::vector<uint64_t> &primes = primes_vector();
	if (boost::numeric_cast<uint64_t>(primes.capacity()) < kReservedPrimesCount)
	{
		primes.reserve(kReservedPrimesCount);
	}
	ExtendPrimesSequence(primes, _dimensions_per_sample);
}


void
HammersleySampler::Fill1DSampleVector(Sample1DContainer_t &_sample_vector,
									  uint64_t const _sample_index)
{
}

void
HammersleySampler::Fill2DSampleVector(Sample2DContainer_t &_sample_vector,
									  uint64_t const _sample_index)
{
}


maths::Decimal
HammersleySampler::SampleDimension_(uint64_t const _dimension, uint64_t const _sample_index) const
{
	return 0._d;
}


} // namespace raytracers