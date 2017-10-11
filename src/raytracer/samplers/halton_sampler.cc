#include "raytracer/samplers/halton_sampler.h"


#include <vector>

#include "boost/numeric/conversion/cast.hpp"

#include "core/profiler.h"



namespace raytracer {



uint64_t
HaltonSampler::InverseRadicalInverse(uint64_t const _base,
									 uint64_t const _inverse,
									 uint64_t const _source_digits_count)
{
	uint64_t result = 0u;
	uint64_t remainder = _inverse;
	for (uint64_t digit_index = 0u; digit_index < _source_digits_count; ++digit_index)
	{
		uint64_t const digit = remainder % _base;
		remainder /= _base;
		result = result * _base + inverse_permutations_()[_base - 1u][digit];
	}
	return result;
}


maths::Decimal
HaltonSampler::RadicalInverse(uint64_t const _base, uint64_t const _value)
{
	maths::Decimal const inv_base{ 1._d / static_cast<maths::Decimal>(_base) };
	uint64_t reversed_digits = 0;
	maths::Decimal inv_base_n = 1._d;
	uint64_t working_value = _value;
	uint64_t count = 0u;
	while (working_value)
	{
		uint64_t const next = working_value / _base;
		uint64_t const digit = working_value - next * _base;
		reversed_digits = reversed_digits * _base + permutations_()[_base - 1u][digit];
		inv_base_n *= inv_base;
		working_value = next;
		++count;
	}
	maths::Decimal result = static_cast<maths::Decimal>(reversed_digits) * inv_base_n;
	YS_ASSERT(result <= maths::almost_one<maths::Decimal>);
	return result;
}


uint64_t
HaltonSampler::ModularInverse(int64_t const _a, int64_t const _m)
{
	int64_t rn = _a, rn1 = _m;
	int64_t sn = 1, sn1 = 0;
	int64_t tn = 0, tn1 = 1;
	int64_t qn{ 0 }, rn2{ 0 }, sn2{ 0 }, tn2{ 0 };
	while (rn1 != 0)
	{
		qn = rn / rn1;
		rn2 = rn - qn*rn1;
		sn2 = sn - qn*sn1;
		tn2 = tn - qn*tn1;
		rn = rn1; rn1 = rn2;
		sn = sn1; sn1 = sn2;
		tn = tn1; tn1 = tn2;
	}
	if (sn < 0)
	{
		sn += _m;
	}
	YS_ASSERT(sn > 0);
	return boost::numeric_cast<uint64_t>(sn);
}



HaltonSampler::PrimesVector_t &
HaltonSampler::primes_()
{
	static PrimesVector_t result{ 2u };
	return result;
}


HaltonSampler::PermutationsVector_t &
HaltonSampler::permutations_()
{
	static PermutationsVector_t result{ {0u}, {0u, 1u} };
	return result;
}


HaltonSampler::PermutationsVector_t &
HaltonSampler::inverse_permutations_()
{
	static PermutationsVector_t result{ {0u}, {0u, 1u} };
	return result;
}


void
HaltonSampler::ExtendPrimesSequence_(PrimesVector_t &_primes,
									 uint64_t const _count)
{
	TIMED_SCOPE(HaltonSampler_ExtendPrimesSequence);
	//
	uint64_t const primes_count = boost::numeric_cast<uint64_t>(_primes.size());
	for (uint64_t prime_index = primes_count; prime_index < _count; ++prime_index)
	{
		uint64_t candidate_prime = _primes.back() + 1u;
		PrimesVector_t::const_iterator pcit = _primes.cbegin();
		uint64_t candidate_divisor = 0u;
		while ((pcit != _primes.cend()) && (candidate_divisor * 2u < candidate_prime))
		{
			candidate_divisor = *pcit;
			if (candidate_prime % candidate_divisor == 0u)
			{
				++candidate_prime;
				pcit = _primes.cbegin();
			}
			else
			{
				++pcit;
			}
		}
		_primes.emplace_back(candidate_prime);
	}
}


void
HaltonSampler::AppendNextFaurePermutations_(PermutationsVector_t &_permutations,
											PermutationsVector_t &_inverse_permutations,
											uint64_t const _next_base)
{
	// "Good Permutations for Extreme Discrepancy" Henri Faure 1991
	_permutations.reserve(_permutations.size() + _next_base);
	_inverse_permutations.reserve(_inverse_permutations.size() + _next_base);
	uint64_t const last_base = static_cast<uint64_t>(_permutations.size());
	for (uint64_t base = last_base + 1u; base <= _next_base; ++base)
	{
		_permutations.push_back(std::move(Permutation_t(base)));
		_inverse_permutations.push_back(std::move(Permutation_t(base)));
		Permutation_t &sigma = _permutations.back();
		Permutation_t &inv_sigma = _inverse_permutations.back();
		if ((base & 1u) == 0u)
		{
			// general case :
			// s.t(l) = c*s(h) + t(k)
			// l = k*b + h; h <= b-1, k <= c-1
			// here, s = t and t = I
			// t.I(l) = 2*t(h) + I(k)
			// b = c
			// k = (0 | 1)
			// if l < c => k = 0, h = l
			// if l >= c => k = 1, h = l - c
			uint64_t const c = base / 2u;
			Permutation_t const &tau = _permutations[c - 1u];
			for (uint64_t l = 0u; l < base; ++l)
			{
				uint64_t const Ik = (l < c) ? 0u : 1u;
				uint64_t const h = (l < c) ? l : (l - c);
				sigma[l] = 2u * tau[h] + Ik;
				inv_sigma[sigma[l]] = l;
			}
		}
		else
		{
			uint64_t const c = (base - 1u) / 2u;
			Permutation_t const &tau = _permutations[2u * c - 1u];
			sigma[c] = c;
			inv_sigma[c] = c;
			for (uint64_t digit = 0u; digit < c; ++digit)
			{
				uint64_t const index = digit;
				uint64_t const increment = (tau[index] < c) ? 0u : 1u;
				sigma[digit] = tau[index] + increment;
				inv_sigma[sigma[digit]] = digit;
			}
			for (uint64_t digit = c + 1; digit < base; ++digit)
			{
				uint64_t const index = digit - 1u;
				uint64_t const increment = (tau[index] < c) ? 0u : 1u;
				sigma[digit] = tau[index] + increment;
				inv_sigma[sigma[digit]] = digit;
			}
		}
	}
}


HaltonSampler::HaltonSampler(uint64_t const _seed,
									 uint64_t const _samples_per_pixel,
									 uint64_t const _dimensions_per_sample,
									 maths::Vec2u const &_tile_resolution) :
	Sampler(_seed, _samples_per_pixel, _dimensions_per_sample),
	tile_resolution_{ _tile_resolution },
	enum_exponents_{ 0u, 0u },
	sample_stride_{ 0u },
	mj_{ 0u, 0u },
	mj_modular_inverses_{ 0u, 0u }
{
	std::vector<uint64_t> &primes = primes_();
	if (boost::numeric_cast<uint64_t>(primes.capacity()) < kReservedPrimesCount_)
	{
		primes.reserve(kReservedPrimesCount_);
	}
	ExtendPrimesSequence_(primes, _dimensions_per_sample * 2u);
	AppendNextFaurePermutations_(permutations_(), inverse_permutations_(), primes_().back());
	//
	// precomputed values for Gruenschloss enumeration method
	while (tile_resolution_.x > boost::numeric_cast<uint64_t>(std::pow(2u, enum_exponents_.x)))
	{
		enum_exponents_.x++;
	}
	tile_resolution_.x = boost::numeric_cast<uint64_t>(std::pow(2u, enum_exponents_.x));
	while (tile_resolution_.y > boost::numeric_cast<uint64_t>(std::pow(3u, enum_exponents_.y)))
	{
		enum_exponents_.y++;
	}
	tile_resolution_.y = boost::numeric_cast<uint64_t>(std::pow(3u, enum_exponents_.y));
	sample_stride_ = maths::FoldProduct(tile_resolution_);
	// mj[i] = sample_stride_ / tile_resolution_[i]
	mj_.x = tile_resolution_.y;
	mj_.y = tile_resolution_.x;
	mj_modular_inverses_.x = ModularInverse(
		boost::numeric_cast<int64_t>(mj_.x),
		boost::numeric_cast<int64_t>(tile_resolution_.x));
	mj_modular_inverses_.y = ModularInverse(
		boost::numeric_cast<int64_t>(mj_.y),
		boost::numeric_cast<int64_t>(tile_resolution_.y));
	//
	LOG_INFO(tools::kChannelGeneral, "HaltonSampler actual tile resolution : " + std::to_string(tile_resolution_.x) + "; " + std::to_string(tile_resolution_.y));
}


void
HaltonSampler::Fill1DPrimarySampleVector(Sample1DContainer_t &_sample_vector,
										 uint64_t const _sample_index)
{
	maths::Vec2f const pixel_sample = SampleCurrentPixel_(_sample_index);
	_sample_vector[0] = pixel_sample.x;
	_sample_vector[1] = pixel_sample.y;
	uint64_t const value_count = boost::numeric_cast<uint64_t>(_sample_vector.size());
	Fill1DRange_(_sample_vector, _sample_index, 2u, value_count);
}


void
HaltonSampler::Fill2DPrimarySampleVector(Sample2DContainer_t &_sample_vector,
										 uint64_t const _sample_index)
{
	maths::Vec2f const pixel_sample = SampleCurrentPixel_(_sample_index);
	_sample_vector[0] = pixel_sample;
	uint64_t const value_count = boost::numeric_cast<uint64_t>(_sample_vector.size());
	Fill2DRange_(_sample_vector, _sample_index, 1u, value_count);
}


void
HaltonSampler::Fill1DSampleVector(Sample1DContainer_t &_sample_vector,
								  uint64_t const _sample_index)
{
	uint64_t const value_count = boost::numeric_cast<uint64_t>(_sample_vector.size());
	Fill1DRange_(_sample_vector, _sample_index, 0u, value_count);
}


void
HaltonSampler::Fill2DSampleVector(Sample2DContainer_t &_sample_vector,
								  uint64_t const _sample_index)
{
	uint64_t const value_count = boost::numeric_cast<uint64_t>(_sample_vector.size());
	Fill2DRange_(_sample_vector, _sample_index, 0u, value_count);
}


void
HaltonSampler::OnArrayReserved_(uint64_t const _dimension_count)
{
	LOG_INFO(tools::kChannelGeneral, "Requested primes sequence extension to " + std::to_string(_dimension_count) + " primes");
	ExtendPrimesSequence_(primes_(), _dimension_count);
	AppendNextFaurePermutations_(permutations_(), inverse_permutations_(), primes_().back());
}


void
HaltonSampler::Fill1DRange_(Sample1DContainer_t &_sample_vector, uint64_t const _sample_index,
							uint64_t const _begin, uint64_t const _end)
{
	YS_ASSERT(_end <= boost::numeric_cast<uint64_t>(_sample_vector.size()));
	for (uint64_t value_index = _begin; value_index < _end; ++value_index)
	{
		_sample_vector[value_index] = SampleDimension_(_sample_index, value_index);
	}
}


void
HaltonSampler::Fill2DRange_(Sample2DContainer_t &_sample_vector, uint64_t const _sample_index,
							uint64_t const _begin, uint64_t const _end)
{
	YS_ASSERT(_end <= boost::numeric_cast<uint64_t>(_sample_vector.size()));
	for (uint64_t value_index = _begin; value_index < _end; ++value_index)
	{
		uint64_t const first_dimension = value_index * 2u;
		_sample_vector[value_index] = maths::Vec2f{
			SampleDimension_(_sample_index, first_dimension),
			SampleDimension_(_sample_index, first_dimension + 1)
		};
	}
}


maths::Vec2f
HaltonSampler::SampleCurrentPixel_(uint64_t const _sample_index) const
{
	// Find the index in the Halton sequence corresponding to _sample_index in current_pixel()
	// Gruenschloss sample enumeration methode is used :
	// http://gruenschloss.org/sample-enum/sample-enum.pdf
	maths::Vec2u const pj{
		current_pixel().x % tile_resolution_.x, current_pixel().y % tile_resolution_.y };
	uint64_t sequence_index = _sample_index;
	if (sample_stride_ > 1u)
	{
		maths::Vec2u const lj {
			InverseRadicalInverse(2u, pj.x, enum_exponents_.x),
			InverseRadicalInverse(3u, pj.y, enum_exponents_.y)
		};
		uint64_t const pixel_first_sequence_index =
			maths::FoldSum(lj * mj_ * mj_modular_inverses_) % sample_stride_;
		sequence_index = pixel_first_sequence_index + _sample_index * sample_stride_;
	}
	//
	maths::Vec2f result{
		RadicalInverse(2u, sequence_index), RadicalInverse(3u, sequence_index)
	};
	YS_ASSERT(result.x < maths::Decimal(pj.x + 1u) / maths::Decimal(tile_resolution_.x));
	YS_ASSERT(result.y < maths::Decimal(pj.y + 1u) / maths::Decimal(tile_resolution_.y));
	YS_ASSERT(result.x >= maths::Decimal(pj.x) / maths::Decimal(tile_resolution_.x));
	YS_ASSERT(result.y >= maths::Decimal(pj.y) / maths::Decimal(tile_resolution_.y));
	result.x = result.x * maths::Decimal(tile_resolution_.x) - maths::Decimal(pj.x);
	result.y = result.y * maths::Decimal(tile_resolution_.y) - maths::Decimal(pj.y);
	return result;
}


maths::Decimal
HaltonSampler::SampleDimension_(uint64_t const _sample_index, uint64_t const _dimension) const
{
	static std::vector<uint64_t> &primes = primes_();
	maths::Decimal result = maths::infinity<maths::Decimal>;
	YS_ASSERT(primes.size() >= _dimension);
	uint64_t const prime_base = primes[_dimension];
	result = RadicalInverse(prime_base, _sample_index + 1u);
	return result;
}


} // namespace raytracers