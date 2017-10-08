#pragma once
#ifndef __YS_HAMMERSLEY_SAMPLER_HPP__
#define __YS_HAMMERSLEY_SAMPLER_HPP__


#include <vector>


#include "maths/maths.h"
#include "maths/vector.h"
#include "raytracer/sampler.h"

namespace raytracer {


// Hammersley sequence generates [0, 1[ values for each dimension in a deterministic way.
// This implementation allows to span these values accross multiple pixels and repeat those tiles
// accross the whole image
class HammersleySampler :
	public Sampler
{
private:
	using PrimesVector_t = std::vector<uint64_t>;
public:
	static maths::Decimal VDCInverse(uint64_t const _base, uint64_t const _value);
private:
	static constexpr uint64_t kReservedPrimesCount = 1024u;
	static PrimesVector_t &primes_vector();
	static void ExtendPrimesSequence(PrimesVector_t &_primes_vector, uint64_t const _count);
	static maths::Decimal RadicalInverse(uint64_t const _prime_index, uint64_t const _value);
public:
	HammersleySampler(uint64_t const _seed,
					  uint64_t const _samples_per_pixel, uint64_t const _dimensions_per_sample,
					  maths::Vec2u const &_tile_resolution);
	void Fill1DSampleVector(Sample1DContainer_t &_sample_vector,
							uint64_t const _sample_index) override;
	void Fill2DSampleVector(Sample2DContainer_t &_sample_vector,
							uint64_t const _sample_index) override;
private:
	maths::Decimal SampleDimension_(uint64_t const _dimension, uint64_t const _sample_index) const;
private:
	maths::Vec2u tile_resolution_;
};


} // namespace raytracer


#endif // __YS_HAMMERSLEY_SAMPLER_HPP__
