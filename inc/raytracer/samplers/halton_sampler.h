#pragma once
#ifndef __YS_HALTON_SAMPLER_HPP__
#define __YS_HALTON_SAMPLER_HPP__


#include <vector>


#include "maths/maths.h"
#include "maths/vector.h"
#include "raytracer/sampler.h"

namespace raytracer {


// Halton sequences generates [0, 1[ values for each dimension in a deterministic way.
// This implementation allows to span these values accross multiple pixels and repeat those tiles
// accross the whole image
class HaltonSampler :
	public Sampler
{
private:
	using PrimesVector_t = std::vector<uint64_t>;
	using OffsetsVector_t = std::vector<uint64_t>;
	using Permutation_t = std::vector<uint64_t>;
	using PermutationsVector_t = std::vector<Permutation_t>;
public:
	static uint64_t InverseRadicalInverse(uint64_t const _base,
										  uint64_t const _inverse,
										  uint64_t const _source_digits_count);
	static maths::Decimal RadicalInverse(uint64_t const _base, uint64_t const _value);
	// actually modular multiplicative inverse
	static uint64_t ModularInverse(int64_t const _a, int64_t const _m);
private:
	static constexpr uint64_t kReservedPrimesCount_ = 1024u;
	static PrimesVector_t &primes_();
	static PermutationsVector_t &permutations_();
	static PermutationsVector_t &inverse_permutations_();
	static void ExtendPrimesSequence_(PrimesVector_t &_primes,
									  uint64_t const _count);
	static void AppendNextFaurePermutations_(PermutationsVector_t &_permutations, 
											 PermutationsVector_t &_inverse_permutations,
											 uint64_t const _next_base);
public:
	HaltonSampler(uint64_t const _seed, 
				  uint64_t const _samples_per_pixel, uint64_t const _dimensions_per_sample,
				  maths::Vec2u const &_tile_resolution);
	void Fill1DPrimarySampleVector(Sample1DContainer_t &_sample_vector,
								   uint64_t const _sample_index) override;
	void Fill2DPrimarySampleVector(Sample2DContainer_t &_sample_vector,
								   uint64_t const _sample_index) override;
	void Fill1DSampleVector(Sample1DContainer_t &_sample_vector,
							uint64_t const _sample_index) override;
	void Fill2DSampleVector(Sample2DContainer_t &_sample_vector,
							uint64_t const _sample_index) override;
private:
	void OnArrayReserved_(uint64_t const _dimension_count) override;
private:
	void Fill1DRange_(Sample1DContainer_t &_sample_vector, uint64_t const _sample_index,
					  uint64_t const _begin, uint64_t const _end);
	void Fill2DRange_(Sample2DContainer_t &_sample_vector, uint64_t const _sample_index,
					  uint64_t const _begin, uint64_t const _end);
	maths::Vec2f SampleCurrentPixel_(uint64_t const _sample_index) const;
	maths::Decimal SampleDimension_(uint64_t const _sample_index, uint64_t const _dimension) const;
private:
	maths::Vec2u tile_resolution_;
	maths::Vec2u enum_exponents_;
	uint64_t sample_stride_;
	maths::Vec2u mj_;
	maths::Vec2u mj_modular_inverses_;
};


} // namespace raytracer


#endif // __YS_HALTON_SAMPLER_HPP__
