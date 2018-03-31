#pragma once
#ifndef __YS_RANDOM_SAMPLER_HPP__
#define __YS_RANDOM_SAMPLER_HPP__

#include "raytracer/sampler.h"


namespace raytracer
{


class RandomSampler final : public Sampler
{
public:
	RandomSampler(uint64_t const _seed,
				  uint64_t const _samples_per_pixel, uint64_t const _dimensions_per_sample);

	void Fill1DSampleVector(Sample1DContainer_t &_sample_vector,
							uint64_t const _sample_index) override;
	void Fill2DSampleVector(Sample2DContainer_t &_sample_vector,
							uint64_t const _sample_index) override;
};


} // namespace raytracer


#endif // __YS_RANDOM_SAMPLER_HPP__
