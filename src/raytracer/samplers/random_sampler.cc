#include "raytracer/samplers/random_sampler.h"


namespace raytracer
{


RandomSampler::RandomSampler(core::RNG &_rng,
							 unsigned const _samples_per_pixel,
							 unsigned const _dimensions_per_sample) :
	Sampler(_rng, _samples_per_pixel, _dimensions_per_sample)
{}

void
RandomSampler::Fill1DSampleVector(Sample1DContainer_t &_sample_vector)
{
	for (Sample1DContainer_t::iterator dit = _sample_vector.begin();
		 dit != _sample_vector.end(); ++dit)
	{
		*dit = rng().GetDecimal();
	}
}


void
RandomSampler::Fill2DSampleVector(Sample2DContainer_t &_sample_vector)
{
	for (Sample2DContainer_t::iterator dit = _sample_vector.begin();
		 dit != _sample_vector.end(); ++dit)
	{
		*dit = { rng().GetDecimal(), rng().GetDecimal() };
	}
}


} // namespace raytracer
