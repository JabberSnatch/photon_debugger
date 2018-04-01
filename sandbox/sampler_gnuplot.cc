#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "raytracer/samplers/halton_sampler.h"
#include "raytracer/samplers/random_sampler.h"


// http://www.info.univ-angers.fr/~gh/tuteurs/tutgnuplot.htm

static constexpr uint64_t kSeed = 0x0000deadbeefcafeu;
static constexpr uint64_t kSampleCountMin = 16u;
static constexpr uint64_t kSampleCountMax = 65536u;


std::ostream& operator<<(std::ostream &_ostream, maths::Vec2f const &_value)
{
	_ostream << _value.x << " " << _value.y;
	return _ostream;
}


void EvalMetrics(raytracer::Sampler &_sampler, uint64_t const _sample_count, std::ostream &_ostream)
{
	_sampler.StartPixel(maths::Vec2u{ 0u, 0u });
	for (uint64_t sample_index = 0u; sample_index < _sample_count; ++sample_index)
	{
		maths::Vec2f sample = _sampler.GetNext<2u>();
		_ostream << sample << std::endl;
		_sampler.StartNextSample();
	}
}


int main()
{
	for (uint64_t sample_count = kSampleCountMin; sample_count <= kSampleCountMax; sample_count *= 2u)
	{
		std::stringstream output_path = std::stringstream{""};
		output_path << "sampler_output" << sample_count << ".txt";
		std::ofstream output_file{ output_path.str() };

		uint64_t const dimensions = 2u;
		{
			output_file << "# HALTON SAMPLING" << std::endl;
			maths::Vec2u const tile_resolution{ 1u, 1u };
			raytracer::HaltonSampler sampler{ kSeed, sample_count, dimensions, tile_resolution };
			EvalMetrics(sampler, sample_count, output_file);
			output_file << std::endl << std::endl;
		}

		{
			output_file << "# RANDOM SAMPLING" << std::endl;
			raytracer::RandomSampler sampler{ kSeed, sample_count, dimensions };
			EvalMetrics(sampler, sample_count, output_file);
			output_file << std::endl << std::endl;
		}
	}

	std::cout << "Hello World" << std::endl;
	return 0;
}
