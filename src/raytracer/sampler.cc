#include "raytracer/sampler.h"

#include <algorithm>

#include "common_macros.h"
#include "core/logger.h"


namespace raytracer
{


Sampler::Sampler(core::RNG &_rng, 
				 unsigned const _samples_per_pixel, unsigned const _dimensions_per_sample) :
	rng_{ _rng },
	samples_per_pixel_{ _samples_per_pixel },
	dimensions_per_sample_{ _dimensions_per_sample },
	current_pixel_{ -1, -1 },
	current_sample_{ 0u },
	current_dimension_1D_{ 0u }, current_dimension_2D_{ 0u },
	extension_sizes_1D_{}, extension_sizes_2D_{},
	arrays_1D_{
		static_cast<size_t>(samples_per_pixel_),
		Sample1DContainer_t(static_cast<size_t>(dimensions_per_sample_)) },
	arrays_2D_{
		static_cast<size_t>(samples_per_pixel_),
		Sample2DContainer_t(static_cast<size_t>(dimensions_per_sample_)) },
	current_extension_1D_{ 0u }, current_extension_2D_{ 0u }
{
}


maths::Decimal
Sampler::Get1D()
{
	SampleVector1DContainer_t::const_iterator const sccit =
		std::next(arrays_1D_.cbegin(), current_sample());
	Sample1DContainer_t const &sample_vector = *sccit;
	if (current_dimension_1D_ < dimensions_per_sample_)
	{
		Sample1DContainer_t::const_iterator const scit =
			std::next(sample_vector.cbegin(), current_dimension_1D_++);
		return *scit;
	}
	else
	{
		LOG_WARNING(tools::kChannelGeneral, "Sampler overtaxing, tried to get a dimension bigger than the number of dimensions per sample.");
		return rng_.GetDecimal();
	}
}


maths::Vec2f
Sampler::Get2D()
{
	SampleVector2DContainer_t::const_iterator const sccit =
		std::next(arrays_2D_.cbegin(), current_sample());
	Sample2DContainer_t const &sample_vector = *sccit;
	if (current_dimension_2D_ < dimensions_per_sample_)
	{
		Sample2DContainer_t::const_iterator const scit =
			std::next(sample_vector.cbegin(), current_dimension_2D_++);
		return *scit;
	}
	else
	{
		LOG_WARNING(tools::kChannelGeneral, "Sampler overtaxing, tried to get a dimension bigger than the number of dimensions per sample.");
		return { rng_.GetDecimal(), rng_.GetDecimal() };
	}
}


void
Sampler::StartPixel(maths::Point2i const &_position)
{
	current_pixel_ = _position;
	current_sample_ = 0u;
	current_dimension_1D_ = current_dimension_2D_ = 0u;
	current_extension_1D_ = current_extension_2D_ = 0u;
	//
	std::for_each(arrays_1D_.begin(), arrays_1D_.end(), [this](Sample1DContainer_t &_sample_vector) {
		Fill1DSampleVector(_sample_vector);
	});
	std::for_each(arrays_2D_.begin(), arrays_2D_.end(), [this](Sample2DContainer_t &_sample_vector) {
		Fill2DSampleVector(_sample_vector);
	});
}


bool
Sampler::StartNextSample()
{
	current_dimension_1D_ = current_dimension_2D_ = 0u;
	current_extension_1D_ = current_extension_2D_ = 0u;
	return ++current_sample_ < samples_per_pixel_;
}


bool
Sampler::SetSampleNumber(unsigned _sample_num)
{
	current_dimension_1D_ = current_dimension_2D_ = 0u;
	current_extension_1D_ = current_extension_2D_ = 0u;
	current_sample_ = _sample_num;
	return current_sample_ < samples_per_pixel_;
}


Sampler::Sample1DContainer_t const &
Sampler::GetArray1D(unsigned const _size)
{
	unsigned const extension_index = extension_1D_index();
	bool const current_index_is_valid = 
		(extension_index < arrays_1D_count());
	YS_ASSERT(current_index_is_valid);
	bool const extension_size_is_valid =
		(static_cast<unsigned>(std::next(arrays_2D_.cbegin(), extension_index)
							   ->size()) == _size);
	YS_ASSERT(extension_size_is_valid);
	SampleVector1DContainer_t::const_iterator const sccit =
		(current_index_is_valid && extension_size_is_valid) ?
		std::next(arrays_1D_.cbegin(), extension_index) : 
		arrays_1D_.cend();
	current_extension_1D_++;
	return *sccit;
}


Sampler::Sample2DContainer_t const &
Sampler::GetArray2D(unsigned const _size)
{
	unsigned const extension_index = extension_2D_index();
	bool const current_index_is_valid =
		(extension_index < arrays_2D_count());
	YS_ASSERT(current_index_is_valid);
	bool const extension_size_is_valid =
		(static_cast<unsigned>(std::next(arrays_2D_.cbegin(), extension_index)
							   ->size()) == _size);
	YS_ASSERT(extension_size_is_valid);
	SampleVector2DContainer_t::const_iterator const sccit =
		(current_index_is_valid && extension_size_is_valid) ?
		std::next(arrays_2D_.cbegin(), extension_index) :
		arrays_2D_.cend();
	current_extension_2D_++;
	return *sccit;
}


void 
Sampler::ReserveArray1D(unsigned _size)
{
	YS_ASSERT(_size != 0u);
	extension_sizes_1D_.push_back(_size);
}


void 
Sampler::ReserveArray2D(const unsigned _size)
{
	YS_ASSERT(_size != 0u);
	extension_sizes_2D_.push_back(_size);
}


unsigned
Sampler::current_sample() const
{
	YS_ASSERT(current_sample_ < samples_per_pixel_);
	return current_sample_;
}


unsigned
Sampler::extension_1D_count() const
{
	return static_cast<unsigned>(extension_sizes_1D_.size());
}


unsigned
Sampler::extension_2D_count() const
{
	return static_cast<unsigned>(extension_sizes_2D_.size());
}


unsigned
Sampler::extension_1D_index() const
{
	unsigned const index = first_extension_index() +
		(extension_1D_count() * current_extension_1D_) + current_sample();
	YS_ASSERT(index < arrays_1D_count());
	return index;
}


unsigned
Sampler::extension_2D_index() const
{
	unsigned const index = first_extension_index() +
		(extension_2D_count() * current_extension_2D_) + current_sample();
	YS_ASSERT(index < arrays_2D_count());
	return index;
}


} // namespace raytracer
