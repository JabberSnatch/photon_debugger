#include "raytracer/sampler.h"

#include <algorithm>

#include "boost/numeric/conversion/cast.hpp"

#include "common_macros.h"
#include "core/logger.h"


namespace raytracer
{


maths::Vec3f
HemisphereMapping(maths::Vec2f const &_uv)
{
	maths::Decimal const z = _uv.u;
	maths::Decimal const sin_theta = std::sqrt(maths::Max(0._d, 1._d - z * z));
	maths::Decimal const phi = 2._d * maths::pi<maths::Decimal> * _uv.v;
	return { sin_theta * std::cos(phi), sin_theta * std::sin(phi), z };
}


Sampler::Sampler(uint64_t const _seed, 
				 uint64_t const _samples_per_pixel, uint64_t const _dimensions_per_sample) :
	rng_{ _seed },
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


template <uint64_t PackSize> Sampler::StorageType_t<PackSize>
Sampler::GetNext()
{
	using SampleVectorContainerPack_t = SampleVectorContainer_t<PackSize>;
	using SampleContainerPack_t = SampleContainer_t<PackSize>;
	typename SampleVectorContainerPack_t::const_iterator const sccit =
		std::next(arrays_<PackSize>().cbegin(), current_sample());
	SampleContainerPack_t const &sample_vector = *sccit;
	if (current_dimension_<PackSize>() < dimensions_per_sample())
	{
		typename SampleContainerPack_t::const_iterator const scit =
			std::next(sample_vector.cbegin(), current_dimension_<PackSize>()++);
		return *scit;
	}
	else
	{
		LOG_WARNING(tools::kChannelGeneral, "Sampler overtaxing, tried to get a dimension bigger than the number of dimensions per sample.");
		return overtaxed_value_<PackSize>();
	}
}
template Sampler::Storage1D_t Sampler::GetNext<1u>();
template Sampler::Storage2D_t Sampler::GetNext<2u>();


void
Sampler::StartPixel(maths::Point2i const &_position)
{
	current_pixel_ = _position;
	current_sample_ = 0u;
	current_dimension_1D_ = current_dimension_2D_ = 0u;
	current_extension_1D_ = current_extension_2D_ = 0u;
	//
	{
		constexpr uint64_t kPackSize = 1u;
		for (SampleVectorContainer_t<kPackSize>::iterator svit = arrays_<kPackSize>().begin();
			 svit != arrays_<kPackSize>().end(); ++svit)
		{
			uint64_t const array_index = boost::numeric_cast<uint64_t>(
				std::distance(arrays_<kPackSize>().begin(), svit));
			SampleContainer_t<kPackSize> &sample_vector = *svit;
			uint64_t const sample_index = SampleIndexFromArrayIndex<kPackSize>(array_index);
			Fill1DSampleVector(sample_vector, sample_index);
		}
	}
	{
		constexpr uint64_t kPackSize = 2u;
		for (SampleVectorContainer_t<kPackSize>::iterator svit = arrays_<kPackSize>().begin();
			 svit != arrays_<kPackSize>().end(); ++svit)
		{
			uint64_t const array_index = boost::numeric_cast<uint64_t>(
				std::distance(arrays_<kPackSize>().begin(), svit));
			SampleContainer_t<kPackSize> &sample_vector = *svit;
			uint64_t const sample_index = SampleIndexFromArrayIndex<kPackSize>(array_index);
			Fill2DSampleVector(sample_vector, sample_index);
		}
	}
}


bool
Sampler::StartNextSample()
{
	current_dimension_1D_ = current_dimension_2D_ = 0u;
	current_extension_1D_ = current_extension_2D_ = 0u;
	return ++current_sample_ < samples_per_pixel_;
}


bool
Sampler::SetSampleNumber(uint64_t _sample_num)
{
	current_dimension_1D_ = current_dimension_2D_ = 0u;
	current_extension_1D_ = current_extension_2D_ = 0u;
	current_sample_ = _sample_num;
	return current_sample_ < samples_per_pixel_;
}


template <uint64_t PackSize> Sampler::SampleContainer_t<PackSize> const &
Sampler::GetArray(uint64_t const _size)
{
	uint64_t const extension_index = current_extension<PackSize>();
	bool const current_index_is_valid = (extension_index < arrays_count<PackSize>());
	YS_ASSERT(current_index_is_valid);
	SampleContainer_t<PackSize> const &extension_array =
		*std::next(arrays_<PackSize>().cbegin(), extension_index);
	uint64_t const extension_size = boost::numeric_cast<uint64_t>(extension_array.size());
	bool const extension_size_is_valid = (extension_size == _size);
	YS_ASSERT(extension_size_is_valid);
	current_extension_<PackSize>()++;
	return extension_array;
}
template Sampler::Sample1DContainer_t const &Sampler::GetArray<1u>(uint64_t const _size);
template Sampler::Sample2DContainer_t const &Sampler::GetArray<2u>(uint64_t const _size);


template <uint64_t PackSize> void
Sampler::ReserveArray(uint64_t const _size)
{
	YS_ASSERT(_size != 0u);
	extension_sizes_<PackSize>().push_back(_size);
	SampleVectorContainer_t<PackSize> const arrays = SampleVectorContainer_t<PackSize>(
		static_cast<size_t>(samples_per_pixel()),
		SampleContainer_t<PackSize>(static_cast<size_t>(_size))
	);
	arrays_<PackSize>().insert(arrays_<PackSize>().end(), arrays.begin(), arrays.end());
}
template void Sampler::ReserveArray<1u>(uint64_t const _size);
template void Sampler::ReserveArray<2u>(uint64_t const _size);


template <uint64_t PackSize>
uint64_t Sampler::SampleIndexFromArrayIndex(uint64_t const _array_index) const
{
	uint64_t result = std::numeric_limits<uint64_t>::max();
	if (_array_index < arrays_count<PackSize>())
	{
		result = _array_index % samples_per_pixel();
	}
	else
	{
		YS_ASSERT(false);
	}
	return result;
}
template uint64_t Sampler::SampleIndexFromArrayIndex<1u>(uint64_t const _array_index) const;
template uint64_t Sampler::SampleIndexFromArrayIndex<2u>(uint64_t const _array_index) const;


uint64_t
Sampler::current_sample() const
{
	YS_ASSERT(current_sample_ < samples_per_pixel_);
	return current_sample_;
}


template <uint64_t PackSize> uint64_t
Sampler::extension_count() const
{
	return boost::numeric_cast<uint64_t>(extension_sizes<PackSize>().size());
}
template uint64_t Sampler::extension_count<1u>() const;
template uint64_t Sampler::extension_count<2u>() const;


template <uint64_t PackSize> uint64_t
Sampler::current_extension() const
{
	uint64_t const index = first_extension_index() +
		(samples_per_pixel() * current_extension_<PackSize>()) + current_sample();
	YS_ASSERT(index < arrays_count<PackSize>());
	return index;
}
template uint64_t Sampler::current_extension<1u>() const;
template uint64_t Sampler::current_extension<2u>() const;


template <uint64_t PackSize> uint64_t
Sampler::arrays_count() const
{
	return boost::numeric_cast<uint64_t>(arrays_<PackSize>().size());
}
template uint64_t Sampler::arrays_count<1u>() const;
template uint64_t Sampler::arrays_count<2u>() const;


template <> Sampler::ExtensionSizeContainer_t const &
Sampler::extension_sizes<1u>() const
{
	return extension_sizes_1D_;
}


template <> Sampler::ExtensionSizeContainer_t const &
Sampler::extension_sizes<2u>() const
{
	return extension_sizes_2D_;
}


template <> Sampler::Storage1D_t
Sampler::overtaxed_value_<1u>()
{
	return rng().GetDecimal();
}


template <> Sampler::Storage2D_t
Sampler::overtaxed_value_<2u>()
{
	return Storage2D_t{ rng().GetDecimal(), rng().GetDecimal() };
}


template <> uint64_t &
Sampler::current_dimension_<1u>()
{
	return current_dimension_1D_;
}


template <> uint64_t &
Sampler::current_dimension_<2u>()
{
	return current_dimension_2D_;
}


template <uint64_t PackSize> Sampler::ExtensionSizeContainer_t &
Sampler::extension_sizes_()
{
	return const_cast<Sampler::ExtensionSizeContainer_t &>(
		const_cast<Sampler const *>(this)->extension_sizes<PackSize>()
	);
}


template <> Sampler::SampleVectorContainer_t<1u> const &
Sampler::arrays_<1u>() const
{
	return arrays_1D_;
}


template <> Sampler::SampleVectorContainer_t<2u> const &
Sampler::arrays_<2u>() const
{
	return arrays_2D_;
}


template <uint64_t PackSize> Sampler::SampleVectorContainer_t<PackSize> &
Sampler::arrays_()
{
	return const_cast<Sampler::SampleVectorContainer_t<PackSize> &>(
		const_cast<Sampler const *>(this)->arrays_<PackSize>()
	);
}


template <> uint64_t const &
Sampler::current_extension_<1u>() const
{
	return current_extension_1D_;
}


template <> uint64_t const &
Sampler::current_extension_<2u>() const
{
	return current_extension_2D_;
}


template <uint64_t PackSize> uint64_t &
Sampler::current_extension_()
{
	return const_cast<uint64_t &>(
		const_cast<Sampler const *>(this)->current_extension_<PackSize>()
	);
}


} // namespace raytracer
