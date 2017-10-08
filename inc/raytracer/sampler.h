#pragma once
#ifndef __YS_SAMPLER_HPP__
#define __YS_SAMPLER_HPP__

#include <vector>

#include "core/rng.h"
#include "maths/maths.h"
#include "maths/point.h"


namespace raytracer
{


maths::Vec3f HemisphereMapping(maths::Vec2f const &_uv);


// Image >=> Pixel >=> Sample >=> Array >=> Dimension
// Arrays are packed by type first and dimension second.
// e.g. for a sampler with N samples per pixel
//	arrays_1D_ :
//		[ SampleVector0, SampleVector1, ... , SampleVectorN,
//		  Ext0Vector0, Ext0Vector1, ... , Ext0VectorN,
//		  Ext1Vector0, Ext1Vector1, ... , Ext1VectorN,
//		  ...
//		  ExtXVector0, ExtXvector1, ... , ExtXVectorN ]
// Any inheriting sampler is expected to follow this layout.

class Sampler
{
private:
	template <uint64_t PackSize> struct PackType { using type = void; };
	template <uint64_t PackSize> struct SampleContainer { using type = void; };
	template <uint64_t PackSize> struct SampleVectorContainer { using type = void; };
	template <> struct PackType<1u>
	{ using type = maths::Decimal; };
	template <> struct PackType<2u>
	{ using type = maths::Vec2f; };
	template <> struct SampleContainer<1u>
	{ using type = std::vector<typename PackType<1u>::type>; };
	template <> struct SampleContainer<2u>
	{ using type = std::vector<typename PackType<2u>::type>; };
	template <> struct SampleVectorContainer<1u>
	{ using type = std::vector<typename SampleContainer<1u>::type>; };
	template <> struct SampleVectorContainer<2u>
	{ using type = std::vector<typename SampleContainer<2u>::type>; };
public:
	template <uint64_t PackSize>
	using StorageType_t = typename PackType<PackSize>::type;
	template <uint64_t PackSize>
	using SampleContainer_t = typename SampleContainer<PackSize>::type;
	template <uint64_t PackSize>
	using SampleVectorContainer_t = typename SampleVectorContainer<PackSize>::type;
public:
	using Storage1D_t = StorageType_t<1u>;
	using Storage2D_t = StorageType_t<2u>;
	using Sample1DContainer_t = SampleContainer_t<1u>;
	using Sample2DContainer_t = SampleContainer_t<2u>;
	using SampleVector1DContainer_t = SampleVectorContainer_t<1u>;
	using SampleVector2DContainer_t = SampleVectorContainer_t<2u>;
	using ExtensionSizeContainer_t = std::vector<size_t>;
public:
	Sampler(uint64_t const _seed,
			uint64_t const _samples_per_pixel, uint64_t const _dimensions_per_sample);
public:
	template <uint64_t PackSize> StorageType_t<PackSize> GetNext();
	void StartPixel(maths::Point2i const &_position);
	bool StartNextSample();
	bool SetSampleNumber(uint64_t _sample_num);
public:
	template <uint64_t PackSize> SampleContainer_t<PackSize> const &GetArray(uint64_t const _size);
	template <uint64_t PackSize> void ReserveArray(uint64_t const _size);
public:
	virtual uint64_t RoundSampleCount(uint64_t const _count) { return _count; }
	virtual void Fill1DSampleVector(Sample1DContainer_t &_sample_vector,
									uint64_t const _sample_index) = 0;
	virtual void Fill2DSampleVector(Sample2DContainer_t &_sample_vector,
									uint64_t const _sample_index) = 0;
public:
	uint64_t samples_per_pixel() const { return samples_per_pixel_; }
protected:
	template <uint64_t PackSize>
	uint64_t SampleIndexFromArrayIndex(uint64_t const _array_index) const;
protected:
	core::RNG &rng() { return rng_; }
	uint64_t dimensions_per_sample() const { return dimensions_per_sample_; }
	maths::Point2i const &current_pixel() const { return current_pixel_; }
	uint64_t current_sample() const;
	uint64_t first_extension_index() const { return samples_per_pixel(); }
	template <uint64_t PackSize> uint64_t extension_count() const;
	template <uint64_t PackSize> uint64_t current_extension() const;
	template <uint64_t PackSize> uint64_t arrays_count() const;
	template <uint64_t PackSize> ExtensionSizeContainer_t const &extension_sizes() const;
private:
	template <uint64_t PackSize> StorageType_t<PackSize> overtaxed_value_();
private:
	core::RNG					rng_;
	uint64_t const				samples_per_pixel_;
	uint64_t const				dimensions_per_sample_;
	maths::Point2i				current_pixel_;
	uint64_t					current_sample_;
private:
	template <uint64_t PackSize> uint64_t &current_dimension_();
	template <uint64_t PackSize> ExtensionSizeContainer_t &extension_sizes_();
	template <uint64_t PackSize> SampleVectorContainer_t<PackSize> const &arrays_() const;
	template <uint64_t PackSize> SampleVectorContainer_t<PackSize> &arrays_();
	template <uint64_t PackSize> uint64_t const &current_extension_() const;
	template <uint64_t PackSize> uint64_t &current_extension_();
private:
	uint64_t					current_dimension_1D_, current_dimension_2D_;
	ExtensionSizeContainer_t	extension_sizes_1D_;
	ExtensionSizeContainer_t	extension_sizes_2D_;
	SampleVector1DContainer_t	arrays_1D_;
	SampleVector2DContainer_t	arrays_2D_;
	uint64_t					current_extension_1D_, current_extension_2D_;
};


} // namespace raytracer


#endif // __YS_SAMPLER_HPP__
