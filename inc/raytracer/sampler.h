#pragma once
#ifndef __YS_SAMPLER_HPP__
#define __YS_SAMPLER_HPP__

#include <vector>

#include "core/rng.h"
#include "maths/maths.h"
#include "maths/point.h"


namespace raytracer
{


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
public:
	using Sample1DContainer_t = std::vector<maths::Decimal>;
	using Sample2DContainer_t = std::vector<maths::Vec2f>;
	using SampleVector1DContainer_t = std::vector<Sample1DContainer_t>;
	using SampleVector2DContainer_t = std::vector<Sample2DContainer_t>;
	using ExtensionSizeContainer_t = std::vector<size_t>;
public:
	Sampler(uint64_t const _seed,
			uint64_t const _samples_per_pixel, uint64_t const _dimensions_per_sample);
	maths::Decimal Get1D();
	maths::Vec2f Get2D();
	void StartPixel(maths::Point2i const &_position);
	bool StartNextSample();
	bool SetSampleNumber(uint64_t _sample_num);
	Sample1DContainer_t const &GetArray1D(uint64_t const _size);
	Sample2DContainer_t const &GetArray2D(uint64_t const _size);
	void ReserveArray1D(uint64_t const _size);
	void ReserveArray2D(uint64_t const _size);
	virtual uint64_t RoundSampleCount(uint64_t const _count) { return _count; }
	virtual void Fill1DSampleVector(Sample1DContainer_t &_sample_vector) = 0;
	virtual void Fill2DSampleVector(Sample2DContainer_t &_sample_vector) = 0;
	uint64_t samples_per_pixel() const { return samples_per_pixel_; }
protected:
	core::RNG &rng() { return rng_; }
	uint64_t dimensions_per_sample() const { return dimensions_per_sample_; }
	maths::Point2i const &current_pixel() const { return current_pixel_; }
	uint64_t current_sample() const;
	uint64_t first_extension_index() const { return samples_per_pixel(); }
	uint64_t extension_1D_count() const;
	uint64_t extension_2D_count() const;
	uint64_t extension_1D_index() const;
	uint64_t extension_2D_index() const;
	uint64_t arrays_1D_count() const { return static_cast<uint64_t>(arrays_1D_.size()); }
	uint64_t arrays_2D_count() const { return static_cast<uint64_t>(arrays_2D_.size()); }
	ExtensionSizeContainer_t const &extension_sizes_1D() const { return extension_sizes_1D_; }
	ExtensionSizeContainer_t const &extension_sizes_2D() const { return extension_sizes_2D_; }
private:
	core::RNG					rng_;
	uint64_t const				samples_per_pixel_;
	uint64_t const				dimensions_per_sample_;
	maths::Point2i				current_pixel_;
	uint64_t					current_sample_;
	uint64_t					current_dimension_1D_, current_dimension_2D_;
	ExtensionSizeContainer_t	extension_sizes_1D_;
	ExtensionSizeContainer_t	extension_sizes_2D_;
	SampleVector1DContainer_t	arrays_1D_;
	SampleVector2DContainer_t	arrays_2D_;
	uint64_t					current_extension_1D_, current_extension_2D_;
};


} // namespace raytracer


#endif // __YS_SAMPLER_HPP__
