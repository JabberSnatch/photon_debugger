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
	Sampler(core::RNG &_rng,
			unsigned const _samples_per_pixel, unsigned const _dimensions_per_sample);
	maths::Decimal Get1D();
	maths::Vec2f Get2D();
	void StartPixel(maths::Point2i const &_position);
	bool StartNextSample();
	bool SetSampleNumber(unsigned _sample_num);
	Sample1DContainer_t const &GetArray1D(unsigned const _size);
	Sample2DContainer_t const &GetArray2D(unsigned const _size);
	void ReserveArray1D(unsigned const _size);
	void ReserveArray2D(unsigned const _size);
	virtual unsigned RoundSampleCount(unsigned const _count) { return _count; }
	virtual void Fill1DSampleVector(Sample1DContainer_t &_sample_vector) = 0;
	virtual void Fill2DSampleVector(Sample2DContainer_t &_sample_vector) = 0;
	unsigned samples_per_pixel() const { return samples_per_pixel_; }
protected:
	core::RNG &rng() { return rng_; }
	unsigned dimensions_per_sample() const { return dimensions_per_sample_; }
	maths::Point2i const &current_pixel() const { return current_pixel_; }
	unsigned current_sample() const;
	unsigned first_extension_index() const { return samples_per_pixel(); }
	unsigned extension_1D_count() const;
	unsigned extension_2D_count() const;
	unsigned extension_1D_index() const;
	unsigned extension_2D_index() const;
	unsigned arrays_1D_count() const { return static_cast<unsigned>(arrays_1D_.size()); }
	unsigned arrays_2D_count() const { return static_cast<unsigned>(arrays_2D_.size()); }
	ExtensionSizeContainer_t const &extension_sizes_1D() const { return extension_sizes_1D_; }
	ExtensionSizeContainer_t const &extension_sizes_2D() const { return extension_sizes_2D_; }
private:
	core::RNG					&rng_;
	unsigned const				samples_per_pixel_;
	unsigned const				dimensions_per_sample_;
	maths::Point2i				current_pixel_;
	unsigned					current_sample_;
	unsigned					current_dimension_1D_, current_dimension_2D_;
	ExtensionSizeContainer_t	extension_sizes_1D_;
	ExtensionSizeContainer_t	extension_sizes_2D_;
	SampleVector1DContainer_t	arrays_1D_;
	SampleVector2DContainer_t	arrays_2D_;
	unsigned					current_extension_1D_, current_extension_2D_;
};


} // namespace raytracer


#endif // __YS_SAMPLER_HPP__
