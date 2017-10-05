#pragma once
#ifndef __YS_HAMMERSLEY_SAMPLER_HPP__
#define __YS_HAMMERSLEY_SAMPLER_HPP__

#include "raytracer/sampler.h"

#include "maths/redecimal.h"

namespace raytracer {


class HammersleySampler :
	public Sampler
{
public:
	static maths::REDecimal vdCInverse(uint64_t const _base, uint64_t const _value);
	static maths::REDecimal RadicalInverse(uint64_t const _base, uint64_t const _value);
};


} // namespace raytracer


#endif // __YS_HAMMERSLEY_SAMPLER_HPP__
