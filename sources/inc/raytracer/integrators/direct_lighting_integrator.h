#pragma once
#ifndef __YS_DIRECT_LIGHTING_INTEGRATOR_HPP__
#define __YS_DIRECT_LIGHTING_INTEGRATOR_HPP__

#include "raytracer/integrator.h"


namespace raytracer {



class DirectLightingIntegrator : public Integrator
{
public:
	DirectLightingIntegrator(Camera& _camera, Film& _film, Sampler& _sampler, uint64_t const _shadow_ray_count);
	void Prepare(PrimitiveContainer_t const &_primitives, LightContainer_t const &_lights) override;
	maths::Vec3f Li(maths::Ray const &_ray,
					raytracer::SurfaceInteraction const &_hit,
					Scene const &_scene) override;
private:
	uint64_t shadow_ray_count_;
};


} // namespace raytracer


#endif // __YS_DIRECT_LIGHTING_INTEGRATOR_HPP__
