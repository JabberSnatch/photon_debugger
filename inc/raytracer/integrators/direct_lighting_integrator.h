#pragma once
#ifndef __YS_DIRECT_LIGHTING_INTEGRATOR_HPP__
#define __YS_DIRECT_LIGHTING_INTEGRATOR_HPP__

#include "raytracer/integrator.h"


namespace raytracer {



class DirectLightingIntegrator : public Integrator
{
public:
	void Prepare(PrimitiveContainer_t const &_primitives, LightContainer_t const &_lights) override;
	maths::Vec3f Li(maths::Ray const &_ray,
					raytracer::SurfaceInteraction const &_hit,
					Scene const &_scene) override;
private:
	std::size_t light_sample_count_;
};


} // namespace raytracer


#endif // __YS_DIRECT_LIGHTING_INTEGRATOR_HPP__
