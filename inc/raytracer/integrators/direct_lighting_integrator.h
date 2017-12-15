#pragma once
#ifndef __YS_DIRECT_LIGHTING_INTEGRATOR_HPP__
#define __YS_DIRECT_LIGHTING_INTEGRATOR_HPP__

#include "raytracer/integrator.h"


namespace raytracer {



class DirectLightingIntegrator : public Integrator
{
public:
	void Prepare() override;
	maths::Vec3f Li(maths::Ray const &_ray,
					raytracer::SurfaceInteraction const &_hit,
					Scene const &_scene) override;
};


} // namespace raytracer


#endif // __YS_DIRECT_LIGHTING_INTEGRATOR_HPP__
