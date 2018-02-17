#pragma once
#ifndef __YS_INTEGRATOR_HPP__
#define __YS_INTEGRATOR_HPP__

#include <vector>

#include "maths/maths.h"
#include "maths/vector.h"


namespace raytracer {

class Camera;
class Film;
class Light;
class Sampler;
class SurfaceInteraction;
class Primitive;


class Integrator
{
public:
	struct Scene
	{
		std::vector<Primitive const *> const &_primitives;
		std::vector<Light const *> const &_lights;
	};
	using PrimitiveContainer_t = std::vector<raytracer::Primitive const*>;
	using LightContainer_t = std::vector<raytracer::Light const*>;
public:
	Integrator(Camera& _camera, Film& _film, Sampler& _sampler);
	virtual ~Integrator() = default;
	virtual void Prepare(PrimitiveContainer_t const &_primitives, LightContainer_t const &_lights) = 0;
	void Integrate(Scene const &_scene, maths::Decimal _t);
	const Camera &camera() const { return camera_; }
	const Film &film() const { return film_; }
protected:
	Sampler &sampler() { return sampler_; }
private:
	virtual maths::Vec3f Li(maths::Ray const &_ray,
							raytracer::SurfaceInteraction const &_hit,
							Scene const &_scene) = 0;
private:
	Camera &camera_;
	Film &film_;
	Sampler &sampler_;
};


class NormalIntegrator final :
	public Integrator
{
public:
	NormalIntegrator(Camera& _camera, Film& _film, Sampler& _sampler, bool const _remap, bool const _absolute);
	void Prepare(PrimitiveContainer_t const &_primitives, LightContainer_t const &_lights) override;
private:
	maths::Vec3f Li(maths::Ray const &_ray,
					raytracer::SurfaceInteraction const &_hit,
					Scene const &_scene) override;
private:
	bool remap_;
	bool absolute_;
};


class AOIntegrator final :
	public Integrator
{
private:
	static constexpr maths::Vec3f kOccludedColor = { maths::zero<maths::Vec3f> };
	static constexpr maths::Vec3f kUnoccludedColor = { maths::one<maths::Vec3f> };
	static constexpr maths::Vec3f kSecondaryRaySelfHitColor = { 0._d, 0._d, 1._d };
	static constexpr maths::Vec3f kPrimaryRaySelfHitColor = { 1._d, 0._d, 0._d };
public:
	AOIntegrator(Camera& _camera, Film& _film, Sampler& _sampler, uint64_t const _sample_count, bool const _use_shading_geometry);
	void Prepare(PrimitiveContainer_t const &_primitives, LightContainer_t const &_lights) override;
private:
	maths::Vec3f Li(maths::Ray const &_ray,
					raytracer::SurfaceInteraction const &_hit,
					Scene const &_scene) override;
private:
	uint64_t sample_count_;
	bool use_shading_geometry_;
};


} // namespace raytracer


#endif // __YS_INTEGRATOR_HPP__
