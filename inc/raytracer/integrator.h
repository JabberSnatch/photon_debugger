#pragma once
#ifndef __YS_INTEGRATOR_HPP__
#define __YS_INTEGRATOR_HPP__

#include <vector>

#include "maths/maths.h"


namespace raytracer {

class Camera;
class Film;
class Sampler;
class Primitive;
class SurfaceInteraction;


class Integrator
{
public:
	virtual ~Integrator() = default;
	virtual void Prepare() = 0;
	void Integrate(std::vector<Primitive*> const &_scene, maths::Decimal _t);
	void SetCamera(Camera *_camera) { camera_ = _camera; }
	void SetFilm(Film *_film) { film_ = _film; }
	void SetSampler(Sampler *_sampler) { sampler_ = _sampler; }
protected:
	Sampler &sampler() { return *sampler_; }
private:
	virtual maths::Vec3f Li(maths::Ray const &_ray,
							raytracer::SurfaceInteraction const &_hit,
							std::vector<Primitive*> const &_scene) = 0;
private:
	Camera *camera_ = nullptr;
	Film *film_ = nullptr;
	Sampler *sampler_ = nullptr;
};


class NormalIntegrator final :
	public Integrator
{
public:
	void Prepare() override;
private:
	maths::Vec3f Li(maths::Ray const &_ray,
					raytracer::SurfaceInteraction const &_hit,
					std::vector<Primitive*> const &_scene) override;
};


class AOIntegrator final :
	public Integrator
{
public:
	AOIntegrator(uint64_t const _sample_count);
	void Prepare() override;
private:
	maths::Vec3f Li(maths::Ray const &_ray,
					raytracer::SurfaceInteraction const &_hit,
					std::vector<Primitive*> const &_scene) override;
private:
	uint64_t sample_count_;
};


} // namespace raytracer


#endif // __YS_INTEGRATOR_HPP__
