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


class Integrator final
{
public:
	void Integrate(std::vector<Primitive*> const &_scene, maths::Decimal _t);
	void SetCamera(Camera *_camera) { camera_ = _camera; }
	void SetFilm(Film *_film) { film_ = _film; }
	void SetSampler(Sampler *_sampler) { sampler_ = _sampler; }
private:
	Camera *camera_ = nullptr;
	Film *film_ = nullptr;
	Sampler *sampler_ = nullptr;
};


} // namespace raytracer


#endif // __YS_INTEGRATOR_HPP__
