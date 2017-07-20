#pragma once
#ifndef __YS_CAMERA_HPP__
#define __YS_CAMERA_HPP__

#include <functional>

#include "maths.h"
#include "raytracer.h"
#include "film.h"
#include "point.h"


namespace raytracer {


class Camera final
{
public:
	using Scene = std::vector<raytracer::Primitive*>;

	Camera() = delete;
	Camera(maths::Point3f const &_position, maths::Point3f const &_target,
		   maths::Vec3f const &_up, maths::Decimal const &_horizontal_fov, Film *_film = nullptr);

	void		Expose(Scene const &_scene, maths::Decimal _t0);
	void		WriteToFile(std::string const &_path) const;
	maths::Ray	Ray(maths::Decimal _u, maths::Decimal _v, maths::Decimal _time) const;

	void		SetFilm(Film *_film) { film_ = _film; }

private:
	bool		ValidateFilm_() const;

	maths::Point3f	position_;		// world position of the center of the lense
	maths::Vec3f	forward_;
	maths::Vec3f	right_;
	maths::Vec3f	up_;
	maths::Decimal	sensor_offset_;	// sensor's offset relative to the position of the lense

	Film			*film_;
};


} // namespace raytracer


#endif // __YS_CAMERA_HPP__
