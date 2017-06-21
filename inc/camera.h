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
	Camera(Film const &_film, maths::Point3f const &_position, maths::Point3f const &_target,
		   maths::Vec3f const &_up, maths::Decimal const &_horizontal_fov);

	void		Expose(Scene const &_shapes, maths::Decimal _t0);
	maths::Ray	Ray(maths::Decimal _u, maths::Decimal _v, maths::Decimal _time) const;

	maths::Point3f	position;		// world position of the center of the lense
	maths::Vec3f	forward;
	maths::Vec3f	right;
	maths::Vec3f	up;
	maths::Decimal	sensor_offset;	// sensor's offset relative to the position of the lense

	Film			film;
};


} // namespace raytracer


#endif // __YS_CAMERA_HPP__
