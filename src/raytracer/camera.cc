#include "raytracer/camera.h"

#include "globals.h"
#include "maths/ray.h"
#include "core/logger.h"
#include "core/profiler.h"
#include "raytracer/primitive.h"
#include "raytracer/surface_interaction.h"
#include "raytracer/render_context.h"
#include "api/param_set.h"

#include <iostream>


namespace raytracer {

Camera::Camera(maths::Point3f const &_position, maths::Point3f const &_target,
			   maths::Vec3f const &_up, maths::Decimal const &_horizontal_fov, Film &_film) :
	position_{ _position },
	forward_{ maths::Normalized(_target - _position) },
	right_{ maths::Normalized(maths::Cross(forward_, _up)) },
	up_{ maths::Cross(right_, forward_) },
	film_ { _film }
{
	maths::Decimal	theta = maths::Radians(_horizontal_fov);
	maths::Decimal	half_width = film_.dimensions().w * 0.5_d;
	// tan(theta) = half_width / sensor_offset <=> tan(theta) / half_width = 1 / sensor_offset
	// <=> half_width / tan(theta) = sensor_offset
	sensor_offset_ = half_width / std::tan(theta * 0.5_d);

#ifdef YS_DEBUG
	// NOTE: This could be factored in a unit test, validating ulp error for various fov and film sizes
	{
		maths::Vec2f	film_dimensions = film_.dimensions();
		maths::Point3f	middle_right = position_ - forward_ * sensor_offset_ + right_ * half_width;
		maths::Point3f	middle_left = position_ - forward_ * sensor_offset_ - right_ * half_width;
		maths::Vec3f	to_right = maths::Normalized(middle_right - position_);
		maths::Vec3f	to_left = maths::Normalized(middle_left - position_);
		maths::Decimal	std_cos = std::cos(theta);
		maths::Decimal	dot_cos = maths::Dot(to_left, to_right);
		maths::DecimalBitsMapper	std_bits{ std_cos };
		maths::DecimalBitsMapper	dot_bits{ dot_cos };
		maths::DecimalBits			ulp_offset = maths::Max(std_bits.bits, dot_bits.bits) - maths::Min(std_bits.bits, dot_bits.bits);
		//YS_ASSERT(ulp_offset < 10);
		//YS_ASSERT(maths::Abs(std_cos - dot_cos) < std::numeric_limits<maths::Decimal>::epsilon());
		if (ulp_offset > 4)
			LOG_WARNING(tools::kChannelGeneral, "A camera has fov error greater or equal than 5 ulp on cos(half_fov)");

		maths::Vec3f		left_ray = maths::Normalized(Ray(0._d, 0.5_d, 0._d).direction);
		maths::Vec3f		right_ray = maths::Normalized(Ray(1._d, 0.5_d, 0._d).direction);
		maths::Vec3f		forward_ray = maths::Normalized(Ray(0.5_d, 0.5_d, 0._d).direction);
		//YS_ASSERT(right_ray == -to_right);
		//YS_ASSERT(left_ray == -to_left);
		//YS_ASSERT(forward_ray == forward);
	}
#endif // YS_DEBUG
}

void
Camera::Expose(Scene const &_scene, maths::Decimal _t)
{
	TIMED_SCOPE(CameraExpose);

	if (!ValidateFilm_())
		return;
	Film	&film = film_;

	// NOTE: Because our ray is computed as if it originated from a real camera
	//		 (object -> focus point -> sensor), our picture is reversed on both dimensions.
	//		 It needs an additional "development" step, where it is reversed again.
	//		 This process is deferred to the film through the image_is_flipped bool.
	film.image_is_flipped = true;

	for (int32_t y = 0; y < film.resolution().h; ++y)
	{
		std::cout << "row " << y << std::endl;
		for (int32_t x = 0; x < film.resolution().w; ++x)
		{
			maths::Decimal	u{ x / maths::Decimal(film.resolution().w - 1) };
			maths::Decimal	v{ y / maths::Decimal(film.resolution().h - 1) };

			raytracer::SurfaceInteraction	closest_hit_info;
			maths::Ray	ray = Ray(u, v, _t);

			maths::Vec3f		color{ 0._d, 0._d, 0._d };
			maths::Vec3f const	up_color{ 0._d, 0._d, 1._d }, down_color{ 0._d, 1._d, 0._d };

			for (raytracer::Primitive const *primitive : _scene)
				primitive->Intersect(ray, closest_hit_info);

			if (closest_hit_info.primitive != nullptr)
				color = (maths::Vec3f)closest_hit_info.shading.normal * 0.5_d + maths::Vec3f(0.5_d);
			else
				color = maths::Lerp(down_color, up_color, .5_d * maths::Normalized(ray.direction).z + .5_d);


			//film.SetPixel(color, { film.resolution().w - 1 - x, film.resolution().h - 1 - y });
			film.SetPixel(color, { x, y });
		}
	}
}


void
Camera::WriteToFile(std::string const &_path) const
{
	if (!ValidateFilm_())
		return;

	film_.WriteToFile(_path);
}


maths::Ray
Camera::Ray(maths::Decimal _u, maths::Decimal _v, maths::Decimal _time) const
{
	YS_ASSERT(_u <= 1._d);
	YS_ASSERT(_u >= 0._d);
	YS_ASSERT(_v <= 1._d);
	YS_ASSERT(_v >= 0._d);

	if (!ValidateFilm_())
		return maths::Ray{};

	// Lower left corner is (0, 0)
	maths::Decimal	x{ _u - .5_d }, y{ _v - .5_d };
	// direction = position - film_point
	// film_point = position - forward * sensor_offset + right * film_width * x + up * film_height * y
	// direction = position - (position - forward * sensor_offset + right * film_width * x + up * film_height * y)
	// direction = forward * sensor_offset - right * film_width * x - up * film_height * y
	maths::Ray		ray{
		position_,
		forward_ * sensor_offset_ - right_ * film_.dimensions().w * x - up_ * film_.dimensions().h * y,
		maths::infinity<maths::Decimal>,
		_time
	};

	return ray;
}


bool
Camera::ValidateFilm_() const
{
	return true;

#if 0
	YS_ASSERT(film_ != nullptr);
	if (film_ == nullptr)
	{
		LOG_WARNING(tools::kChannelGeneral, "Tried using a Camera without an instance of Film. Nothing was done.");
		return false;
	}
	return true;
#endif
}


Camera*
MakeCamera(RenderContext &_context, api::ParamSet const &_params)
{
	maths::Point3f const	position = static_cast<maths::Point3f>(
		_params.FindFloat<3>("position", { 0._d, 0._d, 0._d })
	);
	maths::Point3f const	lookat = static_cast<maths::Point3f>(
		_params.FindFloat<3>("lookat", maths::Vec3f(position) + maths::Vec3f{ 0._d, 1._d, 0._d })
	);
	maths::Vec3f const		up = _params.FindFloat<3>("up", { 0._d, 0._d, 1._d });
	maths::Decimal const	fov = _params.FindFloat("fov", 60._d);

	Camera *camera = _context.AllocCamera();
	new (camera) Camera{ position, lookat, up, fov, *_context.film() };

	return camera;
}

} // namespace raytracer
