#include "camera.h"

#include "globals.h"
#include "ray.h"
#include "shape.h"
#include "surface_interaction.h"
#include "logger.h"
#include "profiler.h"


namespace raytracer {

Camera::Camera(Film const &_film, maths::Point3f const &_position, maths::Point3f const &_target,
			   maths::Vec3f const &_up, maths::Decimal const &_horizontal_fov) :
	position{ _position },
	forward{ maths::Normalized(_target - _position) },
	right{ maths::Normalized(maths::Cross(forward, _up)) },
	up{ maths::Cross(right, forward) },
	film { _film }
{
	maths::Decimal	theta = maths::Radians(_horizontal_fov);
	maths::Decimal	half_width = film.dimensions().w * 0.5_d;
	// tan(theta) = half_width / sensor_offset <=> tan(theta) / half_width = 1 / sensor_offset
	// <=> half_width / tan(theta) = sensor_offset
	sensor_offset = half_width / std::tan(theta * 0.5_d);

#ifdef YS_DEBUG
	// NOTE: This could be factored in a unit test, validating ulp error for various fov and film sizes
	{
		maths::Vec2f	film_dimensions = film.dimensions();
		maths::Point3f	middle_right = position - forward * sensor_offset + right * half_width;
		maths::Point3f	middle_left = position - forward * sensor_offset - right * half_width;
		maths::Vec3f	to_right = maths::Normalized(middle_right - position);
		maths::Vec3f	to_left = maths::Normalized(middle_left - position);
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
Camera::Expose(Scene const &_shapes, maths::Decimal _t)
{
	TIMED_SCOPE(CameraExpose);

	// NOTE: Because our ray is computed as if it originated from a real camera
	//		 (object -> focus point -> sensor), our picture is reversed on both dimensions.
	//		 It needs an additional "development" step, where it is reversed again.
	//		 This process is deferred to the film through the image_is_flipped bool.
	film.image_is_flipped = true;

	for (int32_t y = 0; y < film.resolution().h; ++y)
	{
		for (int32_t x = 0; x < film.resolution().w; ++x)
		{
			maths::Decimal	u{ x / maths::Decimal(film.resolution().w - 1) };
			maths::Decimal	v{ y / maths::Decimal(film.resolution().h - 1) };

			maths::Decimal					closest_hit = maths::infinity<maths::Decimal>;
			raytracer::SurfaceInteraction	closest_hit_info;
			maths::Ray	ray = Ray(u, v, _t);
			
			for (raytracer::Shape const &shape : _shapes)
			{
				raytracer::SurfaceInteraction	hit_info;
				maths::Decimal					t_hit;
				if (shape.Intersect(ray, t_hit, hit_info))
				{
					if (t_hit < closest_hit)
					{
						closest_hit = t_hit;
						closest_hit_info = hit_info;
					}
				}
			}

			maths::Vec3f		color{ 0._d, 0._d, 0._d };
			maths::Vec3f const	up_color{ 0._d, 0._d, 1._d }, down_color{ 0._d, 1._d, 0._d };
			if (closest_hit != maths::infinity<maths::Decimal>)
				color = maths::Abs((maths::Vec3f)closest_hit_info.geometry.normal);
			else
				color = maths::Lerp(down_color, up_color, .5_d * maths::Normalized(ray.direction).z + .5_d);

			//film.SetPixel(color, { film.resolution().w - 1 - x, film.resolution().h - 1 - y });
			film.SetPixel(color, { x, y });
		}
	}
}



maths::Ray
Camera::Ray(maths::Decimal _u, maths::Decimal _v, maths::Decimal _time) const
{
	YS_ASSERT(_u <= 1._d);
	YS_ASSERT(_u >= 0._d);
	YS_ASSERT(_v <= 1._d);
	YS_ASSERT(_v >= 0._d);

	// Lower left corner is (0, 0)
	maths::Decimal	x{ _u - .5_d }, y{ _v - .5_d };
	// direction = position - film_point
	// film_point = position - forward * sensor_offset + right * film_width * x + up * film_height * y
	// direction = position - (position - forward * sensor_offset + right * film_width * x + up * film_height * y)
	// direction = forward * sensor_offset - right * film_width * x - up * film_height * y
	maths::Ray		ray{
		position,
		forward * sensor_offset - right * film.dimensions().w * x - up * film.dimensions().h * y,
		maths::infinity<maths::Decimal>,
		_time
	};

	return ray;
}

} // namespace raytracer
