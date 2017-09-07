#include "raytracer/integrator.h"

#include "core/logger.h"
#include "maths/ray.h"
#include "raytracer/camera.h"
#include "raytracer/film.h"
#include "raytracer/primitive.h"
#include "raytracer/sampler.h"
#include "raytracer/surface_interaction.h"
#include "globals.h"


namespace raytracer {


void
Integrator::Integrate(std::vector<Primitive*> const &_scene, maths::Decimal _t)
{
	TIMED_SCOPE(Integrator_Integrate);

	Film	&film = *film_;

	// NOTE: Because our ray is computed as if it originated from a real camera
	//		 (object -> focus point -> sensor), our picture is reversed on both dimensions.
	//		 It needs an additional "development" step, where it is reversed again.
	//		 This process is deferred to the film through the image_is_flipped bool.
	film.image_is_flipped = true;

	maths::Vec2f const inv_resolution = { 1._d / film.resolution().w, 1._d / film.resolution().h };
	for (int32_t y = 0; y < film.resolution().h; ++y)
	{
		std::cout << "row " << y << std::endl;
		for (int32_t x = 0; x < film.resolution().w; ++x)
		{
			sampler_->StartPixel({ x, y });
			maths::Vec2f const pixel_origin =
				{ static_cast<maths::Decimal>(x), static_cast<maths::Decimal>(y) };
			maths::Vec3f color_accumulator{ maths::zero<maths::Vec3f> };
			for (uint32_t sample_index = 0;
				 sample_index < sampler_->samples_per_pixel();
				 ++sample_index, sampler_->StartNextSample())
			{
				static maths::Vec3f const
					up_color{ 0._d, 0._d, 1._d }, down_color{ 0._d, 1._d, 0._d };
				maths::Vec2f const film_sample = sampler_->Get2D();
				maths::Vec2f const sample_position = pixel_origin + film_sample;
				maths::Vec2f const uv = sample_position * inv_resolution;
				maths::Vec3f color{ 0._d, 0._d, 0._d };
				maths::Ray ray = camera_->Ray(uv.u, uv.v, _t);
				raytracer::SurfaceInteraction closest_hit_info;
				for (raytracer::Primitive const *primitive : _scene)
				{
					primitive->Intersect(ray, closest_hit_info);
				}
				if (closest_hit_info.primitive != nullptr)
				{
					color = static_cast<maths::Vec3f>(closest_hit_info.shading.normal) * 0.5_d
						+ maths::Vec3f(0.5_d);
				}
				else
				{
					color = maths::Lerp(down_color, up_color,
										.5_d * maths::Normalized(ray.direction).z + .5_d);
				}
				color_accumulator += color;
			}
			maths::Vec3f const	Lo =
				color_accumulator / static_cast<maths::Decimal>(sampler_->samples_per_pixel());
			film.SetPixel(Lo, { x, y });
		}
	}
}


} // namespace raytracer
