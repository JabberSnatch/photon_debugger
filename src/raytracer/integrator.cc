#include "raytracer/integrator.h"

#include "core/logger.h"
#include "maths/matrix.h"
#include "maths/ray.h"
#include "maths/transform.h"
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
	for (int64_t y = 0; y < film.resolution().h; ++y)
	{
		std::cout << "row " << y << std::endl;
		for (int64_t x = 0; x < film.resolution().w; ++x)
		{
			sampler_->StartPixel({ static_cast<uint64_t>(x), static_cast<uint64_t>(y) });
			maths::Vec2f const pixel_origin =
				{ static_cast<maths::Decimal>(x), static_cast<maths::Decimal>(y) };
			maths::Vec3f color_accumulator{ maths::zero<maths::Vec3f> };
			for (uint32_t sample_index = 0;
				 sample_index < sampler_->samples_per_pixel();
				 ++sample_index, sampler_->StartNextSample())
			{
				maths::Vec2f const film_sample = sampler_->GetNext<2u>();
				maths::Vec2f const sample_position = pixel_origin + film_sample;
				maths::Vec2f const uv = sample_position * inv_resolution;
				maths::Ray ray = camera_->Ray(uv.u, uv.v, _t);
				raytracer::SurfaceInteraction closest_hit_info;
				for (raytracer::Primitive const *primitive : _scene)
				{
					primitive->Intersect(ray, closest_hit_info);
				}
				maths::Vec3f const color{ Li(ray, closest_hit_info, _scene) };
				color_accumulator += color;
			}
			maths::Vec3f const	Li =
				color_accumulator / static_cast<maths::Decimal>(sampler_->samples_per_pixel());
			film.SetPixel(Li, { x, y });
		}
	}
}


NormalIntegrator::NormalIntegrator(bool const _remap, bool const _absolute) :
	Integrator{},
	remap_{ _remap }, absolute_{ _absolute }
{}


void
NormalIntegrator::Prepare()
{
}


maths::Vec3f
NormalIntegrator::Li(maths::Ray const &_ray,
					 raytracer::SurfaceInteraction const &_hit,
					 std::vector<Primitive*> const &_scene)
{
	static maths::Vec3f const up_color{ 0._d, 0._d, 1._d }, down_color{ 0._d, 1._d, 0._d };
	maths::Vec3f result(0._d);
	if (_hit.primitive != nullptr)
	{
		result = static_cast<maths::Vec3f>(_hit.shading.normal());
		if (remap_)
		{
			result = result * 0.5_d + maths::Vec3f(0.5_d);
		}
		if (absolute_)
		{
			result = maths::Abs(result);
		}
	}
	else
	{
		result = maths::Lerp(down_color, up_color, .5_d * maths::Normalized(_ray.direction).z + .5_d);
	}
	return result;
}


AOIntegrator::AOIntegrator(uint64_t const _sample_count) :
	Integrator(),
	sample_count_{ _sample_count }
{}


void
AOIntegrator::Prepare()
{
	sampler().ReserveArray<2u>(sample_count_);
}


maths::Vec3f
AOIntegrator::Li(maths::Ray const &_ray,
				 raytracer::SurfaceInteraction const &_hit,
				 std::vector<Primitive*> const &_scene)
{
	if (_hit.primitive != nullptr)
	{
		maths::Vec3f occlusion{ maths::zero<maths::Vec3f> };
		Sampler::Sample2DContainer_t const &samples = sampler().GetArray<2u>(sample_count_);
		for (Sampler::Sample2DContainer_t::const_iterator scit = samples.cbegin();
			 scit != samples.cend(); ++scit)
		{
			maths::Vec2f const &sample = *scit;
			maths::Vec3f const sampled_direction = HemisphereMapping(sample);
			maths::Vec3f const shading_normal{ _hit.shading.normal() };
			maths::Vec3f const wi = 
				_hit.shading.dpdu() * sampled_direction.x +
				_hit.shading.dpdv() * sampled_direction.y +
				shading_normal * sampled_direction.z;
			YS_ASSERT(maths::Dot(wi, shading_normal) >= 0._d);
			//
			maths::Decimal const d = maths::Dot(maths::Abs(shading_normal), _hit.position_error);
			maths::Vec3f const offset = (maths::Dot(wi, shading_normal) < 0) ?
				-d * shading_normal : d * shading_normal;
			maths::Point3f origin = _hit.position + offset * 2._d;
			for (int i = 0; i < 3; ++i)
			{
				if (offset[i] > 0) origin[i] = maths::NextDecimalUp(origin[i]);
				else if (offset[i] < 0) origin[i] = maths::NextDecimalDown(origin[i]);
			}
			//
			maths::Ray ray{ origin, wi, maths::infinity<maths::Decimal>, _ray.time };
			raytracer::SurfaceInteraction closest_hit_info;
			for (raytracer::Primitive const *primitive : _scene)
			{
				primitive->Intersect(ray, closest_hit_info);
			}
			if (closest_hit_info.primitive == nullptr)
			{
				occlusion += maths::one<maths::Vec3f>;
			}
			else
			{
				if (closest_hit_info.primitive == _hit.primitive)
				{
					occlusion += maths::Vec3f{ 1._d, 0._d, 0._d };
				}
			}
		}
		return occlusion / static_cast<maths::Decimal>(samples.size());
	}
	else
	{
		return maths::one<maths::Vec3f>;
	}
}


} // namespace raytracer
