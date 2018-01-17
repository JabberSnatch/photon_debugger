
#include "raytracer/integrator.h"


#include <iomanip>
#include <sstream>


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
Integrator::Integrate(Scene const &_scene, maths::Decimal _t)
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
				bool intersected = false;
				for (raytracer::Primitive const *primitive : _scene._primitives)
				{
					bool const ret_intersect = primitive->Intersect(ray, closest_hit_info);
					intersected = ret_intersect || intersected;
				}
				maths::Vec3f const color = Li(ray, closest_hit_info, _scene);
				color_accumulator += color;
			}
			maths::Vec3f const	final_color =
				color_accumulator / static_cast<maths::Decimal>(sampler_->samples_per_pixel());
			film.SetPixel(final_color, { x, y });
		}
	}
}


NormalIntegrator::NormalIntegrator(bool const _remap, bool const _absolute) :
	Integrator{},
	remap_{ _remap }, absolute_{ _absolute }
{}


void
NormalIntegrator::Prepare(PrimitiveContainer_t const &_primitives, LightContainer_t const &_lights)
{
}


maths::Vec3f
NormalIntegrator::Li(maths::Ray const &_ray,
					 raytracer::SurfaceInteraction const &_hit,
					 Scene const &_scene)
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


AOIntegrator::AOIntegrator(uint64_t const _sample_count, bool const _use_shading_geometry) :
	Integrator(),
	sample_count_{ _sample_count },
	use_shading_geometry_{ _use_shading_geometry }
{}


void
AOIntegrator::Prepare(PrimitiveContainer_t const &_primitives, LightContainer_t const &_lights)
{
	sampler().ReserveArray<2u>(sample_count_);
}


maths::Vec3f
AOIntegrator::Li(maths::Ray const &_ray,
				 raytracer::SurfaceInteraction const &_hit,
				 Scene const &_scene)
{
	TIMED_SCOPE(AOIntegrator_Li);
	if (_hit.primitive != nullptr)
	{
		maths::Vec3f occlusion{ maths::zero<maths::Vec3f> };
		Sampler::Sample2DContainer_t const &samples = sampler().GetArray<2u>(sample_count_);
		for (Sampler::Sample2DContainer_t::const_iterator scit = samples.cbegin();
			 scit != samples.cend(); ++scit)
		{
			int64_t const index = std::distance(samples.cbegin(), scit);
			maths::Vec2f const &sample = *scit;
			maths::Vec3f const sampled_direction = HemisphereMapping(sample);
			raytracer::SurfaceInteraction::GeometryProperties const &geometry =
				(!use_shading_geometry_) ? _hit.geometry : _hit.shading;
			maths::Vec3f const normal{ geometry.normal() };
			maths::Vec3f const dpdu{ geometry.dpdu() };
			maths::Vec3f const dpdv{ geometry.dpdv() };
			maths::Vec3f const wi =
				dpdu * sampled_direction.x +
				dpdv * sampled_direction.y +
				normal * sampled_direction.z;
			YS_ASSERT(maths::Dot(wi, normal) > 0._d);
			//
			maths::Point3f origin{ maths::zero<maths::Point3f> };
			bool cast_primary_ray = false;
			bool fixed_shading_normal_self_hitting = false;
			if (use_shading_geometry_)
				// TODO: refactor the raycasting part out of this block
			{	// if we're using shading geometry, we have to account for a wi pointing underneath
				// the geometry's surface
				maths::Vec3f const geometry_normal{ _hit.geometry.normal() };
				if (maths::Dot(wi, geometry_normal) >= 0._d)
				{ // easy case, just use the standard ray
					cast_primary_ray = true;
					origin =
						SurfaceInteraction::OffsetOriginFromErrorBounds(_hit.position,
																		geometry_normal,
																		_hit.position_error);
				}
				else
				{ // cast a ray from underneath the primitive surface
					maths::Point3f const secondary_ray_origin =
						SurfaceInteraction::OffsetOriginFromErrorBounds(_hit.position,
																		-geometry_normal,
																		_hit.position_error);
					maths::Ray secondary_ray{ secondary_ray_origin, wi,
						maths::infinity<maths::Decimal>, _ray.time };
					raytracer::SurfaceInteraction hit_info;
					bool intersected = false;
					for (raytracer::Primitive const *primitive : _scene._primitives)
					{
						bool const ret_intersect = primitive->Intersect(secondary_ray, hit_info);
						intersected = ret_intersect || intersected;
					}
					if (hit_info.primitive != nullptr && intersected)
					{ // found something, we might be on its inside or its outside
						maths::Vec3f const hit_geometry_normal{ hit_info.geometry.normal() };
						if (maths::Dot(hit_geometry_normal, wi) > 0._d)
						{ // inside case, compute the origin of the AO ray
							maths::Vec3f const hit_shading_normal{ hit_info.shading.normal() };
							cast_primary_ray = true;
							fixed_shading_normal_self_hitting = true;
							origin = SurfaceInteraction::OffsetOriginFromErrorBounds(
								hit_info.position, hit_shading_normal, hit_info.position_error);
						}
						else
						{ // outside case, this is a valid AO result
							if (hit_info.primitive != _hit.primitive)
							{
								occlusion += kOccludedColor;
							}
							else
							{ // this is an error
								occlusion += kSecondaryRaySelfHitColor;
								LOG_WARNING(tools::kChannelGeneral, "Secondary ray self-hit");
							}
						}
					}
					else
					{ // nothing found, this is an unoccluded result
						occlusion += kUnoccludedColor;
					}
				} // if (maths::Dot(wi, geometry_normal) < 0._d)
			}
			else
			{ // shading geometry is disabled, no risk of getting a self-hit
				cast_primary_ray = true;
				origin = SurfaceInteraction::OffsetOriginFromErrorBounds(
					_hit.position, normal, _hit.position_error);
			}
			//
			if (cast_primary_ray)
			{
				maths::Ray ray{ origin, wi, maths::infinity<maths::Decimal>, _ray.time };
				raytracer::SurfaceInteraction closest_hit_info;
				bool intersected = false;
				for (raytracer::Primitive const *primitive : _scene._primitives)
				{
					bool const ret_intersect = primitive->Intersect(ray, closest_hit_info);
					intersected = ret_intersect || intersected;
				}
				if (closest_hit_info.primitive == nullptr || !intersected)
				{
					occlusion += kUnoccludedColor;
				}
				else
				{
					if (closest_hit_info.primitive != _hit.primitive)
					{
						occlusion += kOccludedColor;
					}
					else
					{ // this is an error
						occlusion += kPrimaryRaySelfHitColor;
						if (!fixed_shading_normal_self_hitting)
						{
							LOG_WARNING(tools::kChannelGeneral, "Primary ray self-hit");
							LOG_INFO(tools::kChannelGeneral, "	wi.normal : " +
									 std::to_string(maths::Dot(wi, normal)));
							LOG_INFO(tools::kChannelGeneral, "	wi.geometry_normal : " + 
									 std::to_string(maths::Dot(wi, _hit.geometry.normal())));
							auto const precision = std::setprecision(std::numeric_limits<double>::digits10 + 1);
							std::ostringstream camera_hit_position_stream;
							camera_hit_position_stream << "	" << precision <<
								_hit.position.x << "; " <<
								_hit.position.y << "; " <<
								_hit.position.z;
							LOG_INFO(tools::kChannelGeneral, camera_hit_position_stream.str());
							std::ostringstream origin_position_stream;
							origin_position_stream << "	" << precision <<
								origin.x << "; " <<
								origin.y << "; " <<
								origin.z;
							LOG_INFO(tools::kChannelGeneral, origin_position_stream.str());
							std::ostringstream hit_position_stream;
							hit_position_stream << "	" << precision <<
								closest_hit_info.position.x << "; " <<
								closest_hit_info.position.y << "; " <<
								closest_hit_info.position.z;
							LOG_INFO(tools::kChannelGeneral, hit_position_stream.str());
						}
						else
						{
							LOG_ERROR(tools::kChannelGeneral, "Fixed shading normal but still self-hit");
						}
					}
				}
			}
		}
		return occlusion / static_cast<maths::Decimal>(samples.size());
	}
	else
	{
		return kUnoccludedColor;
	}
}


} // namespace raytracer
