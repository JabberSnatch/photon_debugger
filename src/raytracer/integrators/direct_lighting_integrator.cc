#include "raytracer/integrators/direct_lighting_integrator.h"

#include "maths/ray.h"
#include "raytracer/light.h"
#include "raytracer/primitive.h"
#include "raytracer/sampler.h"
#include "raytracer/surface_interaction.h"

#include "common_macros.h"
#include "core/logger.h"

namespace raytracer {


maths::Decimal
PowerHeuristic(uint32_t const _nf, maths::Decimal const _pf,
			   uint32_t const _ng, maths::Decimal const _pg)
{
	maths::Decimal const f = static_cast<maths::Decimal>(_nf) * _pf;
	maths::Decimal const g = static_cast<maths::Decimal>(_ng) * _pg;
	return (f * f) / (f * f + g * g);
}


void
DirectLightingIntegrator::Prepare(PrimitiveContainer_t const &_primitives,
								  LightContainer_t const &_lights)
{
	LOG_INFO(tools::kChannelGeneral, "Preparing DirectLightingIntegartor");
	light_sample_count_ = _primitives.size() + _lights.size();
	sampler().ReserveArray<2u>(light_sample_count_);
}


maths::Vec3f
DirectLightingIntegrator::Li(maths::Ray const &_ray,
							 raytracer::SurfaceInteraction const &_hit,
							 Scene const &_scene)
{
	TIMED_SCOPE(DirectLightingIntegrator_Li);
	// hardcoded perfect diffuse material
	constexpr maths::Decimal material_pdf = 1._d / (2._d * maths::pi<maths::Decimal>);
	//
	if (_hit.primitive == nullptr)
		return maths::Vec3f{ 0.5_d, 0.5_d, 0.5_d };
	Sampler::Sample2DContainer_t const &samples = sampler().GetArray<2u>(light_sample_count_);
	Sampler::Sample2DContainer_t::const_iterator current_sample = samples.cbegin();
	maths::Vec3f result(0._d);
	for (Light const *light : _scene._lights)
	{
		{
			// NOTE: aggregate into a single get array
			// maths::Vec2f const light_sample_ksi = sampler().GetNext<2u>();
			maths::Vec2f const light_sample_ksi = *(current_sample++);
			Light::LiSample const light_sample = light->Sample(_hit, light_sample_ksi);
			maths::Vec3f const light_wi = light_sample.wi();
			maths::Point3f const origin = _hit.OffsetOriginFromErrorBounds(light_wi);
			raytracer::SurfaceInteraction shadow_hit_info{};
			maths::Ray shadow_ray{ origin, light_wi, maths::infinity<maths::Decimal>, _ray.time };
			for (Primitive const *primitive : _scene._primitives)
			{
				primitive->Intersect(shadow_ray, shadow_hit_info);
			}
			if (!shadow_hit_info.primitive)
			{
				maths::Decimal const weight = PowerHeuristic(1u, light_sample.probability,
															 1u, material_pdf);
				maths::Decimal const cos_theta =
					maths::Abs(maths::Dot(light_wi, _hit.shading.normal()));
				maths::Vec3f const contribution = light_sample.li *
					material_pdf * cos_theta * weight / light_sample.probability;
				result += contribution;
			}
			else
			{ // shadowed light, no contribution
			}
		}
		{
			// NOTE: aggregate into a single get array
			// maths::Vec2f const material_sample_ksi = sampler().GetNext<2u>();
			maths::Vec2f const material_sample_ksi = *(current_sample++);
			// hardcoded perfect diffuse material
			maths::Vec3f const material_wi = HemisphereMapping(material_sample_ksi);
			maths::Decimal const light_pdf = light->Pdf(_hit, material_wi);
			if (light_pdf > 0._d)
			{
				maths::Point3f const origin = _hit.OffsetOriginFromErrorBounds(material_wi);
				raytracer::SurfaceInteraction shadow_hit_info{};
				maths::Ray shadow_ray{ origin, material_wi,
					maths::infinity<maths::Decimal>, _ray.time };
				for (Primitive const *primitive : _scene._primitives)
				{
					primitive->Intersect(shadow_ray, shadow_hit_info);
				}
				if (!shadow_hit_info.primitive)
				{
					maths::Decimal const weight = PowerHeuristic(1u, material_pdf,
																 1u, light_pdf);
					maths::Decimal const cos_theta =
						maths::Abs(maths::Dot(material_wi, _hit.shading.normal()));
					maths::Vec3f const contribution = light->Le() *
						light_pdf * cos_theta * weight / material_pdf;
					result += contribution;
				}
				else
				{ // shadowed light, no contribution
				}
			}
			else
			{ // sampled a direction for which the light doesn't contribute
			}
		}
	}
	return result;
}


} // namespace raytracer
