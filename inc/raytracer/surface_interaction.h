#pragma once
#ifndef __YS_SURFACE_INTERACTION_HPP__
#define __YS_SURFACE_INTERACTION_HPP__

#include "raytracer/raytracer.h"
#include "raytracer/shape.h"
#include "maths/vector.h"
#include "maths/point.h"


namespace raytracer {


class SurfaceInteraction
{
public:
	class GeometryProperties
	{
		friend class SurfaceInteraction;
	public:
		GeometryProperties() = default;
		GeometryProperties(maths::Norm3f const &_normal, 
						   maths::Vec3f const &_dpdu, maths::Vec3f const &_dpdv,
						   maths::Norm3f const &_dndu, maths::Norm3f const &_dndv);
		GeometryProperties(maths::Vec3f const &_dpdu, maths::Vec3f const &_dpdv,
						   maths::Norm3f const &_dndu, maths::Norm3f const &_dndv);
		maths::Norm3f const normal() const;
		maths::Norm3f const &normal_quick() const;
		maths::Vec3f const dpdu() const;
		maths::Vec3f const &dpdu_quick() const;
		maths::Vec3f const dpdv() const;
		maths::Vec3f const &dpdv_quick() const;
		maths::Norm3f const dndu() const;
		maths::Norm3f const &dndu_quick() const;
		maths::Norm3f const dndv() const;
		maths::Norm3f const &dndv_quick() const;
		void SetNormal(maths::Norm3f const &_v);
		void SetDpdu(maths::Vec3f const &_v);
		void SetDpdv(maths::Vec3f const &_v);
		void SetDndu(maths::Norm3f const &_v);
		void SetDndv(maths::Norm3f const &_v);
	private:
		maths::Norm3f	normal_;
		maths::Vec3f	dpdu_, dpdv_;		// Tangent, Cotangent
		maths::Norm3f	dndu_, dndv_;
	};
public:
	SurfaceInteraction() = default;
	SurfaceInteraction(maths::Point3f const &_p, maths::Vec3f const &_e,
					   maths::Decimal const &_t, maths::Vec3f const &_wo,
					   Shape const *_shape, maths::Point2f const &_uv,
					   maths::Vec3f const &_dpdu, maths::Vec3f const &_dpdv,
					   maths::Norm3f const &_dndu, maths::Norm3f const &_dndv);
	SurfaceInteraction(maths::Point3f const &_p, maths::Vec3f const &_e,
					   maths::Decimal const &_t, maths::Vec3f const &_wo,
					   Shape const *const _shape, maths::Point2f const &_uv,
					   GeometryProperties const &_geometry,
					   GeometryProperties const &_shading);
	// NOTE: if _is_authoritative the computed shading normal is interpreted as the "real"
	//		 outside of the shape
	void SetShadingGeometry(maths::Vec3f const &_dpdu, maths::Vec3f const &_dpdv,
							maths::Norm3f const &_dndu, maths::Norm3f const &_dndv,
							bool _is_authoritative);
public:
	maths::Point3f		position;
	maths::Vec3f		position_error;
	maths::Decimal		time;
	maths::Vec3f		wo;				// Light flow direction
	Shape const			*shape;
	maths::Point2f		uv;
	GeometryProperties	geometry;		// True geometry properties
	GeometryProperties	shading;		// Shading geometry
	Primitive const		*primitive;
};


} // namespace raytracer


#endif // __YS_SURFACE_INTERACTION_HPP__
