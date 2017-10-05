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
	struct GeometryProperties
	{
		maths::Norm3f	normal;
		maths::Vec3f	dpdu, dpdv;		// Tangent, Cotangent
		maths::Norm3f	dndu, dndv;
	};

public:
	SurfaceInteraction() = default;
	SurfaceInteraction(
		maths::Point3f const &_p, maths::Vec3f const &_e,
		maths::Decimal const &_t, maths::Vec3f const &_wo,
		Shape const *_shape, maths::Point2f const &_uv,
		maths::Vec3f const &_dpdu, maths::Vec3f const &_dpdv,
		maths::Norm3f const &_dndu, maths::Norm3f const &_dndv
	) :
		position{ _p }, error_bounds{ _e },
		time { _t }, wo{ _wo },
		shape{ _shape }, uv{ _uv },
		geometry{ 
			maths::Normalized((maths::Norm3f)maths::Cross(_dpdu, _dpdv)),
			_dpdu, _dpdv, _dndu, _dndv
		},
		shading{ geometry },
		primitive{ nullptr }
	{
		if (shape != nullptr)
			if (shape->flip_normals ^ shape->swaps_handedness)
			{
				geometry.normal = -geometry.normal;
				shading.normal = -shading.normal;
			}
	}

	// NOTE: if _is_authoritative the computed shading normal is interpreted as the "real"
	//		 outside of the shape
	void SetShadingGeometry(maths::Vec3f const &_dpdu, maths::Vec3f const &_dpdv,
							maths::Norm3f const &_dndu, maths::Norm3f const &_dndv,
							bool _is_authoritative)
	{
		shading = GeometryProperties{
			maths::Normalized((maths::Norm3f)maths::Cross(_dpdu, _dpdv)),
			_dpdu, _dpdv, _dndu, _dndv
		};
		if (shape != nullptr)
			if (shape->flip_normals ^ shape->swaps_handedness)
				shading.normal = -shading.normal;
		if (_is_authoritative)
			geometry.normal = maths::FaceForward(geometry.normal, shading.normal);
		else
			shading.normal = maths::FaceForward(shading.normal, geometry.normal);
	}

	maths::Point3f		position;
	maths::Vec3f		error_bounds;
	//maths::Bounds3f		error_bounds;
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
