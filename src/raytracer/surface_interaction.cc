#include "raytracer/surface_interaction.h"


namespace raytracer {


SurfaceInteraction::GeometryProperties::GeometryProperties(
	maths::Norm3f const &_normal,
	maths::Vec3f const &_dpdu, maths::Vec3f const &_dpdv,
	maths::Norm3f const &_dndu, maths::Norm3f const &_dndv
) :
	normal_{ _normal }, dpdu_{ _dpdu }, dpdv_{ _dpdv }, dndu_{ _dndu }, dndv_{ _dndv }
{}

SurfaceInteraction::GeometryProperties::GeometryProperties(
	maths::Vec3f const &_dpdu, maths::Vec3f const &_dpdv,
	maths::Norm3f const &_dndu, maths::Norm3f const &_dndv
) :
	GeometryProperties(
		maths::Normalized((maths::Norm3f)maths::Cross(_dpdu, _dpdv)),
		_dpdu, _dpdv, _dndu, _dndv)
{}


maths::Norm3f const
SurfaceInteraction::GeometryProperties::normal() const
{
	return maths::Normalized(normal_);
}
maths::Norm3f const&
SurfaceInteraction::GeometryProperties::normal_quick() const
{
	return normal_;
}
maths::Vec3f const
SurfaceInteraction::GeometryProperties::dpdu() const
{
	return maths::Normalized(dpdu_);
}
maths::Vec3f const&
SurfaceInteraction::GeometryProperties::dpdu_quick() const
{
	return dpdu_;
}
maths::Vec3f const
SurfaceInteraction::GeometryProperties::dpdv() const
{
	return maths::Normalized(dpdv_);
}
maths::Vec3f const&
SurfaceInteraction::GeometryProperties::dpdv_quick() const
{
	return dpdv_;
}
maths::Norm3f const
SurfaceInteraction::GeometryProperties::dndu() const
{
	return maths::Normalized(dndu_);
}
maths::Norm3f const&
SurfaceInteraction::GeometryProperties::dndu_quick() const
{
	return dndu_;
}
maths::Norm3f const
SurfaceInteraction::GeometryProperties::dndv() const
{
	return maths::Normalized(dndv_);
}
maths::Norm3f const&
SurfaceInteraction::GeometryProperties::dndv_quick() const
{
	return dndv_;
}
void
SurfaceInteraction::GeometryProperties::SetNormal(maths::Norm3f const &_v)
{
	normal_ = _v;
}
void
SurfaceInteraction::GeometryProperties::SetDpdu(maths::Vec3f const &_v)
{
	dpdu_ = _v;
}
void
SurfaceInteraction::GeometryProperties::SetDpdv(maths::Vec3f const &_v)
{
	dpdv_ = _v;
}
void
SurfaceInteraction::GeometryProperties::SetDndu(maths::Norm3f const &_v)
{
	dndu_ = _v;
}
void
SurfaceInteraction::GeometryProperties::SetDndv(maths::Norm3f const &_v)
{
	dndv_ = _v;
}



SurfaceInteraction::SurfaceInteraction(
	maths::Point3f const &_p, maths::Vec3f const &_e,
	maths::Decimal const &_t, maths::Vec3f const &_wo,
	Shape const *_shape, maths::Point2f const &_uv,
	maths::Vec3f const &_dpdu, maths::Vec3f const &_dpdv,
	maths::Norm3f const &_dndu, maths::Norm3f const &_dndv
) :
	position{ _p }, position_error{ _e },
	time{ _t }, wo{ _wo },
	shape{ _shape }, uv{ _uv },
	geometry{ _dpdu, _dpdv, _dndu, _dndv },
	shading{ geometry },
	primitive{ nullptr }
{
	if (shape != nullptr)
	{
		if (shape->flip_normals ^ shape->swaps_handedness)
		{
			geometry.normal_ = -geometry.normal_quick();
			shading.normal_ = -shading.normal_quick();
		}
	}
}


SurfaceInteraction::SurfaceInteraction(
	maths::Point3f const &_p, maths::Vec3f const &_e,
	maths::Decimal const &_t, maths::Vec3f const &_wo,
	Shape const *const _shape, maths::Point2f const &_uv,
	GeometryProperties const &_geometry,
	GeometryProperties const &_shading
) :
	position{ _p }, position_error{ _e },
	time{ _t }, wo{ _wo },
	shape{ _shape }, uv{ _uv },
	geometry{ _geometry }, shading{ _shading },
	primitive{ nullptr }
{}


void
SurfaceInteraction::SetShadingGeometry(maths::Vec3f const &_dpdu, maths::Vec3f const &_dpdv, 
									   maths::Norm3f const &_dndu, maths::Norm3f const &_dndv, 
									   bool _is_authoritative)
{
	shading = GeometryProperties{ _dpdu, _dpdv, _dndu, _dndv };
	if (shape != nullptr)
	{
		if (shape->flip_normals ^ shape->swaps_handedness)
			shading.normal_ = -shading.normal_quick();
	}
	if (_is_authoritative)
		geometry.normal_ = maths::FaceForward(geometry.normal_quick(), shading.normal_quick());
	else
		shading.normal_ = maths::FaceForward(shading.normal_quick(), geometry.normal_quick());
}


} // namespace raytracer
