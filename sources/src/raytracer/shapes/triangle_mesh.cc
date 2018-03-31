#include "raytracer/shapes/triangle_mesh.h"

#include <numeric>

#include "assimp/Importer.hpp"
#include "assimp/postprocess.h"
#include "assimp/scene.h"

#include "maths/ray.h"
#include "maths/transform.h"
#include "raytracer/shapes/triangle.h"
#include "raytracer/surface_interaction.h"


namespace raytracer {


template class TriangleMesh<InstancingPolicyClass::SharedSource>;
template class TriangleMesh<InstancingPolicyClass::Transformed>;


template <typename InstancingPolicy>
TriangleMesh<InstancingPolicy>::TriangleMesh(maths::Transform const &_world_transform,
											 bool _flip_normals,
											 TriangleMeshData const &_data) :
	Shape(_world_transform, _flip_normals),
	data_{ _data }
{}


template <typename InstancingPolicy>
TriangleMesh<InstancingPolicy>::TriangleMesh(
	maths::Transform const &_world_transform,
	bool _flip_normals,
	TriangleMesh<InstancingPolicy> const &_sibling_instance) :
	Shape(_world_transform, _flip_normals),
	data_{ _sibling_instance.data_ }
{}


template <>
bool
TriangleMesh<InstancingPolicyClass::Transformed>::Intersect(maths::Ray const &_ray,
															maths::Decimal &_tHit,
															SurfaceInteraction &_hit_info) const
{
	maths::Ray bvh_ray{ _ray };
	bool const result = data_.bvh().Intersect(bvh_ray, _hit_info);
	_tHit = (result) ? bvh_ray.tMax : _tHit;
	return result;
}

template <>
bool
TriangleMesh<InstancingPolicyClass::SharedSource>::Intersect(maths::Ray const &_ray,
															 maths::Decimal &_tHit,
															 SurfaceInteraction &_hit_info) const
{
	maths::Ray bvh_ray{ world_transform(_ray, maths::Transform::kInverse) };
	bool const result = data_.bvh().Intersect(bvh_ray, _hit_info);
	if (result)
	{
		_tHit = bvh_ray.tMax;
		_hit_info = world_transform(_hit_info, maths::Transform::kForward);
	}
	return result;
}


template <>
bool
TriangleMesh<InstancingPolicyClass::Transformed>::DoesIntersect(maths::Ray const &_ray) const
{
	return data_.bvh().DoesIntersect(_ray);
}

template <>
bool
TriangleMesh<InstancingPolicyClass::SharedSource>::DoesIntersect(maths::Ray const &_ray) const
{
	return data_.bvh().DoesIntersect(world_transform(_ray, maths::Transform::kInverse));
}


template <typename InstancingPolicy>
maths::Decimal
TriangleMesh<InstancingPolicy>::Area() const
{
	TriangleMeshData::TriangleContainer_t const &triangles = data_.triangles();
	return std::accumulate(triangles.cbegin(), triangles.cend(), 0._d,
						   [](maths::Decimal const &_acc, Triangle const* const _triangle) {
		return _acc + _triangle->Area();
	});
}


template <typename InstancingPolicy>
Shape::SurfacePoint
TriangleMesh<InstancingPolicy>::SampleSurface(maths::Vec2f const &_ksi) const
{
	YS_ASSERT(false);
	return Shape::SurfacePoint();
}


template <>
maths::Bounds3f
TriangleMesh<InstancingPolicyClass::Transformed>::ObjectBounds() const
{
	return world_transform(data_.bounds(), maths::Transform::kInverse);
}

template <>
maths::Bounds3f
TriangleMesh<InstancingPolicyClass::SharedSource>::ObjectBounds() const
{
	return data_.bounds();
}


template <>
maths::Bounds3f
TriangleMesh<InstancingPolicyClass::Transformed>::WorldBounds() const
{
	return data_.bounds();
}

template <>
maths::Bounds3f
TriangleMesh<InstancingPolicyClass::SharedSource>::WorldBounds() const
{
	return world_transform(data_.bounds(), maths::Transform::kForward);
}


} // namespace raytracer
