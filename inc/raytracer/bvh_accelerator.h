#pragma once
#ifndef __YS_BVH_ACCELERATOR_HPP__
#define __YS_BVH_ACCELERATOR_HPP__

#include <vector>

#include "core/memory_region.h"
#include "raytracer/raytracer.h"
#include "maths/maths.h"

#include "raytracer/primitive.h"

#include "maths/bounds.h"
#include "maths/point.h"

namespace raytracer {


class BvhAccelerator
	: public Primitive
{
	struct BvhPrimitiveDesc;
	struct BvhNode;
	struct LinearBvhNode;

public:
	using PrimitiveArray_t = std::vector<Primitive*>;

	BvhAccelerator(PrimitiveArray_t const &_primitives,
				   uint32_t _node_max_size);
	~BvhAccelerator();

	
	BvhNode		*BuildRecursive(core::MemoryRegion &_region,
								std::vector<BvhPrimitiveDesc> &_primitive_desc,
								uint32_t _first, uint32_t _end, int &_node_count,
								PrimitiveArray_t &_ordered_primitives);

	bool	Intersect(maths::Ray &_ray, SurfaceInteraction &_hit_info) const override;
	bool	DoesIntersect(maths::Ray const &_ray) const override;

	maths::Bounds3f	WorldBounds() const override;

private:
	void		BuildLeafNode_(BvhNode *_node, 
							   maths::Bounds3f const &_bounds,
							   std::vector<BvhPrimitiveDesc> &_primitive_desc,
							   PrimitiveArray_t &_ordered_primitives,
							   uint32_t _first, uint32_t _last, uint32_t _primitive_count) const;

	uint32_t	FlattenBvhRecursive_(BvhNode const &_node, uint32_t &_offset);

	PrimitiveArray_t	primitives_;
	LinearBvhNode		*nodes_;
	uint32_t const		node_max_size_;


	struct BvhPrimitiveDesc
	{
		BvhPrimitiveDesc(uint32_t _primitive_index, maths::Bounds3f const &_bounds) :
			primitive_index{ _primitive_index }, bounds{ _bounds },
			centroid{ maths::Blend<maths::Point3f>::Do({
				{_bounds.min, 0.5_d}, {_bounds.max, 0.5_d}
			}) }
		{}
		uint32_t			primitive_index;
		maths::Bounds3f		bounds;
		maths::Point3f		centroid;
	};

	struct BvhNode
	{
		BvhNode(uint32_t _first_index, uint32_t _count, maths::Bounds3f const &_bounds);
		BvhNode(uint32_t _axis, BvhNode &_left, BvhNode &_right);

		maths::Bounds3f		bounds;
		BvhNode				*children[2];
		uint32_t			first_primitive_index, primitive_count;
		uint32_t			split_axis;
	};

	struct LinearBvhNode
	{
		maths::Bounds3f		bounds;				// 2 * 3 * 4 || 2 * 3 * 8
		union									// 4 || 8
		{
			uint32_t	first_primitive_index;
			uint32_t	right_child_offset;
		};
		uint16_t		primitive_count;		// 2
		uint8_t			split_axis;				// 1
#ifdef YS_DECIMAL_IS_DOUBLE						// 2 * 3 * 8 + 4 + 2 + 1 = 55 -> extra 9 bytes
		uint8_t			padding[9];
#else											// 2 * 3 * 4 + 4 + 2 + 1 = 31 -> extra 1 bytes
		uint8_t			padding[1];
#endif // !YS_DECIMAL_IS_DOUBLE
	};
	static_assert(sizeof(LinearBvhNode) == 32 || sizeof(LinearBvhNode) == 64);

	// Used for the implementation of surface area heuristic partition algorithm
	struct SahBucketDesc
	{
		static constexpr uint32_t kBucketCount = 12;
		uint32_t			primitive_count = 0;
		maths::Bounds3f		bounds{};
	};
};


} // namespace raytracer


#endif // __YS_BVH_ACCELERATOR_HPP__
