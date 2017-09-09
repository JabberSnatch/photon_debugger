#include "raytracer/bvh_accelerator.h"

#include "globals.h"
#include "core/logger.h"
#include "core/memory_region.h"
#include "core/alloc.h"

#include "maths/ray.h"

namespace raytracer {


BvhAccelerator::BvhAccelerator(BvhAccelerator::PrimitiveArray_t const &_primitives,
							   uint32_t _node_max_size) :
	primitives_{ _primitives },
	node_max_size_{ _node_max_size }
{
	YS_ASSERT(primitives_.size() <= maths::highest_value<uint32_t>);
	YS_ASSERT(_node_max_size <= maths::highest_value<uint16_t>);
	if (primitives_.size() == 0)
	{
		LOG_WARNING(tools::kChannelGeneral, "Requested the construction of a BVH using an empty list of primitives.");
		return;
	}

	std::vector<BvhPrimitiveDesc>	primitive_desc;
	primitive_desc.reserve(primitives_.size());
	for (size_t i = 0; i < primitives_.size(); ++i)
		primitive_desc.emplace_back(static_cast<uint32_t>(i), primitives_[i]->WorldBounds());

	core::MemoryRegion					allocator{ 1024 * 1024 };
	int									node_count = 0;
	BvhAccelerator::PrimitiveArray_t	ordered_primitives;
	BvhNode								*root;
	root = BuildRecursive(allocator, primitive_desc,
						  0u, uint32_t(primitives_.size()), node_count,
						  ordered_primitives);
	primitives_.swap(ordered_primitives);

	nodes_ = core::AllocAligned<LinearBvhNode>(node_count);
	uint32_t	offset = 0u;
	FlattenBvhRecursive_(*root, offset);
}
BvhAccelerator::~BvhAccelerator()
{
	if (nodes_ != nullptr)
		core::FreeAligned(nodes_);
}

BvhAccelerator::BvhNode	*
BvhAccelerator::BuildRecursive(core::MemoryRegion &_region,
							   std::vector<BvhPrimitiveDesc> &_primitive_desc,
							   uint32_t _first, uint32_t _last, int &_node_count,
							   PrimitiveArray_t &_ordered_primitives)
{
	TIMED_SCOPE(BvhAccelerator_BuildRecursive);

	using raytracer::BvhAccelerator;
	YS_ASSERT(_first < _last);
	YS_ASSERT(_ordered_primitives.size() <= maths::highest_value<uint32_t>);

	BvhNode				*node = reinterpret_cast<BvhNode*>(_region.Alloc(sizeof(BvhNode)));
	++_node_count;

	maths::Bounds3f		bounds{};
	for (size_t i = _first; i < _last; ++i)
		bounds = maths::Union(bounds, _primitive_desc[i].bounds);

	uint32_t			primitive_count = _last - _first;
	if (primitive_count <= node_max_size_)
	{
		BuildLeafNode_(node, bounds, _primitive_desc, _ordered_primitives,
					   _first, _last, primitive_count);
		return node;
	}
	else
	{
		maths::Bounds3f	centroids_bounds{};
		for (uint32_t i = _first; i < _last; ++i)
			centroids_bounds = maths::Union(centroids_bounds, _primitive_desc[i].centroid);
		uint32_t		axis = centroids_bounds.MaximumExtent();
		uint32_t		middle = (_last + _first) / 2;
		if (centroids_bounds.max[axis] == centroids_bounds.min[axis])
		{
			BuildLeafNode_(node, bounds, _primitive_desc, _ordered_primitives,
						   _first, _last, primitive_count);
			return node;
		}
		else
		{
			SahBucketDesc	buckets[SahBucketDesc::kBucketCount];
			for (uint32_t i = _first; i < _last; ++i)
			{
				uint32_t bucket_index = uint32_t(SahBucketDesc::kBucketCount *
					centroids_bounds.Offset(_primitive_desc[i].centroid)[axis]);
				if (bucket_index == SahBucketDesc::kBucketCount)
					bucket_index = SahBucketDesc::kBucketCount - 1;
				buckets[bucket_index].primitive_count++;
				buckets[bucket_index].bounds = maths::Union(buckets[bucket_index].bounds,
															_primitive_desc[i].bounds);
			}

			maths::Decimal	lowest_cost = maths::infinity<maths::Decimal>;
			uint32_t		lowest_cost_index = maths::highest_value<uint32_t>;
			for (uint32_t i = 0; i < SahBucketDesc::kBucketCount - 1; ++i)
			{
				maths::Bounds3f		b0, b1;
				uint32_t			count0{ 0u }, count1{ 0u };
				for (uint32_t j = 0; j <= i; ++j)
				{
					b0 = maths::Union(b0, buckets[j].bounds);
					count0 += buckets[j].primitive_count;
				}
				for (uint32_t j = i + 1; j < SahBucketDesc::kBucketCount; ++j)
				{
					b1 = maths::Union(b1, buckets[j].bounds);
					count1 += buckets[j].primitive_count;
				}
				maths::Decimal	cost = .125_d + 
					(count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) / bounds.SurfaceArea();

				if (cost < lowest_cost)
				{
					lowest_cost = cost;
					lowest_cost_index = i;
				}
			}

			maths::Decimal	leaf_cost{ static_cast<maths::Decimal>(primitive_count) };
			if (primitive_count > node_max_size_ || lowest_cost < leaf_cost)
			{
				BvhPrimitiveDesc	*pivot_primitive = 
					std::partition(&_primitive_desc[_first], &_primitive_desc[_last - 1] + 1,
								   [=](BvhPrimitiveDesc const &desc)
				{
					maths::Decimal const	desc_offset = centroids_bounds.Offset(desc.centroid)[axis];
					uint32_t bucket = uint32_t(SahBucketDesc::kBucketCount * desc_offset);
					if (bucket == SahBucketDesc::kBucketCount)
						bucket = SahBucketDesc::kBucketCount - 1;
					return bucket <= lowest_cost_index;
				});

				middle = static_cast<uint32_t>(pivot_primitive - _primitive_desc.data());
			}
			else
			{
				BuildLeafNode_(node, bounds, _primitive_desc, _ordered_primitives,
							   _first, _last, primitive_count);
				return node;
			}

			BvhNode	*left = BuildRecursive(_region, _primitive_desc, _first, middle,
										   _node_count, _ordered_primitives);
			BvhNode *right = BuildRecursive(_region, _primitive_desc, middle, _last,
											_node_count, _ordered_primitives);
			new (node) BvhNode(axis, *left, *right);
			return node;
		}
	}

	YS_ASSERT(false);
	LOG_ERROR(tools::kChannelGeneral, "An error occured while building a BVH. BuildRecursive returned nullptr.");
	return nullptr;
}

bool
BvhAccelerator::Intersect(maths::Ray &_ray, SurfaceInteraction &_hit_info) const
{
	TIMED_SCOPE(BvhAccelerator_Intersect);

	bool	hit = false;
	maths::Vec3f	direction_inverse = maths::one<maths::Vec3f> / _ray.direction;
	// We will use these values as indices in Ray::DoesIntersect, which is why they are int.
	maths::Vector<int, 3>	is_negative{
		direction_inverse.x < 0._d, direction_inverse.y < 0._d, direction_inverse.z < 0._d
	};
	
	uint32_t	to_visit_offset{ 0 }, current_node_index{ 0 };
	uint32_t	nodes_to_visit[64];
	for (;;)
	{
		LinearBvhNode const &node = nodes_[current_node_index];
		if (_ray.DoesIntersect(node.bounds, direction_inverse, is_negative))
		{
			if (node.primitive_count > 0)
			{
				for (uint16_t i = 0; i < node.primitive_count; ++i)
					hit |= primitives_[node.first_primitive_index + i]->Intersect(_ray, _hit_info);
				if (to_visit_offset == 0) break;
				current_node_index = nodes_to_visit[--to_visit_offset];
			}
			else
			{
				if (is_negative[node.split_axis])
				{
					nodes_to_visit[to_visit_offset++] = current_node_index + 1;
					current_node_index = node.right_child_offset;
				}
				else
				{
					nodes_to_visit[to_visit_offset++] = node.right_child_offset;
					current_node_index = current_node_index + 1;
				}
			}
		}
		else
		{
			if (to_visit_offset == 0) break;
			current_node_index = nodes_to_visit[--to_visit_offset];
		}
	}

	return hit;
}

bool
BvhAccelerator::DoesIntersect(maths::Ray const &_ray) const
{
	maths::Vec3f	direction_inverse = maths::one<maths::Vec3f> / _ray.direction;
	// We will use these values as indices in Ray::DoesIntersect, which is why they are int.
	maths::Vector<int, 3>	is_negative{
		direction_inverse.x < 0._d, direction_inverse.y < 0._d, direction_inverse.z < 0._d
	};

	uint32_t	to_visit_offset{ 0 }, current_node_index{ 0 };
	uint32_t	nodes_to_visit[64];
	for (;;)
	{
		LinearBvhNode const &node = nodes_[current_node_index];
		if (_ray.DoesIntersect(node.bounds, direction_inverse, is_negative))
		{
			if (node.primitive_count > 0)
			{
				for (uint16_t i = 0; i < node.primitive_count; ++i)
					if (primitives_[node.first_primitive_index + i]->DoesIntersect(_ray))
						return true;
				if (to_visit_offset == 0) break;
				current_node_index = nodes_to_visit[--to_visit_offset];
			}
			else
			{
				if (is_negative[node.split_axis])
				{
					nodes_to_visit[to_visit_offset++] = current_node_index + 1;
					current_node_index = node.right_child_offset;
				}
				else
				{
					nodes_to_visit[to_visit_offset++] = node.right_child_offset;
					current_node_index = current_node_index + 1;
				}
			}
		}
		else
		{
			if (to_visit_offset == 0) break;
			current_node_index = nodes_to_visit[--to_visit_offset];
		}
	}

	return false;
}


maths::Bounds3f
BvhAccelerator::WorldBounds() const
{
	maths::Bounds3f		world_bounds{};
	for (auto &&primitive : primitives_)
		world_bounds = maths::Union(world_bounds, primitive->WorldBounds());
	return world_bounds;
}


void
BvhAccelerator::BuildLeafNode_(BvhNode *_node,
							   maths::Bounds3f const &_bounds,
							   std::vector<BvhPrimitiveDesc> &_primitive_desc,
							   PrimitiveArray_t &_ordered_primitives,
							   uint32_t _first, uint32_t _last, uint32_t _primitive_count) const
{
	uint32_t	first_primitive_index = static_cast<uint32_t>(_ordered_primitives.size());
	for (uint32_t i = _first; i < _last; ++i)
	{
		uint32_t	primitive_index = _primitive_desc[i].primitive_index;
		_ordered_primitives.push_back(primitives_[primitive_index]);
	}
	new (_node) BvhNode(first_primitive_index, _primitive_count, _bounds);
}

uint32_t
BvhAccelerator::FlattenBvhRecursive_(BvhNode const &_node, uint32_t &_offset)
{
	LinearBvhNode &linear_node = nodes_[_offset];
	linear_node.bounds = _node.bounds;
	linear_node.primitive_count = _node.primitive_count;
	
	uint32_t	self_offset = _offset++;
	if (_node.primitive_count > 0)
		linear_node.first_primitive_index = _node.first_primitive_index;
	else
	{
		linear_node.split_axis = _node.split_axis;
		FlattenBvhRecursive_(*_node.children[0], _offset);
		linear_node.right_child_offset = FlattenBvhRecursive_(*_node.children[1], _offset);
	}

	return self_offset;
}


BvhAccelerator::BvhNode::BvhNode(uint32_t _first_index,
								 uint32_t _count,
								 maths::Bounds3f const &_bounds) :
	bounds{ _bounds },
	first_primitive_index{ _first_index }, primitive_count{ _count },
	children{ nullptr, nullptr }, split_axis{ 0 }
{}
BvhAccelerator::BvhNode::BvhNode(uint32_t _axis,
								 BvhNode &_left,
								 BvhNode &_right) :
	bounds{ maths::Union(_left.bounds, _right.bounds) },
	first_primitive_index{ 0 }, primitive_count{ 0 },
	children{ &_left, &_right }, split_axis{ _axis }
{}

} // namespace raytracer
