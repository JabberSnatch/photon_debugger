#pragma once
#ifndef __YS_RESOURCE_CONTEXT_HPP__
#define __YS_RESOURCE_CONTEXT_HPP__

#include <unordered_set>
#include <vector>

#include "api/param_set.h"
#include "api/transform_cache.h"
#include "core/memory_region.h"
#include "core/noncopyable.h"
#include "core/nonmovable.h"


namespace raytracer {
class Shape;
} // namespace raytracer


namespace api {


class ResourceContext final :
	private core::noncopyable,
	private core::nonmovable
{
public:
	enum class ObjectType
	{
		kFilm = 0,
		kCamera,
		kShape,
		kLight,
		kSampler,
		kIntegrator,
		kCount
	};
	template <typename T> static constexpr ObjectType GetType();
public:
	struct ObjectDescriptor
	{
		std::string const	unique_id;
		ObjectType const	type_id;
		ParamSet const		&param_set;
		std::string const	subtype_id;
	};
	using ObjectDescriptorContainer_t = std::vector<ObjectDescriptor const *>;
private:
	struct ObjectInstance
	{
		ObjectInstance(std::string const *const _unique_id, void *const _instance) :
			unique_id{ _unique_id }, instance{ _instance }
		{}
		std::string const	*unique_id;
		void				*instance;
	};
	using ObjectInstanceContainer_t = std::vector<ObjectInstance>;
private:
	using UsedShapePtrContainer_t = std::unordered_set<raytracer::Shape const*>;
public:
	explicit ResourceContext(std::string const &_workdir);
	bool IsUniqueIdFree(std::string const &_unique_id) const;
	void PushDescriptor(std::string const &_unique_id,
						ObjectType _type,
						ParamSet const &_param_set,
						std::string const &_subtype_id = "");
	ObjectDescriptor const &GetAnyDescOfType(ObjectType const _type) const;
	ObjectDescriptorContainer_t GetAllDescsOfType(ObjectType const _type) const;
	template <typename T> T& Fetch(std::string const &_unique_id);
	void SetWorkdir(std::string const &_workdir);
	std::string const	&workdir() const;
	core::MemoryRegion	&mem_region();
	TransformCache		&transform_cache();
public:
	void FlagLightShape(raytracer::Shape const &_shape);
	bool IsShapeLight(raytracer::Shape const &_shape) const;
private:
	template <typename T> T* MakeObject_(ObjectDescriptor const &_object_desc);
private:
	std::string						workdir_{ "" };
	core::MemoryRegion				mem_region_{};
	TransformCache					transform_cache_{};
	ObjectDescriptorContainer_t		object_descriptors_{};
	ObjectInstanceContainer_t		object_instances_{};
	UsedShapePtrContainer_t			light_shapes_{};
};

} // namespace api


#endif // __YS_RESOURCE_CONTEXT_HPP__
