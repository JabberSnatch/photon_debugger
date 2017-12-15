#include "api/resource_context.h"

#include <algorithm>
#include <string>
#include <sstream>

#include <boost/numeric/conversion/cast.hpp>

#include "api/factory_functions.h"
#include "raytracer/light.h"

namespace api {


template <typename T>
constexpr ResourceContext::ObjectType
ResourceContext::GetType()
{
	static_assert(false, "Unimplemented type");
	return ObjectType::kCount;
}
template <>
constexpr ResourceContext::ObjectType
ResourceContext::GetType<raytracer::Film>()
{
	return ObjectType::kFilm;
}
template <>
constexpr ResourceContext::ObjectType
ResourceContext::GetType<raytracer::Camera>()
{
	return ObjectType::kCamera;
}
template <>
constexpr ResourceContext::ObjectType
ResourceContext::GetType<raytracer::Shape>()
{
	return ObjectType::kShape;
}
template <>
constexpr ResourceContext::ObjectType
ResourceContext::GetType<raytracer::Light>()
{
	return ObjectType::kLight;
}
template <>
constexpr ResourceContext::ObjectType
ResourceContext::GetType<raytracer::Sampler>()
{
	return ObjectType::kSampler;
}
template <>
constexpr ResourceContext::ObjectType
ResourceContext::GetType<raytracer::Integrator>()
{
	return ObjectType::kIntegrator;
}


ResourceContext::ResourceContext()
{
	std::fill(descriptor_counts_.begin(), descriptor_counts_.end(), 0u);
}


bool
ResourceContext::IsUniqueIdFree(std::string const &_unique_id) const
{
	ObjectDescriptorContainer_t::const_iterator odcit = std::find_if(
		object_descriptors_.cbegin(), object_descriptors_.cend(),
		[&_unique_id](ObjectDescriptor const *_object_desc)
	{
		return _object_desc->unique_id == _unique_id;
	});
	return odcit == object_descriptors_.cend();
}


void
ResourceContext::PushDescriptor(std::string const &_unique_id,
								ObjectType _type,
								ParamSet const &_param_set,
								std::string const &_subtype_id)
{
	if (IsUniqueIdFree(_unique_id))
	{
		object_descriptors_.emplace_back(new (mem_region_)
			ObjectDescriptor{ _unique_id, _type, _param_set, _subtype_id });
		++(descriptor_counts_[static_cast<uint32_t>(_type)]);
	}
	else
	{
		LOG_WARNING(tools::kChannelGeneral, "Object ID " + _unique_id + " is referenced by two separate descriptors.");
	}
}


ResourceContext::ObjectDescriptor const &
ResourceContext::GetAnyDescOfType(ObjectType const _type) const
{
	ObjectDescriptorContainer_t::const_iterator odcit = std::find_if(
		object_descriptors_.cbegin(), object_descriptors_.cend(),
		[&_type](ObjectDescriptor const *_object_desc)
	{
		return _object_desc->type_id == _type;
	});
	if (odcit == object_descriptors_.cend())
	{
		LOG_ERROR(tools::kChannelGeneral, "No descriptor of type " + std::to_string(static_cast<int>(_type)) + " exists. Shutdown.");
	}
	return **odcit;
}


ResourceContext::ObjectDescriptorContainer_t
ResourceContext::GetAllDescsOfType(ObjectType const _type) const
{
	ObjectDescriptorContainer_t result{};
	std::copy_if(object_descriptors_.cbegin(), object_descriptors_.cend(),
				 std::back_inserter(result),
				 [_type](ObjectDescriptor const *_object_desc)
	{
		return _object_desc->type_id == _type;
	});
	YS_ASSERT(result.size() ==
			  boost::numeric_cast<size_t>(descriptor_counts_[static_cast<uint32_t>(_type)]));
	return result;
}


template <typename T>
T&
ResourceContext::Fetch(std::string const &_unique_id)
{
	T* result = nullptr;
	ObjectInstanceContainer_t::const_iterator oicit = std::find_if(
		object_instances_.cbegin(), object_instances_.cend(),
		[&_unique_id](ObjectInstance const &_object_instance)
	{
		return *(_object_instance.unique_id) == _unique_id;
	});
	if (oicit == object_instances_.cend())
	{
		ObjectDescriptorContainer_t::const_iterator odcit = std::find_if(
			object_descriptors_.cbegin(), object_descriptors_.cend(),
			[&_unique_id](ObjectDescriptor const *_object_desc)
		{
			return _object_desc->unique_id == _unique_id;
		});
		if (odcit != object_descriptors_.cend())
		{
			ObjectDescriptor const &object_desc = **odcit;
			if (GetType<T>() == object_desc.type_id)
			{
				result = MakeObject_<T>(object_desc);
				if (result)
				{
					object_instances_.emplace_back(&(object_desc.unique_id), result);
				}
				else
				{
					LOG_ERROR(tools::kChannelGeneral, "Unexpected error, no object created.");
				}
			}
			else
			{
				LOG_ERROR(tools::kChannelGeneral, "Fetch expected an object of type " + std::to_string(static_cast<int> (GetType<T>())) + " but got " + std::to_string(static_cast<int>(object_desc.type_id)) + " instead.	 Object ID : " + object_desc.unique_id);
			}
		}
		else
		{
			LOG_ERROR(tools::kChannelGeneral, "Object ID " + _unique_id + " doesn't exist.");
		}
	}
	else
	{
		result = reinterpret_cast<T*>(oicit->instance);
	}
	return *result;
}
template
raytracer::Film&
ResourceContext::Fetch<raytracer::Film>(std::string const &_unique_id);
template
raytracer::Camera&
ResourceContext::Fetch<raytracer::Camera>(std::string const &_unique_id);
template
raytracer::Shape&
ResourceContext::Fetch<raytracer::Shape>(std::string const &_unique_id);
template
raytracer::Sampler&
ResourceContext::Fetch<raytracer::Sampler>(std::string const &_unique_id);
template
raytracer::Integrator&
ResourceContext::Fetch<raytracer::Integrator>(std::string const &_unique_id);


void
ResourceContext::SetWorkdir(std::string const &_workdir)
{
	workdir_ = _workdir;
}


std::string
ResourceContext::MakeUniqueID(ObjectType const _type) const
{
	std::stringstream result_stream;
	if (_type == ObjectType::kFilm)
	{
		result_stream << "film";
	}
	else if (_type == ObjectType::kCamera)
	{
		result_stream << "camera";
	}
	else if (_type == ObjectType::kShape)
	{
		result_stream << "shape";
	}
	else if (_type == ObjectType::kLight)
	{
		result_stream << "light";
	}
	else if (_type == ObjectType::kSampler)
	{
		result_stream << "sampler";
	}
	else if (_type == ObjectType::kIntegrator)
	{
		result_stream << "integrator";
	}
	else
	{
		YS_ASSERT(false);
		LOG_ERROR(tools::kChannelGeneral, "Unhandled ObjectType");
	}
	uint32_t const existing_count = descriptor_counts_[static_cast<uint32_t>(_type)];
	result_stream << existing_count;
	return result_stream.str();
}


std::string const &
ResourceContext::workdir() const
{
	return workdir_;
}


core::MemoryRegion &
ResourceContext::mem_region()
{
	return mem_region_;
}


TransformCache &
ResourceContext::transform_cache()
{
	return transform_cache_;
}



template <typename T>
T*
ResourceContext::MakeObject_(ObjectDescriptor const &_object_desc)
{
	static_assert(false, "Unimplemented type");
	return nullptr;
}
template <>
raytracer::Film*
ResourceContext::MakeObject_(ObjectDescriptor const &_object_desc)
{
	return MakeFilm(*this, _object_desc.param_set);
}
template <>
raytracer::Camera*
ResourceContext::MakeObject_(ObjectDescriptor const &_object_desc)
{
	return MakeCamera(*this, _object_desc.param_set);
}
template <>
raytracer::Shape*
ResourceContext::MakeObject_(ObjectDescriptor const &_object_desc)
{
	MakeShapeCallback_t const &callback = LookupShapeFunc(_object_desc.subtype_id);
	return callback(*this, _object_desc.param_set);
}
template <>
raytracer::Sampler*
ResourceContext::MakeObject_(ObjectDescriptor const &_object_desc)
{
	MakeSamplerCallback_t const &callback = LookupSamplerFunc(_object_desc.subtype_id);
	return callback(*this, _object_desc.param_set);
}
template <>
raytracer::Integrator*
ResourceContext::MakeObject_(ObjectDescriptor const &_object_desc)
{
	MakeIntegratorCallback_t const &callback = LookupIntegratorFunc(_object_desc.subtype_id);
	return callback(*this, _object_desc.param_set);
}


} // namespace api
