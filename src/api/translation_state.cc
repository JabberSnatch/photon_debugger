#include "api/translation_state.h"

#include <sstream>

#include "boost/filesystem.hpp"

#include "api/factory_functions.h"
#include "api/render_context.h"
#include "maths/transform.h"


namespace api {


TranslationState::TranslationState():
	workdir_{ boost::filesystem::current_path().string() },
	output_path_{ boost::filesystem::absolute("image.png", workdir_).string() },
	resource_context_{}, render_context_{},
	scope_depth_{ 1 }, transform_stack_{ maths::Transform{} },
	cached_object_id_{ "" }, parameters_{ nullptr }
{
	ResetResourceCounters();
}
void
TranslationState::SceneBegin()
{
	render_context_.Clear();
	parameters_ = nullptr;
	YS_ASSERT(scope_depth_ == 1);
	YS_ASSERT(transform_stack_.size() == 1 && transform_stack_.back() == maths::Transform{});
	YS_ASSERT(cached_object_id_.empty());
}
void
TranslationState::SceneEnd()
{
	SceneSetup_();
}
void
TranslationState::ScopeBegin()
{
	parameters_ = new (resource_context_.mem_region()) ParamSet();
	transform_stack_.push_back(transform_stack_.back());
	scope_depth_++;
}
void
TranslationState::ScopeEnd()
{
	scope_depth_--;
	transform_stack_.pop_back();
	cached_object_id_.clear();
}
void
TranslationState::Workdir(std::string const &_absolute_path)
{
	workdir_ = boost::filesystem::absolute(_absolute_path).string();
	render_context_.workdir = workdir_;
}
void
TranslationState::Output(std::string const &_relative_path)
{
	output_path_ = boost::filesystem::absolute(_relative_path, workdir_).string();
}
void
TranslationState::ObjectId(std::string const &_object_id)
{
	if (cached_object_id_.empty())
	{
		cached_object_id_ = _object_id;
	}
	else
	{
		// Ideally, the parser should be in charge of catching this error
		LOG_ERROR(tools::kChannelParsing, "An object has more than one ID field, ignored the latest");
	}
}
void
TranslationState::Film()
{
	PushObjectDesc_(ResourceContext::ObjectType::kFilm, "");
}
void
TranslationState::Camera()
{
	PushObjectDesc_(ResourceContext::ObjectType::kCamera, "");
}
void
TranslationState::Shape(std::string const &_type)
{
	maths::Transform const &transform =
		resource_context_.transform_cache().Lookup(transform_stack_.back());
	param_set().PushTransform("world_transform", transform);
	PushObjectDesc_(ResourceContext::ObjectType::kShape, _type);
}
void
TranslationState::Light(std::string const &_type)
{
	PushObjectDesc_(ResourceContext::ObjectType::kLight, _type);
}
void
TranslationState::Sampler(std::string const &_type)
{
	PushObjectDesc_(ResourceContext::ObjectType::kSampler, _type);
}
void
TranslationState::Integrator(std::string const &_type)
{
	PushObjectDesc_(ResourceContext::ObjectType::kIntegrator, _type);
}
void
TranslationState::Identity()
{
	transform_stack_.back() = maths::Transform{};
}
void
TranslationState::Translate(maths::Vec3f const &_t)
{
	transform_stack_.back() = transform_stack_.back() * maths::Translate(_t);
}
void
TranslationState::Rotate(maths::Decimal _angle, maths::Vec3f const &_axis)
{
	transform_stack_.back() = transform_stack_.back() * maths::Rotate(_angle, _axis);
}
void
TranslationState::Scale(maths::Decimal _x, maths::Decimal _y, maths::Decimal _z)
{
	transform_stack_.back() = transform_stack_.back() * maths::Scale(_x, _y, _z);
}


void
TranslationState::SceneSetup_()
{
	ResourceContext::ObjectDescriptor const &integrator_desc =
		resource_context_.GetAnyDescOfType(ResourceContext::ObjectType::kIntegrator);
	ResourceContext::ObjectDescriptor const &camera_desc =
		resource_context_.GetAnyDescOfType(ResourceContext::ObjectType::kCamera);
	ResourceContext::ObjectDescriptor const &sampler_desc =
		resource_context_.GetAnyDescOfType(ResourceContext::ObjectType::kSampler);
	ResourceContext::ObjectDescriptorContainer_t shape_descs =
		resource_context_.GetAllDescsOfType(ResourceContext::ObjectType::kShape);
	ResourceContext::ObjectDescriptorContainer_t light_descs =
		resource_context_.GetAllDescsOfType(ResourceContext::ObjectType::kLight);
	//
	raytracer::Integrator	&integrator =
		resource_context_.Fetch<raytracer::Integrator>(integrator_desc.unique_id);
	raytracer::Camera		&camera =
		resource_context_.Fetch<raytracer::Camera>(camera_desc.unique_id);
	raytracer::Sampler		&sampler =
		resource_context_.Fetch<raytracer::Sampler>(sampler_desc.unique_id);
	//
	RenderContext::LightContainer_t lights{};
	lights.reserve(light_descs.size());
	std::transform(light_descs.cbegin(), light_descs.cend(), std::back_inserter(lights),
				   [this](ResourceContext::ObjectDescriptor const *_object_desc) {
					   raytracer::Light const &light =
						   resource_context_.Fetch<raytracer::Light>(_object_desc->unique_id);
					   return &light;
				   });
	//
	ResourceContext::ObjectDescriptorContainer_t::iterator valid_shapes_end =
		std::remove_if(shape_descs.begin(), shape_descs.end(),
					   [this](ResourceContext::ObjectDescriptor const *_object_desc) {
						   raytracer::Shape const &shape =
							   resource_context_.Fetch<raytracer::Shape>(_object_desc->unique_id);
						   return resource_context_.IsShapeLight(shape);
					   });
	shape_descs.erase(valid_shapes_end, shape_descs.end());
	RenderContext::PrimitiveContainer_t primitives{};
	primitives.reserve(shape_descs.size());
	std::transform(shape_descs.cbegin(), shape_descs.cend(), std::back_inserter(primitives),
				   [this](ResourceContext::ObjectDescriptor const *_object_desc) {
					   raytracer::Shape const &shape =
						   resource_context_.Fetch<raytracer::Shape>(_object_desc->unique_id);
					   return new (resource_context_.mem_region()) raytracer::GeometryPrimitive(shape);
				   });
	//
	constexpr uint32_t kBvhNodeMaxSize = 20;
	if (primitives.size() > kBvhNodeMaxSize)
	{
		LOG_INFO(tools::kChannelGeneral, "Primitive count exceeded threshold, building BVH");
		raytracer::Primitive *const bvh =
			new (resource_context_.mem_region()) raytracer::BvhAccelerator(primitives, 
																		   kBvhNodeMaxSize);
		primitives.clear();
		primitives.emplace_back(bvh);
	}
	else
	{
	}
	//
	render_context_ = api::RenderContext(camera, sampler, integrator, primitives, lights);
}

void
TranslationState::PushObjectDesc_(ResourceContext::ObjectType const _type,
								  std::string const &_subtype_id)
{
	std::string const unique_id = cached_object_id_.empty() ?
		MakeUniqueResourceID_(_type) :
		cached_object_id_;
	resource_context_.PushDescriptor(unique_id, _type, *parameters_, _subtype_id);
	++(resource_counters_[static_cast<uint32_t>(_type)]);
}


void
TranslationState::ResetResourceCounters()
{
	std::fill(resource_counters_.begin(), resource_counters_.end(), 0u);
}


std::string
TranslationState::MakeUniqueResourceID_(ResourceContext::ObjectType const _type) const
{
	std::stringstream result_stream;
	if (_type == ResourceContext::ObjectType::kFilm)
	{
		result_stream << "film";
	}
	else if (_type == ResourceContext::ObjectType::kCamera)
	{
		result_stream << "camera";
	}
	else if (_type == ResourceContext::ObjectType::kShape)
	{
		result_stream << "shape";
	}
	else if (_type == ResourceContext::ObjectType::kLight)
	{
		result_stream << "light";
	}
	else if (_type == ResourceContext::ObjectType::kSampler)
	{
		result_stream << "sampler";
	}
	else if (_type == ResourceContext::ObjectType::kIntegrator)
	{
		result_stream << "integrator";
	}
	else
	{
		YS_ASSERT(false);
		LOG_ERROR(tools::kChannelGeneral, "Unhandled ObjectType");
	}
	uint32_t const existing_count = resource_counters_[static_cast<uint32_t>(_type)];
	result_stream << existing_count;
	return result_stream.str();
}


} // namespace api
