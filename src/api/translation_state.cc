#include "api/translation_state.h"

#include "boost/filesystem.hpp"

#include "api/factory_functions.h"
#include "maths/transform.h"


namespace api {


TranslationState::TranslationState():
	workdir_{ boost::filesystem::current_path().string() },
	output_path_{ boost::filesystem::absolute("image.png", workdir_).string() },
	render_context_{},
	scope_depth_{ 1 }, transform_stack_{ maths::Transform{} },
	mem_region_{}, parameters_{new (mem_region_) ParamSet()},
	scene_desc_{}
{}
void
TranslationState::Workdir(std::string const &_absolute_path)
{
	workdir_ = boost::filesystem::absolute(_absolute_path).string();
}
void
TranslationState::Output(std::string const &_relative_path)
{
	output_path_ = boost::filesystem::absolute(_relative_path, workdir_).string();
}
void
TranslationState::Film()
{
	PushObjectDesc_(ObjectIdentifier::kFilm, "");
}
void
TranslationState::Camera()
{
	PushObjectDesc_(ObjectIdentifier::kCamera, "");
}
void
TranslationState::Shape(std::string const &_type)
{
	maths::Transform const &transform =
		render_context_.transform_cache().Lookup(transform_stack_.back());
	param_set().PushTransform("world_transform", transform);
	PushObjectDesc_(ObjectIdentifier::kShape, _type);
}
void
TranslationState::Sampler(std::string const &_type)
{
	PushObjectDesc_(ObjectIdentifier::kSampler, _type);
}
void
TranslationState::Integrator(std::string const &_type)
{
	PushObjectDesc_(ObjectIdentifier::kIntegrator, _type);
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
TranslationState::ScopeBegin()
{
	transform_stack_.push_back(transform_stack_.back());
	scope_depth_++;
}
void
TranslationState::ScopeEnd()
{
	scope_depth_--;
	transform_stack_.pop_back();
	param_set().Clear();
}
void
TranslationState::SceneBegin()
{
	render_context_.Clear();
}
void
TranslationState::SceneEnd()
{
	SceneSetup_();
}


void
TranslationState::SceneSetup_()
{
	{ // NOTE: the integrator has to be created first, see render_context.h:45,50,55
		ObjectDescriptorContainer_t &descriptors = object_desc_vector_(ObjectIdentifier::kIntegrator);
		YS_ASSERT(!descriptors.empty());
		ObjectDescriptor_t const &descriptor = descriptors.back();
		std::string const &type = std::get<0>(descriptor);
		ParamSet const &parameters = *std::get<1>(descriptor);
		MakeIntegratorCallback_t const &callback = LookupIntegratorFunc(type);
		raytracer::Integrator *integrator = callback(render_context_, parameters);
		render_context_.SetIntegrator(integrator);
	}
	{
		ObjectDescriptorContainer_t &descriptors = object_desc_vector_(ObjectIdentifier::kFilm);
		YS_ASSERT(!descriptors.empty());
		ParamSet const &parameters = *std::get<1>(descriptors.back());
		raytracer::Film *film = MakeFilm(render_context_, parameters);
		render_context_.SetFilm(film);
	}
	{
		ObjectDescriptorContainer_t &descriptors = object_desc_vector_(ObjectIdentifier::kCamera);
		YS_ASSERT(!descriptors.empty());
		ParamSet const &parameters = *std::get<1>(descriptors.back());
		raytracer::Camera *camera = MakeCamera(render_context_, parameters);
		render_context_.SetCamera(camera);
	}
	{
		ObjectDescriptorContainer_t &descriptors = object_desc_vector_(ObjectIdentifier::kSampler);
		YS_ASSERT(!descriptors.empty());
		ObjectDescriptor_t const &descriptor = descriptors.back();
		std::string const &type = std::get<0>(descriptor);
		ParamSet const &parameters = *std::get<1>(descriptor);
		MakeSamplerCallback_t const &callback = LookupSamplerFunc(type);
		raytracer::Sampler *sampler = callback(render_context_, parameters);
		render_context_.SetSampler(sampler);
	}
	{
		ObjectDescriptorContainer_t &descriptors = object_desc_vector_(ObjectIdentifier::kShape);
		for (ObjectDescriptorContainer_t::const_iterator odcit = descriptors.cbegin();
			 odcit != descriptors.cend(); ++odcit)
		{
			ObjectDescriptor_t const &descriptor = *odcit;
			std::string const &type = std::get<0>(descriptor);
			ParamSet const &parameters = *std::get<1>(descriptor);
			MakeShapeCallback_t	const &callback = LookupShapeFunc(type);
			std::vector<raytracer::Shape*> shapes = callback(render_context_, parameters);
			for (size_t i = 0; i < shapes.size(); ++i)
			{
				raytracer::GeometryPrimitive *prim =
				new (render_context_.mem_region()) raytracer::GeometryPrimitive(*shapes[i]);
				render_context_.AddPrimitive(prim);
			}
		}
	}
}

void
TranslationState::PushObjectDesc_(ObjectIdentifier const _id, std::string const &_name)
{
	ObjectDescriptorContainer_t &descriptors = object_desc_vector_(_id);
	descriptors.push_back(std::make_tuple(_name, parameters_));
	parameters_ = new (mem_region_) ParamSet();
}

TranslationState::ObjectDescriptorContainer_t &
TranslationState::object_desc_vector_(ObjectIdentifier const _id)
{
	return scene_desc_[static_cast<unsigned>(_id)];
}



} // namespace api
