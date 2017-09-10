#include "api/translation_state.h"

#include "boost/filesystem.hpp"

#include "api/factory_functions.h"
#include "maths/transform.h"


namespace api {


TranslationState::TranslationState() :
	workdir_{ boost::filesystem::current_path().string() },
	output_path_{ boost::filesystem::absolute("image.png", workdir_).string() },
	render_context_{}, parameters_{},
	scope_depth_{ 1 }, transform_stack_{ maths::Transform{} }
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
	raytracer::Film *film = MakeFilm(render_context_, parameters_);
	render_context_.SetFilm(film);
}
void
TranslationState::Camera()
{
	raytracer::Camera *camera = MakeCamera(render_context_, parameters_);
	render_context_.SetCamera(camera);
}
void
TranslationState::Shape(std::string const &_type)
{
	maths::Transform const &transform =
		render_context_.transform_cache().Lookup(transform_stack_.back());
	parameters_.PushTransform("world_transform", transform);
	MakeShapeCallback_t	const &callback = LookupShapeFunc(_type);
	std::vector<raytracer::Shape*> shapes = callback(render_context_, parameters_);
	for (size_t i = 0; i < shapes.size(); ++i)
	{
		raytracer::GeometryPrimitive *prim =
			new (render_context_.mem_region()) raytracer::GeometryPrimitive(*shapes[i]);
		render_context_.AddPrimitive(prim);
	}
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
	parameters_.Clear();
}
void
TranslationState::SceneBegin()
{
	render_context_.Clear();
}
void
TranslationState::SceneEnd()
{
	if (render_context_.GoodForRender())
		render_context_.RenderAndWrite(output_path_);
}


} // namespace api
