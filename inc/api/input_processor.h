#pragma once
#ifndef __YS_INPUT_PARSER_HPP__
#define __YS_INPUT_PARSER_HPP__

#include <string>
#include <map>
#include <vector>
#include <functional>

#include "raytracer/raytracer.h"

#include "core/memory_region.h"
#include "api/param_set.h"
#include "maths/transform.h"
#include "raytracer/render_context.h"

namespace api {

// Camera, Film, Shape
// Scene : list of shapes
// Shape -> Transform + flip_normals + individual properties
// For instance : Sphere{Transform, flip_normals, radius, min_z, max_z, max_phi}
// Each parameter can have default value.
// Specifying that we just want a sphere would be fine.
// Values should be assigned explicitly
//
// scope_id [ value ] et ref[ id ] last
//
// Film scope_id [ default_film ] { resolution [ 1920 1080 ] }
// Camera {
// position float [ 0 -5 5 ] lookat float [ 0 1 .5 ]
// up float [ 0 0 1 ]
// fov float [ 60 ]
// film ref[ default_film ]
// } 
// Translate 10 0 0
// Rotate 30 0 0
// Shape sphere {
// radius float [ 5 ]
// min_z float [ -10 ]
// max_z float [ 10 ]
// max_phi float [ 360 ]
// flip_normals bool [ false ]
// Translate 5 0 0
// Rotate 25 0 25
// Scale 3 3 1
// }
// Shape sphere {
// radius float [ 5 ]
// min_z float [ -3 ]
// max_z float [ 4 ]
// max_phi float [ 360 ]
// }


bool	ProcessInputFile(std::string const &_path);



enum TokenId
{
	kNone,

	kOutput,
	kFilm,
	kCamera,
	kShape,
	kScopeBegin,
	kScopeEnd,
	kParamBegin,
	kParamEnd,
	kType,
	kTransformIdentity,
	kTranslate,
	kRotate,
	kScale,
	kString,
	kNumber,

	kValueGroup,
	kDefinitionGroup,
	kParamGroup,
	kPropertiesGroup,
	kAttributeGroup,
	kTranslateGroup,
	kRotateGroup,
	kScaleGroup,
	kFilmGroup,
	kCameraGroup,
	kShapeGroup,
	kOutputGroup,
	kSceneGroup,

	kEnd,
	kDefault,
};


struct Token
{
	TokenId		id;
	std::string	text;
};


class TransformCache final
{
public:
	TransformCache() = default;
	~TransformCache() = default;
	TransformCache(TransformCache const &) = delete;
	TransformCache(TransformCache &&) = delete;
	TransformCache &operator=(TransformCache const &) = delete;
	TransformCache &operator=(TransformCache &&) = delete;

	maths::Transform const &Lookup(maths::Transform const &_t);

private:
	std::map<maths::Transform, maths::Transform *>	lookup_table_;
	core::MemoryRegion<>	mem_region_;
};


// NOTE: Implementation doesn't enforce syntax compliance, be very careful about using it out of the
//		 input processor.
class TranslationState final
{
public:
	using MakeShapeCallback_t = std::function <
		std::vector<raytracer::Shape*>(raytracer::RenderContext &_context, 
									   maths::Transform const & _t, bool _flip_normals,
									   ParamSet const & _params)
	>;

	TranslationState();
	~TranslationState() = default;
	TranslationState(TranslationState const &) = delete;
	TranslationState(TranslationState&&) = delete;
	TranslationState &operator=(TranslationState const &) = delete;
	TranslationState &operator=(TranslationState &&) = delete;

	void			Workdir(std::string const &_absolute_path);
	void			Output(std::string const &_relative_path);
	void			Film();
	void			Camera();
	void			Shape(std::string const &_type);
	void			Identity();
	void			Translate(maths::Vec3f const &_t);
	void			Rotate(maths::Decimal _angle, maths::Vec3f const &_axis);
	void			Scale(maths::Decimal _x, maths::Decimal _y, maths::Decimal _z);
	void			ScopeBegin();
	void			ScopeEnd();
	void			SceneBegin();
	void			SceneEnd();

	ParamSet		&param_set() { return parameters_; }

	static bool		LookupShapeFunc(std::string const &_id, MakeShapeCallback_t &_func);

private:
	std::string					workdir_;
	std::string					output_path_;

	raytracer::RenderContext	render_context_;
	ParamSet					parameters_;
	TransformCache				transform_cache_;

	uint32_t							scope_depth_;
	std::vector<maths::Transform>		transform_stack_;

	static std::unordered_map<std::string, MakeShapeCallback_t> const	shape_callbacks_;
};

} // namespace api


#endif // __YS_INPUT_PARSER_HPP__
