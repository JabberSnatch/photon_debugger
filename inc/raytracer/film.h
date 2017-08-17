#ifndef __YS_FILM_HPP__
#define __YS_FILM_HPP__

#include <vector>

#include "raytracer/raytracer.h"
#include "api/api.h"
#include "maths/maths.h"
#include "maths/vector.h"


namespace raytracer
{


/// Stores a buffer of 3 float vectors.
/// Their values can be anything, as they will be "tonemapped" when written to a file.
/// No tonemapping function has been implemented yet. It's only a clamp.
class Film
{
public:
	Film() = delete;
	Film(int32_t _width, int32_t _height, maths::Decimal _side);

	void	SetPixel(maths::Vec3f const &_value, maths::Vec2i32 const &_pos);
	void	WriteToFile(std::string const &_path) const;

	maths::Vec2i32 const	&resolution() const { return resolution_; }
	maths::Vec2f const		&dimensions() const { return dimensions_; }
	maths::Decimal			aspect() const { return aspect_; }

	// NOTE: Even though this function could be static, it is intended to define an interface for
	//		 different kinds of films, as they could use different tonemapping functions.
	maths::Vec3f	MapToLimitedRange(maths::Vec3f const &_color) const;

	bool						image_is_flipped = false;

private:
	std::vector<maths::Vec3f>	pixels_;
	maths::Vec2i32				resolution_;
	maths::Decimal				aspect_;
	maths::Vec2f				dimensions_;
};


Film *MakeFilm(RenderContext &_context, api::ParamSet const &_params);

} // raytracer


#endif // __YS_FILM_HPP__
