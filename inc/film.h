#ifndef __YS_FILM_HPP__
#define __YS_FILM_HPP__

#include <vector>
#include <cstdint>

#include "vector.h"


namespace raytracer
{


/// Stores a buffer of 3 float vectors.
/// Their values can be anything, as they will be "tonemapped" when written to a file.
/// No tonemapping function has been implemented yet. It's only a clamp.
class Film
{
public:
	Film() = delete;
	Film(int32_t _width, int32_t _height);

	void	SetPixel(maths::Vec3f &&_value, maths::Vec2i32 &&_pos);
	void	WriteToFile(std::string const &_path) const;

	// NOTE: Even though this function could be static, it is intended to define an interface for
	//		 different kinds of films, as they could use different tonemapping functions.
	maths::Vec3f	MapToLimitedRange(maths::Vec3f const &_color) const;

private:
	std::vector<maths::Vec3f>	pixels_;
	maths::Vec2i32				size_;
};


} // raytracer


#endif // __YS_FILM_HPP__
