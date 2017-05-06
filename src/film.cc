#include "film.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "globals.h"
#include "ys_assert.h"
#include "vector.h"

namespace raytracer
{

Film::Film(int32_t _width, int32_t _height, maths::Decimal _side) :
	pixels_(_width * _height, maths::Vec3f{ 0.f }),
	resolution_{ _width, _height },
	aspect_{ resolution_.w / maths::Decimal(resolution_.h) },
	dimensions_{ _side, _side / aspect_ }
{
	YS_ASSERT(_width > 0 && _height > 0);
}

void
Film::SetPixel(maths::Vec3f const &_value, maths::Vec2i32 const &_pos)
{
	YS_ASSERT(_pos.x >= 0 && _pos.y >= 0);
	pixels_[_pos.x + _pos.y * resolution_.w] = _value;
}

void
Film::WriteToFile(std::string const &_path) const
{
	TIMED_SCOPE_FULL(WriteToFile, globals::profiler);

	int blue_count = 0, green_count = 0;
	uint8_t	*buffer = new uint8_t[resolution_.w * resolution_.h * 3];
	
	uint8_t	*head = buffer;
	for (int32_t j = resolution_.h - 1; j >= 0; --j)
	{
		for (int32_t i = 0; i < resolution_.w; ++i)
		{
			maths::Vec3f pixel_value = MapToLimitedRange(pixels_[j * resolution_.w + i]);
			
			YS_ASSERT(pixel_value.r >= 0.f && pixel_value.g >= 0.f && pixel_value.b >= 0.f);
			YS_ASSERT(pixel_value.r <= 1.f && pixel_value.g <= 1.f && pixel_value.b <= 1.f);

			maths::Vector<uint8_t, 3>	rgb_color = (maths::Vector<uint8_t, 3>)(pixel_value * 255._d);
			*head = rgb_color.r;
			*(head + 1) = rgb_color.g;
			*(head + 2) = rgb_color.b;
			head += 3;
		}
	}

	stbi_write_png(_path.c_str(),
				   resolution_.w, resolution_.h, 3, buffer, resolution_.w * 3);

	delete[] buffer;
}

maths::Vec3f
Film::MapToLimitedRange(maths::Vec3f const &_color) const
{
	return maths::Clamp(_color, maths::zero<maths::Decimal>, maths::one<maths::Decimal>);
}

} // raytracer
