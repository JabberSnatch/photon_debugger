#include "raytracer/film.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "globals.h"
#include "maths/vector.h"
#include "core/profiler.h"
#include "raytracer/render_context.h"
#include "api/param_set.h"

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
	using BmpColor = maths::Vector<uint8_t, 3>;
	TIMED_SCOPE(WriteToFile);

	int blue_count = 0, green_count = 0;
	uint8_t	*buffer = new uint8_t[resolution_.w * resolution_.h * 3];
	
	uint8_t	*head = buffer;
	for (int32_t j = resolution_.h - 1; j >= 0; --j)
	{
		for (int32_t i = 0; i < resolution_.w; ++i)
		{
			int32_t	x{ i }, y{ j };
			// NOTE: Depending on our camera implementation, we might want to unflip a picture
			//		 before writing it to a file. We assume that both axes are flipped
			//		 which might not always be the case.
			if (image_is_flipped)
			{
				x = resolution_.w - 1 - x;
				y = resolution_.h - 1 - y;
			}

			maths::Vec3f pixel_value = MapToLimitedRange(pixels_[y * resolution_.w + x]);
			
			YS_ASSERT(pixel_value.r >= 0.f && pixel_value.g >= 0.f && pixel_value.b >= 0.f);
			YS_ASSERT(pixel_value.r <= 1.f && pixel_value.g <= 1.f && pixel_value.b <= 1.f);

			BmpColor	rgb_color = (BmpColor)(pixel_value * 255._d);
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
