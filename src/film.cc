#include "film.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "profiler.h"
#include "ys_assert.h"

namespace raytracer
{

Film::Film(int32_t _width, int32_t _height) :
	pixels_(_width * _height, maths::Vec3f::zero()),
	size_{ _width, _height }
{
	YS_ASSERT(_width > 0 && _height > 0);
}

void
Film::SetPixel(maths::Vec3f &&_value, maths::Vec2i32 &&_pos)
{
	YS_ASSERT(_pos.x >= 0 && _pos.y >= 0);
	pixels_[_pos.x + _pos.y * size_.w] = _value;
}

void
Film::WriteToFile(std::string const &_path) const
{
	//TIMED_SCOPE(WriteToFile, gProfiler);
	tools::TimeProbe	WriteToFile_probe{ gProfiler.GetTimer("WriteToFile_probe") };


	uint8_t	*buffer = new uint8_t[size_.w * size_.h * 3];
	
	uint8_t	*head = buffer;
	for (int32_t j = size_.h - 1; j >= 0; --j)
	{
		for (int32_t i = 0; i < size_.w; ++i)
		{
			maths::Vec3f pixel_value = MapToLimitedRange(pixels_[j * size_.w + i]);
			
			YS_ASSERT(pixel_value.r >= 0.f && pixel_value.g >= 0.f && pixel_value.b >= 0.f);
			YS_ASSERT(pixel_value.r < 1.f && pixel_value.g < 1.f && pixel_value.b < 1.f);

			*head++ = uint8_t(255.99f * pixel_value.r);
			*head++ = uint8_t(255.99f * pixel_value.g);
			*head++ = uint8_t(255.99f * pixel_value.b);
		}
	}

	stbi_write_png(_path.c_str(),
				   size_.w, size_.h, 3, buffer, size_.h * 3);

	delete[] buffer;
}

maths::Vec3f
Film::MapToLimitedRange(maths::Vec3f const &_color) const
{
	return maths::Clamp(_color, 0.f, maths::almost_onef());
}

} // raytracer
