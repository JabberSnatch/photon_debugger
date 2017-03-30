#ifndef __YS_VECTOR_HPP__
#define __YS_VECTOR_HPP__


namespace maths
{


struct Vec3f final
{
	constexpr Vec3f() :
		Vec3f{ zero() }
	{}
	constexpr Vec3f(float _e0, float _e1, float _e2) :
		e{ _e0, _e1, _e2 }
	{}

	static constexpr Vec3f zero() { return Vec3f{ 0.f, 0.f, 0.f }; }
	static constexpr Vec3f one() { return Vec3f{ 1.f, 1.f, 1.f }; }

	constexpr const Vec3f& operator+() const { return *this; }
	constexpr Vec3f operator-() const { return Vec3f{ -x, -y, -z }; }
	constexpr float operator[](int _i) const { return e[_i]; }
	constexpr float& operator[](int _i) { return e[_i]; }

	constexpr Vec3f& operator+=(Vec3f const &_rhs);
	constexpr Vec3f& operator-=(Vec3f const &_rhs);
	constexpr Vec3f& operator*=(Vec3f const &_rhs);
	constexpr Vec3f& operator/=(Vec3f const &_rhs);
	constexpr Vec3f& operator*=(float _rhs);
	constexpr Vec3f& operator/=(float _rhs);

	constexpr Vec3f normalized() const;
	constexpr float sqr_length() const;
	inline float length() const;

	union
	{
		float	e[3];
		struct
		{
			float x;
			float y;
			float z;
		};
		struct
		{
			float r;
			float g;
			float b;
		};
	};
};

constexpr float Dot(Vec3f const &_lhs, Vec3f const &_rhs);
constexpr Vec3f Cross(Vec3f const &_lhs, Vec3f const &_rhs);
constexpr Vec3f Reflect(Vec3f const &_v, Vec3f const &_n);

constexpr Vec3f operator+(Vec3f const &_lhs, Vec3f const &_rhs);
constexpr Vec3f operator-(Vec3f const &_lhs, Vec3f const &_rhs);
constexpr Vec3f operator*(Vec3f const &_lhs, Vec3f const &_rhs);
constexpr Vec3f operator/(Vec3f const &_lhs, Vec3f const &_rhs);
constexpr Vec3f operator*(float _s, Vec3f const &_v);
constexpr Vec3f operator*(Vec3f const &_v, float _s);
constexpr Vec3f operator/(Vec3f const &_v, float _s);


} // maths


#include "vector.inl"

#endif // __YS_VECTOR_HPP__
