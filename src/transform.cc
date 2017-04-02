#include "transform.h"

namespace maths
{


bool
Transform::IsIdentity() const
{
	YS_ASSERT(mInv_ == m_);
	return m_ == Mat4x4f::Identity();
}


namespace transform
{

Transform
RotateX(float _theta)
{
	float sin_theta = std::sin(Radians<float>(_theta));
	float cos_theta = std::cos(Radians<float>(_theta));
	Mat4x4f m {	1.f, 0.f, 0.f, 0.f,
				0.f, cos_theta, -sin_theta, 0.f,
				0.f, sin_theta, cos_theta, 0.f,
				0.f, 0.f, 0.f, 1.f };
	return Transform{ m, matrix::Transpose(m) };
}
Transform
RotateY(float _theta)
{
	float sin_theta = std::sin(Radians(_theta));
	float cos_theta = std::cos(Radians(_theta));
	Mat4x4f m {	cos_theta, 0.f, sin_theta, 0.f,
				0.f, 1.f, 0.f, 0.f,
				-sin_theta, 0.f, cos_theta, 0.f,
				0.f, 0.f, 0.f, 1.f };
	return Transform{ m, matrix::Transpose(m) };
}
Transform
RotateZ(float _theta)
{
	float sin_theta = std::sin(Radians(_theta));
	float cos_theta = std::cos(Radians(_theta));
	Mat4x4f m { cos_theta, -sin_theta, 0.f, 0.f,
				sin_theta, cos_theta, 0.f, 0.f,
				0.f, 0.f, 1.f, 0.f,
				0.f, 0.f, 0.f, 1.f };
	return Transform{ m, matrix::Transpose(m) };
}
Transform
Rotate(float _theta, Vec3f const &_axis)
{
	Vec3f a = vector::Normalized(_axis);
	float sin_theta = std::sin(Radians(_theta));
	float cos_theta = std::cos(Radians(_theta));
	Mat4x4f m{};

	m[0][0] = a.x * a.x + (1.f - a.x * a.x) * cos_theta;
	m[0][1] = a.x * a.y * (1.f - cos_theta) - a.z * sin_theta;
	m[0][2] = a.x * a.z * (1.f - cos_theta) + a.y * sin_theta;

	m[1][0] = a.y * a.x * (1.f - cos_theta) + a.z * sin_theta;
	m[1][1] = a.y * a.y + (1.f - a.y * a.y) * cos_theta;
	m[1][2] = a.y * a.z * (1.f - cos_theta) - a.x * sin_theta;

	m[2][0] = a.z * a.x * (1.f - cos_theta) - a.y * sin_theta;
	m[2][1] = a.z * a.y * (1.f - cos_theta) + a.x * sin_theta;
	m[2][2] = a.z * a.z + (1.f - a.z * a.z) * cos_theta;

	return Transform{ m, matrix::Transpose(m) };
}

Transform
LookAt(Vec3f const &_position, Vec3f const &_target, Vec3f const &_up)
{
	Vec3f direction = vector::Normalized(_target - _position);
	Vec3f left = vector::Normalized(vector::Cross(vector::Normalized(_up), direction));
	Vec3f ortho_up = vector::Cross(direction, left);

	Mat4x4f camera_to_world{};
	camera_to_world.SetColumn(0, { left, 0.f });
	camera_to_world.SetColumn(1, { ortho_up, 0.f });
	camera_to_world.SetColumn(2, { direction, 0.f });
	camera_to_world.SetColumn(3, { _position, 1.f });

	return Transform{ matrix::Inverse(camera_to_world), camera_to_world };
}

} // namespace transform
} // namespace maths
