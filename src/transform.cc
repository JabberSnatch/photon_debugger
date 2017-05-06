#include "transform.h"

#include "vector.h"
#include "surface_interaction.h"


namespace maths
{


bool
Transform::IsIdentity() const
{
	YS_ASSERT(mInv_ == m_);
	return m_ == Mat4x4f::Identity();
}

bool
Transform::SwapsHandedness() const
{
	Decimal determinant =
		m_[0][0] * (m_[1][1] * m_[2][2] - m_[1][2] * m_[2][1]) -
		m_[0][1] * (m_[1][0] * m_[2][2] - m_[1][2] * m_[2][0]) +
		m_[0][2] * (m_[1][0] * m_[2][1] - m_[1][1] * m_[2][0]);
	return determinant < 0._d;
}

raytracer::SurfaceInteraction
Transform::operator()(raytracer::SurfaceInteraction const &_v, OpDirection _dir) const
{
	raytracer::SurfaceInteraction result;
	result.position = (*this)(_v.position, _dir);
	result.time = _v.time;
	result.wo = (*this)(_v.wo, _dir);
	result.shape = _v.shape;
	result.uv = _v.uv;
	result.geometry.normal = (*this)(_v.geometry.normal, _dir);
	result.geometry.dpdu = (*this)(_v.geometry.dpdu, _dir);
	result.geometry.dpdv = (*this)(_v.geometry.dpdv, _dir);
	result.geometry.dndu = (*this)(_v.geometry.dndu, _dir);
	result.geometry.dndv = (*this)(_v.geometry.dndv, _dir);
	result.shading.normal = (*this)(_v.shading.normal, _dir);
	result.shading.dpdu = (*this)(_v.shading.dpdu, _dir);
	result.shading.dpdv = (*this)(_v.shading.dpdv, _dir);
	result.shading.dndu = (*this)(_v.shading.dndu, _dir);
	result.shading.dndv = (*this)(_v.shading.dndv, _dir);
	return result;
}


Transform
RotateX(Decimal _theta)
{
	Decimal sin_theta = std::sin(Radians(_theta));
	Decimal cos_theta = std::cos(Radians(_theta));
	Mat4x4f m {	1.f, 0.f, 0.f, 0.f,
				0.f, cos_theta, -sin_theta, 0.f,
				0.f, sin_theta, cos_theta, 0.f,
				0.f, 0.f, 0.f, 1.f };
	return Transform{ m, Transpose(m) };
}
Transform
RotateY(Decimal _theta)
{
	Decimal sin_theta = std::sin(Radians(_theta));
	Decimal cos_theta = std::cos(Radians(_theta));
	Mat4x4f m {	cos_theta, 0.f, sin_theta, 0.f,
				0.f, 1.f, 0.f, 0.f,
				-sin_theta, 0.f, cos_theta, 0.f,
				0.f, 0.f, 0.f, 1.f };
	return Transform{ m, Transpose(m) };
}
Transform
RotateZ(Decimal _theta)
{
	Decimal sin_theta = std::sin(Radians(_theta));
	Decimal cos_theta = std::cos(Radians(_theta));
	Mat4x4f m { cos_theta, -sin_theta, 0.f, 0.f,
				sin_theta, cos_theta, 0.f, 0.f,
				0.f, 0.f, 1.f, 0.f,
				0.f, 0.f, 0.f, 1.f };
	return Transform{ m, Transpose(m) };
}
Transform
Rotate(Decimal _theta, Vec3f const &_axis)
{
	Vec3f a = Normalized(_axis);
	Decimal sin_theta = std::sin(Radians(_theta));
	Decimal cos_theta = std::cos(Radians(_theta));
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

	return Transform{ m, Transpose(m) };
}

Transform
LookAt(Vec3f const &_position, Vec3f const &_target, Vec3f const &_up)
{
	Vec3f direction = Normalized(_target - _position);
	Vec3f left = Normalized(Cross(Normalized(_up), direction));
	Vec3f ortho_up = Cross(direction, left);

	Mat4x4f camera_to_world{};
	camera_to_world.SetColumn(0, { left, 0.f });
	camera_to_world.SetColumn(1, { ortho_up, 0.f });
	camera_to_world.SetColumn(2, { direction, 0.f });
	camera_to_world.SetColumn(3, { _position, 1.f });

	return Transform{ Inverse(camera_to_world), camera_to_world };
}

} // namespace maths
