#include "raytracer/shapes/triangle.h"

#include "common_macros.h"
#include "globals.h"
#include "core/profiler.h"
#include "core/logger.h"
#include "maths/bounds.h"
#include "maths/transform.h"
#include "maths/vector.h"
#include "raytracer/triangle_mesh_data.h"
#include "raytracer/surface_interaction.h"

namespace raytracer {

Triangle::Triangle(maths::Transform const &_world_transform, bool _flip_normals,
				   TriangleMeshRawData const &_mesh_data, int32_t _face_index) :
	Shape(_world_transform, _flip_normals),
	mesh_data_{ _mesh_data }
{
	YS_ASSERT(_face_index < mesh_data_.triangle_count);
	vertex_index_ = &mesh_data_.indices[TriangleMeshRawData::IndexOffset(_face_index)];
}

bool
Triangle::Intersect(maths::Ray const &_ray,
					maths::Decimal &_tHit,
					SurfaceInteraction &_hit_info) const
{
	TIMED_SCOPE(Triangle_Intersect);

	maths::Point3f const	&v0 = mesh_data_.vertices[vertex_index_[0]];
	maths::Point3f const	&v1 = mesh_data_.vertices[vertex_index_[1]];
	maths::Point3f const	&v2 = mesh_data_.vertices[vertex_index_[2]];

	maths::Point3f	p0 = v0;
	maths::Point3f	p1 = v1;
	maths::Point3f	p2 = v2;

	/*	This is fun to read, but it's probably less efficient than an unrolled loop.
	maths::Vector<maths::Vec3f, 3>	vertices{ maths::Vec3f(p0), maths::Vec3f(p1), maths::Vec3f(p2) };
	vertices -= maths::Vec3f(_ray.origin);
	*/

	maths::Vec3f const	ro = maths::Vec3f(_ray.origin);
	p0 -= ro; p1 -= ro; p2 -= ro;

	uint32_t const		kz = maths::MaximumDimension(maths::Abs(_ray.direction));
	uint32_t const		kx = kz + 1 == 3 ? 0 : kz + 1;
	uint32_t const		ky = kx + 1 == 3 ? 0 : kx + 1;
	maths::Vec3f const	rd = maths::Swizzle(_ray.direction, kx, ky, kz);
	p0 = maths::Swizzle(p0, kx, ky, kz);
	p1 = maths::Swizzle(p1, kx, ky, kz);
	p2 = maths::Swizzle(p2, kx, ky, kz);

	maths::Decimal const	sx = -rd.x / rd.z;
	maths::Decimal const	sy = -rd.y / rd.z;
	maths::Decimal const	sz = 1._d / rd.z;
	maths::Vec3f const		shear{ sx, sy, 0._d };
	p0 += shear * p0.z;
	p1 += shear * p1.z;
	p2 += shear * p2.z;

	// Edge function : e(p, p0, p1) = (p1.x - p0.x)(p.y - p0.y) - (p.x - p0.x)(p1.y - p0.y)
	// Returns a positive value if p is on the left side of line p0p1 and a negative value if it's
	// on the right side
	// Thanks to the earlier coordinate system change, p = (0, 0)
	// (p1.x - p0.x)(-p0.y) - (-p0.x)(p1.y - p0.y)
	// p0.y*p0.x - p1.x*p0.y - (p0.y*p0.x - p0.x*p1.y)
	// p0.x*p1.y - p1.x*p0.y
	maths::Decimal	e0 = p1.x*p2.y - p2.x*p1.y;
	maths::Decimal	e1 = p2.x*p0.y - p0.x*p2.y;
	maths::Decimal	e2 = p0.x*p1.y - p1.x*p0.y;
	// NOTE: If any of these is zero, we don't know for sure whether we hit the triangle or not.
	if (e0 == 0._d || e1 == 0._d || e2 == 0._d)
		return false;

	if ((e0 > 0._d || e1 > 0._d || e2 > 0._d) && (e0 < 0._d || e1 < 0._d || e2 < 0._d))
		return false;
	maths::Decimal const	e_sum = e0 + e1 + e2;
	if (e_sum == 0._d)
		return false;

	p0.z *= sz; p1.z *= sz; p2.z *= sz;
	maths::Decimal const	t_scaled = e0*p0.z + e1*p1.z + e2*p2.z;
	if (e_sum < 0._d && (t_scaled >= 0._d || t_scaled < _ray.tMax * e_sum))
		return false;
	if (e_sum > 0._d && (t_scaled <= 0._d || t_scaled > _ray.tMax * e_sum))
		return false;

	maths::Decimal const	e_sum_inverse = 1._d / e_sum;
	maths::Decimal const	b0 = e0 * e_sum_inverse;
	maths::Decimal const	b1 = e1 * e_sum_inverse;
	maths::Decimal const	b2 = e2 * e_sum_inverse;
	maths::Decimal const	t = t_scaled * e_sum_inverse;

	// Should definitely investigate more on this. 
	// Reference is PBR, 3.9.6, Pharr et al
	maths::Decimal const	maxXt = maths::MaximumComponent(maths::Abs(
		maths::Vec3f{ p0.x, p1.x, p2.x }));
	maths::Decimal const	maxYt = maths::MaximumComponent(maths::Abs(
		maths::Vec3f{ p0.y, p1.y, p2.y }));
	maths::Decimal const	maxZt = maths::MaximumComponent(maths::Abs(
		maths::Vec3f{ p0.z, p1.z, p2.z }));
	maths::Decimal const	deltaX = maths::gamma(5u) * (maxXt + maxZt);
	maths::Decimal const	deltaY = maths::gamma(5u) * (maxYt + maxZt);
	maths::Decimal const	deltaZ = maths::gamma(3u) * maxZt;
	maths::Decimal const	deltaE = 2._d * (maths::gamma(2u) * maxXt * maxYt +
											 deltaY * maxXt +
											 deltaX * maxYt);
	maths::Decimal const	maxE = maths::MaximumComponent(maths::Abs(
		maths::Vec3f{ e0, e1, e2 }));
	maths::Decimal const	deltaT = 3._d * (maths::gamma(3u) * maxE * maxZt +
											 deltaE * maxZt +
											 deltaZ * maxE) * maths::Abs(e_sum_inverse);
	if (t <= deltaT)
		return false;

	maths::Vec3f			dpdu, dpdv;
	maths::Point2f const	uv0{ uv(0) }, uv1{ uv(1) }, uv2{ uv(2) };
	maths::Vec2f const		duv02 = uv0 - uv2, duv12 = uv1 - uv2;
	maths::Vec3f const		dp02 = v0 - v2, dp12 = v1 - v2;

	maths::Decimal const	x_abs_sum =
		maths::Abs(b0 * p0.x) + maths::Abs(b1 * p1.x) + maths::Abs(b2 * p2.x);
	maths::Decimal const	y_abs_sum =
		maths::Abs(b0 * p0.y) + maths::Abs(b1 * p1.y) + maths::Abs(b2 * p2.y);
	maths::Decimal const	z_abs_sum =
		maths::Abs(b0 * p0.z) + maths::Abs(b1 * p1.z) + maths::Abs(b2 * p2.z);
	maths::Vec3f const error_bounds =
		maths::gamma(7) * maths::Vec3f{ x_abs_sum, y_abs_sum, z_abs_sum };

	maths::Decimal const	matrix_determinant = duv02.x * duv12.y - duv02.y * duv12.x;
	if (matrix_determinant != 0._d)
	{
		maths::Decimal const	matrix_determinant_inverse = 1._d / matrix_determinant;
		dpdu = (duv12.y * dp02 - duv02.y * dp12) * matrix_determinant_inverse;
		dpdv = (duv02.x * dp12 - duv12.x * dp02) * matrix_determinant_inverse;
	}
	else
		maths::OrthonormalBasis(maths::Normalized(maths::Cross(v2 - v0, v1 - v0)), dpdu, dpdv);

	maths::Point3f const	hit_point = b0 * v0 + maths::Vec3f(b1 * v1) + maths::Vec3f(b2 * v2);
	maths::Point2f const	hit_uv = b0 * uv0 + maths::Vec2f(b1 * uv1) + maths::Vec2f(b2 * uv2);

	maths::Norm3f const flat_normal = maths::Norm3f{ maths::Normalized(maths::Cross(dp02, dp12)) };
	maths::Norm3f geometry_normal = flat_normal;
	SurfaceInteraction::GeometryProperties shading{
		maths::Normalized(dpdu), maths::Normalized(dpdv), maths::Norm3f(0._d), maths::Norm3f(0._d)
	};

	bool const	mesh_has_normals = mesh_data_.has_normals();
	bool const	mesh_has_tangents = mesh_data_.has_tangents();
	if (mesh_has_normals || mesh_has_tangents)
	{
		//
		if (mesh_has_normals)
		{
			maths::Norm3f const		&n0 = mesh_data_.normals[vertex_index_[0]];
			maths::Norm3f const		&n1 = mesh_data_.normals[vertex_index_[1]];
			maths::Norm3f const		&n2 = mesh_data_.normals[vertex_index_[2]];
			shading.SetNormal(maths::Normalized(b0 * n0 + b1 * n1 + b2 * n2));
			//
			if (matrix_determinant != 0._d)
			{
				maths::Norm3f const		dn02 = n0 - n2;
				maths::Norm3f const		dn12 = n1 - n2;
				maths::Decimal const	matrix_determinant_inverse = 1._d / matrix_determinant;
				shading.SetDndu((duv12.y * dn02 - duv02.y * dn12) * matrix_determinant_inverse);
				shading.SetDndv((duv02.x * dn12 - duv12.x * dn02) * matrix_determinant_inverse);
			}
		}
		if (mesh_has_tangents)
		{
			shading.SetDpdu(maths::Normalized(
				b0 * mesh_data_.tangents[vertex_index_[0]] +
				b1 * mesh_data_.tangents[vertex_index_[1]] +
				b2 * mesh_data_.tangents[vertex_index_[2]]));
		}
		// NOTE: I suppose this line can be moved inside the condition, but still have to check
		shading.SetDpdv(maths::Cross(shading.normal_quick(), shading.dpdu_quick()));
		if (maths::SqrLength(shading.dpdu_quick()) > 0._d)
		{
			shading.SetDpdv(shading.dpdv());
			shading.SetDpdu(maths::Cross(shading.dpdv_quick(), shading.normal_quick()));
		}
		else
		{
			maths::Vec3f shading_tangent = shading.dpdu_quick();
			maths::Vec3f shading_bitangent = shading.dpdv_quick();
			maths::OrthonormalBasis((maths::Vec3f)shading.normal_quick(),
									shading_tangent,
									shading_bitangent);
			shading.SetDpdu(shading_tangent);
			shading.SetDpdv(shading_bitangent);
		}
		geometry_normal = maths::FaceForward(flat_normal, shading.normal_quick());
	} // (mesh_has_normals || mesh_has_tangents)
	else if (flip_normals ^ swaps_handedness)
	{
		geometry_normal = -geometry_normal;
		shading.SetNormal(geometry_normal);
	}

	_tHit = t;

	SurfaceInteraction::GeometryProperties const geometry{
		geometry_normal, dpdu, dpdv, maths::Norm3f(0._d), maths::Norm3f(0._d)
	};
	_hit_info = SurfaceInteraction{
		hit_point, error_bounds, t, -_ray.direction, this, hit_uv, geometry, shading
	};
	//if (maths::Dot(_ray.direction, _hit_info.shading.normal) > 0._d)
	//	return false;

	return true;
}


bool
Triangle::DoesIntersect(maths::Ray const &_ray) const
{
	LOG_ERROR(tools::kChannelGeneral, "Call to unimplemented Triangle::DoesIntersect. Returning false.");
	YS_ASSERT(false);
	return false;
}


maths::Vec2f
UniformSampleTriangle(maths::Vec2f const &_ksi)
{
	maths::Decimal const ksi1_sqrt = std::sqrt(maths::ExtendToOne(_ksi.x));
	return maths::Vec2f(1._d - ksi1_sqrt, maths::ExtendToOne(_ksi.y) * ksi1_sqrt);
}

Shape::SurfacePoint
Triangle::SampleSurface(maths::Vec2f const &_ksi) const
{
	maths::Point3f position{};
	maths::Norm3f normal{};
	maths::Vec3f position_error{};
	maths::Vec2f const barycentric = UniformSampleTriangle(_ksi);
	maths::Decimal const barycentric_z = 1._d - barycentric.x - barycentric.y;
	maths::Point3f const &p0 = mesh_data_.vertices[vertex_index_[0]];
	maths::Point3f const &p1 = mesh_data_.vertices[vertex_index_[1]];
	maths::Point3f const &p2 = mesh_data_.vertices[vertex_index_[2]];
	maths::Vec3f const b0{ barycentric.x * p0 };
	maths::Vec3f const b1{ barycentric.y * p1 };
	maths::Vec3f const b2{ barycentric_z * p2 };
	position = maths::Point3f{ b0 + b1 + b2 };
	if (mesh_data_.has_normals())
	{
		normal = maths::Normalized(
			barycentric.x * mesh_data_.normals[vertex_index_[0]] +
			barycentric.y * mesh_data_.normals[vertex_index_[1]] +
			barycentric_z * mesh_data_.normals[vertex_index_[2]]);
	}
	else
	{
		normal = maths::Norm3f{ maths::Normalized(maths::Cross(p1 - p0, p2 - p0)) };
	}
	if (flip_normals)
	{
		normal = -normal;
	}
	position_error = maths::gamma(6u) * (maths::Abs(b0) + maths::Abs(b1) + maths::Abs(b2));
	return SurfacePoint{ position, normal, position_error };
}


maths::Decimal
Triangle::Area() const
{
	maths::Point3f const	&p0 = mesh_data_.vertices[vertex_index_[0]];
	maths::Point3f const	&p1 = mesh_data_.vertices[vertex_index_[1]];
	maths::Point3f const	&p2 = mesh_data_.vertices[vertex_index_[2]];
	return .5_d * maths::Length(maths::Cross(p1 - p0, p2 - p0));
}

maths::Bounds3f
Triangle::ObjectBounds() const
{
	return maths::Union(
		maths::Bounds3f{
			world_transform(mesh_data_.vertices[vertex_index_[0]], maths::Transform::kInverse),
			world_transform(mesh_data_.vertices[vertex_index_[1]], maths::Transform::kInverse)
		},
		world_transform(mesh_data_.vertices[vertex_index_[2]], maths::Transform::kInverse)
	);
}
maths::Bounds3f
Triangle::WorldBounds() const
{
	return maths::Union(
		maths::Bounds3f{
			mesh_data_.vertices[vertex_index_[0]],
			mesh_data_.vertices[vertex_index_[1]]
		},
		mesh_data_.vertices[vertex_index_[2]]
	);
}

maths::Point2f
Triangle::uv(uint32_t _index) const
{
	YS_ASSERT(_index < 3u);
	if (mesh_data_.has_uvs())
		return mesh_data_.uvs[vertex_index_[_index]];
	else
		return maths::Point2f(_index < 1u ? 0._d : 1._d, _index < 2u ? 0._d : 1._d);
}

} // namespace raytracer
