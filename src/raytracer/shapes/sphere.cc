#include "raytracer/shapes/sphere.h"

#include "globals.h"
#include "core/logger.h"
#include "maths/bounds.h"
#include "maths/transform.h"
#include "maths/redecimal.h"
#include "raytracer/surface_interaction.h"
#include "core/profiler.h"


namespace raytracer {

Sphere::Sphere(maths::Transform const &_world_transform, bool _flip_normals,
			   maths::Decimal _radius, maths::Decimal _z_min, maths::Decimal _z_max,
			   maths::Decimal _phi_max) :
	Shape(_world_transform, _flip_normals),
	radius{ _radius },
	z_min{ maths::Clamp(maths::Min(_z_min, _z_max), -radius, radius) },
	z_max{ maths::Clamp(maths::Max(_z_min, _z_max), -radius, radius) },
	theta_min{ std::acos(maths::Clamp(z_min / radius, -1._d, 1._d)) },
	theta_max{ std::acos(maths::Clamp(z_max / radius, -1._d, 1._d)) },
	phi_max{ maths::Radians(maths::Clamp(_phi_max, 0._d, 360._d)) }
{}

bool
Sphere::Intersect(maths::Ray const &_ray,
				  maths::Decimal &_tHit,
				  SurfaceInteraction &_hit_info) const
{
	TIMED_SCOPE(Sphere_Intersect);

	maths::Vec3f origin_error{ maths::zero<maths::Vec3f> }, direction_error{ maths::zero<maths::Vec3f> };
	maths::Ray ray = world_transform(_ray,
									 origin_error, direction_error,
									 maths::Transform::kInverse);

	// xx + yy + zz - rr = 0
	// (ox + tdx)^2 + (oy + tdy)^2 + (oz + tdz)^2 = rr
	// oxox + 2oxdxt + dxdxt^2 + oyoy + 2oydyt + dydyt^2 + ozoz + 2ozdzt + dzdzt^2 = rr
	// (dxdx + dydy + dzdz)t^2 + 2(oxdx + oydy + ozdz)t + (oxox + oyoy + ozoz - rr) = 0
	maths::REDecimal const	ox{ ray.origin.x, origin_error.x },
							oy{ ray.origin.y, origin_error.y },
							oz{ ray.origin.z, origin_error.z };
	maths::REDecimal const	dx{ ray.direction.x, direction_error.x },
							dy{ ray.direction.y, direction_error.y },
							dz{ ray.direction.z, direction_error.z };
	maths::REDecimal const A = dx * dx + dy * dy + dz * dz;
	maths::REDecimal const B = maths::REDecimal(2._d) * (ox * dx + oy * dy + oz * dz);
	maths::REDecimal const re_radius = maths::REDecimal(radius);
	maths::REDecimal const C = ox * ox + oy * oy + oz * oz - re_radius * re_radius;
	maths::REDecimal	t0, t1;
	if (!maths::Quadratic(A, B, C, t0, t1))
		return false;

	if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0._d)
		return false;
	
	// NOTE: The way this is done feels weird to me. We want to test if t0 is valid, then t1, but
	//		 it is implemented in a strange nested way.
	//		 I'd much rather have it implemented in a loop of some sort, testing all of the
	//		 conditions in a sequence.
	maths::REDecimal tHit = t0;
	if (tHit.LowerBound() <= 0._d)
	{
		tHit = t1;
		if (tHit.UpperBound() > ray.tMax)
			return false;
	}

	maths::Point3f		pHit{ ray(tHit.value) };
	pHit *= radius / maths::Distance(pHit, maths::zero<maths::Point3f>);
	if (pHit.x == 0._d && pHit.y == 0._d)
		pHit.x = 1e-5f * radius;
	maths::Decimal		phi = std::atan2(pHit.x, pHit.y);
	if (phi < 0._d)
		phi += 2._d * maths::pi<maths::Decimal>;

	if ((z_min > -radius && pHit.z < z_min) ||
		(z_max < radius && pHit.z > z_max) ||
		phi > phi_max)
	{
		if (tHit == t1) return false;
		if (t1.UpperBound() > ray.tMax) return false;
		tHit = t1;

		pHit = ray(tHit.value);
		pHit *= radius / maths::Distance(pHit, maths::zero<maths::Point3f>);
		if (pHit.x == 0._d && pHit.y == 0._d)
			pHit.x = 1e-5f * radius;
		phi = std::atan2(pHit.x, pHit.y);
		if (phi < 0._d)
			phi += 2._d * maths::pi<maths::Decimal>;

		if ((z_min > -radius && pHit.z < z_min) ||
			(z_max < radius && pHit.z > z_max) ||
			phi > phi_max)
			return false;
	}

	maths::Decimal const theta_delta = theta_max - theta_min;
	maths::Decimal const u = phi / phi_max;
	maths::Decimal const theta = std::acos(maths::Clamp(pHit.z / radius, -1._d, 1._d));
	maths::Decimal const v = (theta - theta_min) / theta_delta;

	maths::Decimal const z_radius = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
	maths::Decimal const inv_z_radius = 1._d / z_radius;
	maths::Decimal const cos_phi = pHit.x * inv_z_radius;
	maths::Decimal const sin_phi = pHit.y * inv_z_radius;
	maths::Vec3f const dpdu(-phi_max * pHit.y, phi_max * pHit.x, 0._d);
	maths::Vec3f const dpdv{
		theta_delta *
		maths::Vec3f{ pHit.z * cos_phi, pHit.z * sin_phi, -radius * std::sin(theta) }
	};

	maths::Vec3f const error_bounds = maths::gamma(5) * maths::Abs(maths::Vec3f(pHit));

	// NOTE: Using differential geometry first fundamental form to compute dndu and dndv
	//		 see Gray(1991) for a reference on differential geometry.
	maths::Vec3f const d2pduu = -phi_max * phi_max * maths::Vec3f{ pHit.x, pHit.y, 0._d };
	maths::Vec3f const d2pduv =
		theta_delta * pHit.z * phi_max * maths::Vec3f{ -sin_phi, cos_phi, 0._d };
	maths::Vec3f const d2pdvv = -theta_delta * theta_delta * (maths::Vec3f)pHit;

	maths::Decimal const E = Dot(dpdu, dpdu);
	maths::Decimal const F = Dot(dpdu, dpdv);
	maths::Decimal const G = Dot(dpdv, dpdv);
	maths::Vec3f const N = maths::Normalized(maths::Cross(dpdu, dpdv));
	maths::Decimal const e = Dot(N, d2pduu);
	maths::Decimal const f = Dot(N, d2pduv);
	maths::Decimal const g = Dot(N, d2pdvv);

	maths::Decimal const inv_EGF2 = 1._d / (E * G - F * F);
	maths::Norm3f const dndu = maths::Norm3f{
		(f * F - e * G) * inv_EGF2 * dpdu + (e * F - f * E) * inv_EGF2 * dpdv
	};
	maths::Norm3f const dndv = maths::Norm3f{
		(g * F - f * G) * inv_EGF2 * dpdu + (f * F - g * E) * inv_EGF2 * dpdv
	};

	_hit_info = world_transform(SurfaceInteraction(
		pHit, error_bounds, ray.time, -ray.direction, this, maths::Point2f(u, v), dpdu, dpdv, dndu, dndv
	), maths::Transform::kForward);

	_tHit = tHit.value;

	return true;
}
bool
Sphere::DoesIntersect(maths::Ray const &_ray) const
{
	LOG_ERROR(tools::kChannelGeneral, "Call to unimplemented Triangle::DoesIntersect. Returning false.");
	YS_ASSERT(false);
	return false;
}

maths::Decimal
Sphere::Area() const
{
	return phi_max * radius * (z_max - z_min);
}

maths::Bounds3f
Sphere::ObjectBounds() const
{
	return maths::Bounds3f{
		{-radius, -radius, z_min},
		{radius, radius, z_max}
	};
}


} // namespace raytracer
