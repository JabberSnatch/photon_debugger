
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>

#include "maths/maths.h"
#include "maths/matrix.h"
#include "maths/point.h"
#include "maths/ray.h"
#include "maths/transform.h"
#include "maths/vector.h"

#include "raytracer/primitive.h"
#include "raytracer/triangle_mesh_data.h"
#include "raytracer/shapes/triangle_mesh.h"
#include "raytracer/surface_interaction.h"



using InstancingPolicy = raytracer::InstancingPolicyClass::Transformed;
using TriangleMesh = raytracer::TriangleMesh<InstancingPolicy>;


constexpr maths::Decimal kQuadExtentX = 2._d;
constexpr maths::Decimal kQuadExtentY = 2._d;
constexpr maths::Decimal kQuadAbsBoundX = kQuadExtentX / 2._d;
constexpr maths::Decimal kQuadAbsBoundY = kQuadExtentY / 2._d;

static TriangleMesh const &
test_quad()
{
	using raytracer::TriangleMeshRawData;
	using raytracer::TriangleMeshData;

	static TriangleMeshRawData::IndicesContainer_t const indices{
		0, 1, 2,
		2, 3, 0
	};
	static TriangleMeshRawData::VerticesContainer_t const vertices{
		{ kQuadAbsBoundX, kQuadAbsBoundY, 0._d },
		{ -kQuadAbsBoundX, kQuadAbsBoundY, 0._d },
		{ -kQuadAbsBoundX, -kQuadAbsBoundY, 0._d },
		{ kQuadAbsBoundX, -kQuadAbsBoundY, 0._d }
	};
	static TriangleMeshRawData const raw_data{
		2,
		indices,
		vertices
	};

	static maths::Transform const world_transform{
		maths::Mat4x4f::Identity()
	};

	static const bool flip_normals = false;
	static TriangleMeshData const data{
		world_transform,
		flip_normals,
		raw_data,
		InstancingPolicy{}
	};

	static TriangleMesh const result{
		world_transform,
		flip_normals,
		data
	};

	return result;
}


template <typename T, uint32_t n>
std::ostream &operator<<(std::ostream &_ostream, std::array<T, n> const &_array)
{
	std::array<T, n>::const_iterator const last_elem = std::prev(_array.end());
	std::for_each(_array.cbegin(), last_elem, [&_ostream](T const &_v) {
		_ostream << _v << " ";
	});
	_ostream << *(last_elem);
	return _ostream;
}


std::ostream &OutputHitInfo(std::ostream &_ostream,
							raytracer::SurfaceInteraction const &_hit_info,
							bool const _intersect,
							maths::Decimal const _tHit)
{
	_ostream << _intersect << "; " << _tHit << "; "
			 << _hit_info.position.e << "; "
			 << _hit_info.geometry.normal().e;
	return _ostream;
}



constexpr maths::Decimal ktMax = maths::infinity<maths::Decimal>;
constexpr maths::Decimal kTime = 0._d;

void TestRoutine(raytracer::Shape const &_shape, maths::Ray const &_primary_ray)
{
	maths::Decimal primary_t = -1._d;
	raytracer::SurfaceInteraction primary_hit{};
	bool const primary_intersect = _shape.Intersect(_primary_ray, primary_t, primary_hit);

	std::cout << "PRIMARY.." << std::endl;
	OutputHitInfo(std::cout, primary_hit, primary_intersect, primary_t) << std::endl;

	if (primary_intersect)
	{
		std::cout << "OK" << std::endl;

		maths::Vec3f const primary_offset_w{ primary_hit.geometry.normal() };
		maths::Point3f const secondary_origin =
			primary_hit.OffsetOriginFromErrorBounds(primary_offset_w);

		maths::Ray const secondary_inward_ray{ secondary_origin,
											   _primary_ray.direction,
											   ktMax, kTime };
		maths::Ray const secondary_outward_ray{ secondary_origin,
												-_primary_ray.direction,
												ktMax, kTime };

		maths::Decimal secondary_inward_t{ -1._d }, secondary_outward_t{ -1._d };
		raytracer::SurfaceInteraction secondary_inward_hit{}, secondary_outward_hit{};

		bool const inward_intersect = _shape.Intersect(secondary_inward_ray,
													   secondary_inward_t,
													   secondary_inward_hit);
		bool const outward_intersect = _shape.Intersect(secondary_outward_ray,
														secondary_outward_t,
														secondary_outward_hit);

		std::cout << "SECONDARY INWARD.." << std::endl;
		OutputHitInfo(std::cout, secondary_inward_hit, inward_intersect, secondary_inward_t)
			<< std::endl;
		if (inward_intersect)
			std::cout << "OK" << std::endl;
		else
			std::cout << "FAIL: No inward intersection" << std::endl;

		std::cout << "SECONDARY OUTWARD.." << std::endl;
		OutputHitInfo(std::cout, secondary_outward_hit, outward_intersect, secondary_outward_t)
			<< std::endl;
		if (!outward_intersect)
			std::cout << "OK" << std::endl;
		else
			std::cout << "FAIL: Outward intersection" << std::endl;
	}
	else
	{
		std::cout << "FAIL: No intersection" << std::endl;
	}
};


maths::Point3f PointFromSpherical(maths::Decimal const _theta,
								  maths::Decimal const _phi,
								  maths::Decimal const _distance)
{
	return maths::Point3f{
		_distance * std::sin(_theta) * std::cos(_phi),
		_distance * std::sin(_theta) * std::sin(_phi),
		_distance * std::cos(_theta)
	};
}

constexpr maths::Decimal kThetaMax = maths::pi<maths::Decimal> / 2._d;
constexpr maths::Decimal kPhiMax = maths::pi<maths::Decimal> / 4._d;

constexpr int kThetaSubdiv = 2;
constexpr int kPhiSubdiv = 1;
constexpr maths::Decimal kDistance = 2._d;

constexpr int kXSubdiv = 4;
constexpr int kYSubdiv = 4;

int main()
{
	TriangleMesh const &quad = test_quad();

	maths::Decimal const theta_step_width = kThetaMax / kThetaSubdiv;
	maths::Decimal const phi_step_width = kPhiMax / kPhiSubdiv;
	maths::Decimal const x_step_width = kQuadExtentX / kXSubdiv;
	maths::Decimal const y_step_width = kQuadExtentY / kYSubdiv;
	for (int theta_step = 0; theta_step < kThetaSubdiv; ++theta_step)
	{
		maths::Decimal const theta = theta_step * theta_step_width;
		for (int phi_step = 0; phi_step <= kPhiSubdiv; ++phi_step)
		{
			maths::Decimal const phi = phi_step * phi_step_width;
			maths::Point3f const origin = PointFromSpherical(theta, phi, kDistance);
			for (int x_step = 0; x_step <= kXSubdiv; ++x_step)
			{
				maths::Decimal const x_target = -kQuadAbsBoundX + x_step_width * x_step;
				for (int y_step = 0; y_step <= kYSubdiv; ++y_step)
				{
					maths::Decimal const y_target = -kQuadAbsBoundY + y_step_width * y_step;
					maths::Point3f const target{ x_target, y_target, 0._d };
					maths::Vec3f const direction{ maths::Normalized(target - origin) };
					std::cout << origin.e << std::endl;
					std::cout << target.e << std::endl;
					std::cout << direction.e << std::endl;

					maths::Ray const primary_ray{ origin, direction, ktMax, kTime };
					TestRoutine(quad, primary_ray);
					
					std::cout << std::endl;
				}
			}
		}
	}

	std::system("pause");
	return 0;
}
