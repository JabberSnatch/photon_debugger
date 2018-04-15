
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

#include <boost/program_options.hpp>

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

struct IntersectInfo;
struct TestResults;
struct TestContext;

TestResults TestRoutine(raytracer::Shape const &_shape, maths::Ray const &_primary_ray);

maths::Point3f PointFromSpherical(maths::Decimal const _theta,
								  maths::Decimal const _phi,
								  maths::Decimal const _distance);


namespace std {
template <typename T, uint32_t n>
std::ostream &operator<<(std::ostream &_ostream, std::array<T, n> const &_array);
std::ostream &operator<<(std::ostream &_ostream, IntersectInfo const &_info);
template <typename T, uint32_t n>
std::ostream &operator<<(std::ostream &_ostream, maths::Vector<T, n> const &_vector)
{
	return _ostream << _vector.e;
}
}

// BOOST EXTENSION
#include <regex>
#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
namespace boost::program_options
{
template <typename T, uint32_t kSize>
void validate(boost::any &_v,
			  std::vector<std::string> const &_values,
			  maths::Vector<T, kSize> *_target_type, int)
{
	bool validate_error = false;
	maths::Vector<T, kSize> result{};

	validators::check_first_occurrence(_v);

	std::string const &str = validators::get_single_string(_values);

	static std::string const number_pattern{ "((\\d)+(\\.)?(\\d)*)" };

	static std::regex const vector_regex{ "(" + number_pattern + "x)*" + number_pattern,
										  std::regex::icase };
	static std::regex const integer_regex{ number_pattern, std::regex::icase };

	std::smatch match{};
	validate_error = !std::regex_match(str, match, vector_regex);
	if (!validate_error)
	{
		std::sregex_iterator const regex_begin = std::sregex_iterator{ str.cbegin(), str.cend(), integer_regex };
		std::sregex_iterator const regex_end{};

		std::ptrdiff_t const match_count = std::distance(regex_begin, regex_end);
		validate_error = match_count != kSize;
		if (!validate_error)
		{
			for (std::sregex_iterator rit = regex_begin; rit != regex_end; ++rit)
			{
				int const index = static_cast<int>(std::distance(regex_begin, rit));
				try {
					result[index] = boost::lexical_cast<T>(rit->str());
				}
				catch(std::exception&)
				{
					std::cout << "wrong number type : " << rit->str() << std::endl;
					validate_error = true;
				}
			}
		}
		else
		{
			std::cout << "wrong size" << std::endl;
		}
	}
	else
	{
		std::cout << "wrong format" << std::endl;
	}

	if (validate_error)
	{
		throw validation_error(validation_error::invalid_option_value);
	}

	_v = boost::any(result);
}
}
// <<<<<<< BOOST EXTENSION

struct IntersectInfo
{
	maths::Decimal t = -1._d;
	raytracer::SurfaceInteraction hit{};
	bool intersect = false;
};


struct TestResults
{
	IntersectInfo primary{};
	IntersectInfo inward{};
	IntersectInfo outward{};
};


constexpr maths::Decimal	ktMax = maths::infinity<maths::Decimal>;
constexpr maths::Decimal	kTime = 0._d;

constexpr maths::Decimal	kQuadExtentX = 1._d;
constexpr maths::Decimal	kQuadExtentY = 1._d;
constexpr maths::Decimal	kQuadAbsBoundX = kQuadExtentX / 2._d;
constexpr maths::Decimal	kQuadAbsBoundY = kQuadExtentY / 2._d;

constexpr maths::Decimal	kGeomCoeff = 1._d;

constexpr maths::Decimal	kThetaMax = 90._d;
constexpr maths::Decimal	kPhiMax = 45._d;

constexpr int				kThetaSubdiv = 4;
constexpr int				kPhiSubdiv = 2;
constexpr maths::Decimal	kDistance = 2._d;

constexpr int				kXSubdiv = 16;
constexpr int				kYSubdiv = 16;

constexpr int				kExecCount = kThetaSubdiv * kPhiSubdiv * kXSubdiv * kYSubdiv;


struct TestContext
{
	maths::Decimal			tMax{ ktMax };
	maths::Decimal			time{ kTime };

	maths::Vec2f			target_extents{ kQuadExtentX, kQuadExtentY };
	maths::Vector2<int>		target_subdiv{ kXSubdiv, kYSubdiv };
	maths::Decimal			quad_coeff{ kGeomCoeff };

	maths::Vec2f			origin_extents{ kThetaMax, kPhiMax };
	maths::Vector2<int>		origin_subdiv{ kThetaSubdiv, kPhiSubdiv };
	maths::Decimal			origin_distance{ kDistance };

	int						exec_count() const
	{ return maths::FoldProduct(target_subdiv * origin_subdiv); }

	maths::Vec2f			target_step_size() const
	{ return target_extents / static_cast<maths::Vec2f>(target_subdiv); }
	maths::Vec2f			origin_step_size() const
	{ return origin_extents / static_cast<maths::Vec2f>(origin_subdiv); }
	maths::Vec2f			target_abs_bounds() const
	{ return target_extents / maths::Vec2f{ 2._d }; }

	TriangleMesh const		&quad();

private:
	std::unique_ptr<raytracer::TriangleMeshRawData>	quad_raw_data_cache_;
	std::unique_ptr<raytracer::TriangleMeshData>	quad_data_cache_;
	std::unique_ptr<TriangleMesh>					quad_cache_;
};


int main(int argc, char **argv)
{
	using ResultsContainer_t = std::vector<TestResults>;

	TestContext test_context{};

	{
		namespace bpo = boost::program_options;

		bpo::options_description option_desc{ "available options" };
		option_desc.add_options()
			("help", "help message")
			("target-extents",
			 bpo::value(&test_context.target_extents)->default_value(test_context.target_extents),
			 "float")
			("target-subdiv",
			 bpo::value(&test_context.target_subdiv)->default_value(test_context.target_subdiv),
			 "int")
			("quad-coeff",
			 bpo::value(&test_context.quad_coeff)->default_value(test_context.quad_coeff),
			 "float")
			("origin-extents",
			 bpo::value(&test_context.origin_extents)->default_value(test_context.origin_extents),
			 "float - extents of the ray's origin expressed in degrees in spherical coordinates (theta, phi)")
			("origin-subdiv",
			 bpo::value(&test_context.origin_subdiv)->default_value(test_context.origin_subdiv),
			 "int")
			("origin-distance",
			 bpo::value(&test_context.origin_distance)->default_value(test_context.origin_distance),
			 "float");

		bpo::variables_map vm{};
		try
		{
			bpo::store(bpo::parse_command_line(argc, argv, option_desc), vm);
			bpo::notify(vm);
		}
		catch (bpo::error &e)
		{
			std::cout << e.what() << std::endl;
		}

		if (vm.count("help"))
		{
			std::cout << option_desc << std::endl;
			return 1;
		}

		std::cout << test_context.target_extents << std::endl;
		std::cout << test_context.target_subdiv << std::endl;
		std::cout << test_context.quad_coeff << std::endl;
		std::cout << test_context.origin_extents << std::endl;
		std::cout << test_context.origin_subdiv << std::endl;
		std::cout << test_context.origin_distance << std::endl;
		std::cout << test_context.exec_count() << std::endl << std::endl;;
	}

	TriangleMesh const &quad = test_context.quad();

	ResultsContainer_t results_container{};
	results_container.reserve(test_context.exec_count());

	{
		maths::Vec2f const origin_step_size = test_context.origin_step_size();
		maths::Vec2f const target_step_size = test_context.target_step_size();
		maths::Vec2f const target_abs_bounds = test_context.target_abs_bounds();

		for (int theta_step = 0; theta_step < test_context.origin_subdiv.x; ++theta_step)
		{
			maths::Decimal const theta = maths::Radians(theta_step * origin_step_size.x);

			for (int phi_step = 0; phi_step < test_context.origin_subdiv.y; ++phi_step)
			{
				maths::Decimal const phi = maths::Radians(phi_step * origin_step_size.y);
				maths::Point3f const origin = PointFromSpherical(theta, phi, test_context.origin_distance);

				for (int x_step = 0; x_step < test_context.target_subdiv.x; ++x_step)
				{
					maths::Decimal const x_target =
						-target_abs_bounds.x +
						target_step_size.x * x_step +
						target_step_size.x * .5_d;

					for (int y_step = 0; y_step < test_context.target_subdiv.y; ++y_step)
					{
						maths::Decimal const y_target =
							-target_abs_bounds.y +
							target_step_size.y * y_step +
							target_step_size.y * .5_d;
						maths::Point3f const target{ x_target, y_target, 0._d };
						maths::Vec3f const direction{ maths::Normalized(target - origin) };
#ifdef TESTROUTINE_LOGGING
						std::cout << origin.e << std::endl;
						std::cout << target.e << std::endl;
						std::cout << direction.e << std::endl;
#endif

						maths::Ray const primary_ray{ origin, direction, ktMax, kTime };
						results_container.emplace_back(TestRoutine(quad, primary_ray));

#ifdef TESTROUTINE_LOGGING
						std::cout << std::endl;
#endif
					}
				}
			}
		}
	}

	{
		int64_t const primary_fail_count = std::count_if(
			results_container.cbegin(), results_container.cend(),
			[](TestResults const &_results)
			{
				return !_results.primary.intersect;
			});
		int64_t const inward_fail_count = std::count_if(
			results_container.cbegin(), results_container.cend(),
			[](TestResults const &_results)
			{
				return !_results.inward.intersect;
			});
		int64_t const outward_fail_count = std::count_if(
			results_container.cbegin(), results_container.cend(),
			[](TestResults const &_results)
			{
				return _results.outward.intersect;
			});

		maths::Decimal const primary_fail_ratio =
			static_cast<maths::Decimal>(primary_fail_count) /
			static_cast<maths::Decimal>(test_context.exec_count());
		maths::Decimal const inward_fail_ratio =
			static_cast<maths::Decimal>(inward_fail_count) /
			static_cast<maths::Decimal>(test_context.exec_count());
		maths::Decimal const outward_fail_ratio =
			static_cast<maths::Decimal>(outward_fail_count) /
			static_cast<maths::Decimal>(test_context.exec_count());

		std::cout << "SUMMARY" << std::endl;
		std::cout << results_container.size() << std::endl;
		std::cout << "TESTED " << test_context.exec_count() << " CASES" << std::endl;
		std::cout << "PRIMARY FAILS .. " << primary_fail_count << " .. "
				  << primary_fail_ratio << std::endl;
		std::cout << "INWARD FAILS .. " << inward_fail_count << " .. "
				  << inward_fail_ratio << std::endl;
		std::cout << "OUTWARD FAILS .. " << outward_fail_count << " .. "
				  << outward_fail_ratio << std::endl;

#if 0
		std::cout << "INWARD FAILS LIST .. " << std::endl;
		std::for_each(results_container.cbegin(), results_container.cend(),
					  [](TestResults const &_results)
					  {
						  if (!_results.inward.intersect)
						  {
							  std::cout << _results.primary << std::endl << std::endl;
						  }
					  });
		std::cout << "INWARD SUCCESS SAMPLE .. " << std::endl;
		std::find_if(std::next(results_container.cbegin(), results_container.size() / 2), results_container.cend(),
					 [](TestResults const &_results)
					 {
						 if (_results.inward.intersect)
						 {
							 std::cout << _results.primary << std::endl << std::endl;
							 return true;
						 }
						 return false;
					 });
#endif
	}

	return 0;
}


TestResults TestRoutine(raytracer::Shape const &_shape, maths::Ray const &_primary_ray)
{
	TestResults result{};

	result.primary.intersect = _shape.Intersect(_primary_ray, result.primary.t, result.primary.hit);

#ifdef TESTROUTINE_LOGGING
	std::cout << "PRIMARY.." << std::endl;
	std::cout << result.primary << std::endl;
#endif

	if (result.primary.intersect)
	{
#ifdef TESTROUTINE_LOGGING
		std::cout << "OK" << std::endl;
#endif

		maths::Vec3f const primary_offset_w{ result.primary.hit.geometry.normal() };

		maths::Point3f const secondary_origin =
			result.primary.hit.OffsetOriginFromErrorBounds(primary_offset_w);

		maths::Ray const secondary_inward_ray{ secondary_origin,
											   _primary_ray.direction,
											   ktMax, kTime };
		maths::Ray const secondary_outward_ray{ secondary_origin,
												-_primary_ray.direction,
												ktMax, kTime };

		result.inward.intersect = _shape.Intersect(secondary_inward_ray,
												   result.inward.t,
												   result.inward.hit);
		result.outward.intersect = _shape.Intersect(secondary_outward_ray,
													result.outward.t,
													result.outward.hit);

#ifdef TESTROUTINE_LOGGING
		std::cout << "SECONDARY INWARD.." << std::endl;
		std::cout << result.inward << std::endl;
		if (result.inward.intersect)
			std::cout << "OK" << std::endl;
		else
			std::cout << "FAIL: No inward intersection" << std::endl;

		std::cout << "SECONDARY OUTWARD.." << std::endl;
		std::cout << result.outward << std::endl;
		if (!result.outward.intersect)
			std::cout << "OK" << std::endl;
		else
			std::cout << "FAIL: Outward intersection" << std::endl;
#endif
	}
	else
	{
#ifdef TESTROUTINE_LOGGING
		std::cout << "FAIL: No intersection" << std::endl;
#endif
	}

	return result;
};


TriangleMesh const &
TestContext::quad()
{
	using raytracer::TriangleMeshRawData;
	using raytracer::TriangleMeshData;

	maths::Vec2f bounds = target_abs_bounds() * quad_coeff;

	static TriangleMeshRawData::IndicesContainer_t const indices{
		0, 1, 2,
		2, 3, 0
	};
	TriangleMeshRawData::VerticesContainer_t const vertices{
		maths::Point3f{ bounds.x, bounds.y, 0._d },
		maths::Point3f{ -bounds.x, bounds.y, 0._d },
		maths::Point3f{ -bounds.x, -bounds.y, 0._d },
		maths::Point3f{ bounds.x, -bounds.y, 0._d }
	};
	quad_raw_data_cache_ = std::make_unique<TriangleMeshRawData>(
		2,
		indices,
		vertices
		);

	static maths::Transform const world_transform{
		maths::Mat4x4f::Identity()
	};

	static const bool flip_normals = false;

	quad_data_cache_ = std::make_unique<TriangleMeshData>(
		world_transform,
		flip_normals,
		*quad_raw_data_cache_,
		InstancingPolicy{}
		);

	quad_cache_ = std::make_unique<TriangleMesh>(
		world_transform,
		flip_normals,
		*quad_data_cache_);

	return *quad_cache_;
}


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


namespace std {
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

std::ostream &operator<<(std::ostream &_ostream, IntersectInfo const &_info)
{
	_ostream << _info.intersect << "; " << _info.t << std::endl
			 << _info.hit.position.e << "; " << _info.hit.position_error << std::endl
			 << _info.hit.geometry.normal().e << std::endl
			 << _info.hit.wo;
	return _ostream;
}
}
