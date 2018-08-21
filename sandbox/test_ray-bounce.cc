
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

#include <boost/program_options.hpp>

#include "globals.h"

#include "core/logger.h"

#include "maths/maths.h"
#include "maths/matrix.h"
#include "maths/point.h"
#include "maths/ray.h"
#include "maths/redecimal.h"
#include "maths/transform.h"
#include "maths/vector.h"

#include "raytracer/primitive.h"
#include "raytracer/triangle_mesh_data.h"
#include "raytracer/shapes/triangle_mesh.h"
#include "raytracer/surface_interaction.h"


// NOTE: pathological case :
// --method experimental --target-offset 1048576x1048576x1
// 2^20 gives a decent amount of error, going beyond is all errors


using InstancingPolicy = raytracer::InstancingPolicyClass::Transformed;
using TriangleMesh = raytracer::TriangleMesh<InstancingPolicy>;

struct IntersectInfo;
struct TestResults;
struct TestContext;

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


namespace watertight_intersect {

template <typename InternalFloat>
struct Triangle {
public:
	using InternalVec3 = maths::Vector3<InternalFloat>;
public:
	bool Intersect(maths::Ray const &_ray,
				   maths::Decimal &_tHit,
				   raytracer::SurfaceInteraction &_hit_info) const;

	std::array<InternalVec3, 3u> vertices;
};

template <typename InternalFloat>
struct TriangleList {
public:
	using Triangle_t = Triangle<InternalFloat>;
public:
	TriangleList(std::initializer_list<Triangle_t> _init_list) : triangles{ _init_list } {}
	bool Intersect(maths::Ray const &_ray,
				   maths::Decimal &_tHit,
				   raytracer::SurfaceInteraction &_hit_info) const;
	std::vector<Triangle_t> triangles;
};

} // namespace watertight_intersect


struct IntersectInfo
{
	maths::Decimal t = -maths::infinity<maths::Decimal>;
	raytracer::SurfaceInteraction hit{};
	bool intersect = false;
};


struct TestResults
{
	maths::Decimal NdotL = -maths::infinity<maths::Decimal>;
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
constexpr maths::Decimal	kDistance = 1._d;

constexpr int				kXSubdiv = 16;
constexpr int				kYSubdiv = 16;

constexpr int				kExecCount = kThetaSubdiv * kPhiSubdiv * kXSubdiv * kYSubdiv;


struct TestContext
{
	bool					short_output{ false };

	maths::Decimal			tMax{ ktMax };
	maths::Decimal			time{ kTime };

	maths::Vec3f			target_offset{ 1._d, 1._d, 1._d };
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

	std::string				method{ "baseline" };

	TriangleMesh const		&quad();
	template <typename InternalFloat>
	watertight_intersect::TriangleList<InternalFloat> wi_quad();

public:
	using ResultsContainer_t = std::vector<TestResults>;
	ResultsContainer_t RunTest();
	void PrintResults(ResultsContainer_t const &_results_container) const;
private:
	template <typename ShapeType>
	ResultsContainer_t RunTest_impl(ShapeType const &_shape) const;
	template <typename ShapeType>
	static TestResults TestRoutine(ShapeType const &_shape, maths::Ray const &_primary_ray);

private:
	std::unique_ptr<raytracer::TriangleMeshRawData>	quad_raw_data_cache_;
	std::unique_ptr<raytracer::TriangleMeshData>	quad_data_cache_;
	std::unique_ptr<TriangleMesh>					quad_cache_;
};


int main(int argc, char **argv)
{


	globals::logger.BindPath(tools::kChannelGeneral, "raybounce_general.txt");

	TestContext test_context{};

	{
		namespace bpo = boost::program_options;

		bpo::options_description option_desc{ "available options" };
		option_desc.add_options()
			("help", "help message")
			("target-offset",
			 bpo::value(&test_context.target_offset)->default_value(test_context.target_offset),
			 "float")
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
			 "float")
			("method",
			 bpo::value(&test_context.method)->default_value(test_context.method),
			 "string")
			("short-output",
			 bpo::value(&test_context.short_output)->default_value(test_context.short_output)->implicit_value(true),
			 "bool")
			;

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

		if (!test_context.short_output)
		{
			std::cout << test_context.target_offset << std::endl;
			std::cout << test_context.target_extents << std::endl;
			std::cout << test_context.target_subdiv << std::endl;
			std::cout << test_context.quad_coeff << std::endl;
			std::cout << test_context.origin_extents << std::endl;
			std::cout << test_context.origin_subdiv << std::endl;
			std::cout << test_context.origin_distance << std::endl;
			std::cout << test_context.exec_count() << std::endl;
			std::cout << test_context.method << std::endl;
			std::cout << std::endl;
		}
	}

	using StdChrono = std::chrono::high_resolution_clock;
	StdChrono::time_point const start = StdChrono::now();
	TestContext::ResultsContainer_t const results_container = test_context.RunTest();
	StdChrono::duration const delta = StdChrono::now() - start;
	if (!test_context.short_output)
	{
		std::cout << "RUNNING TIME .. "
				  << std::chrono::duration_cast<std::chrono::milliseconds>(delta).count() << " ms" << std::endl;
	}
	test_context.PrintResults(results_container);

	return 0;
}


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
		maths::Point3f{ bounds.x, bounds.y, 0._d } + target_offset,
		maths::Point3f{ -bounds.x, bounds.y, 0._d } + target_offset,
		maths::Point3f{ -bounds.x, -bounds.y, 0._d } + target_offset,
		maths::Point3f{ bounds.x, -bounds.y, 0._d } + target_offset
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

template <typename InternalFloat>
watertight_intersect::TriangleList<InternalFloat>
TestContext::wi_quad()
{
	using TriangleType = watertight_intersect::Triangle<InternalFloat>;
	if (!quad_raw_data_cache_)
		quad();

	std::array<TriangleType::InternalVec3, 4> vertices{};
	for (int i = 0; i < 4; ++i)
	{
		maths::Point3f const &source = quad_raw_data_cache_->vertices[i];
		vertices[i] = TriangleType::InternalVec3{ source.x, source.y, source.z };
	}

	watertight_intersect::TriangleList<InternalFloat> const result{
		TriangleType{ vertices[0], vertices[1], vertices[2] },
		TriangleType{ vertices[2], vertices[3], vertices[0] }
	};

	return result;
}

TestContext::ResultsContainer_t
TestContext::RunTest()
{
	if (method == "baseline")
		return RunTest_impl(quad());
	else if (method == "exp_redecimal")
		return RunTest_impl(wi_quad<maths::REDecimal>());
	else if (method == "exp_single")
		return RunTest_impl(wi_quad<float>());
	else if (method == "exp_double")
		return RunTest_impl(wi_quad<double>());
	else
		return ResultsContainer_t{};
}

template <typename ShapeType>
TestContext::ResultsContainer_t
TestContext::RunTest_impl(ShapeType const &_shape) const
{
	ResultsContainer_t results_container{};
	results_container.reserve(exec_count());

	maths::Vec2f const origin_step_size_cache = origin_step_size();
	maths::Vec2f const target_step_size_cache = target_step_size();
	maths::Vec2f const target_abs_bounds_cache = target_abs_bounds();

	for (int theta_step = 0; theta_step < origin_subdiv.x; ++theta_step)
	{
		maths::Decimal const theta = maths::Radians(theta_step * origin_step_size_cache.x);

		for (int phi_step = 0; phi_step < origin_subdiv.y; ++phi_step)
		{
			maths::Decimal const phi = maths::Radians(phi_step * origin_step_size_cache.y);
			maths::Point3f const origin = PointFromSpherical(theta, phi, origin_distance) + target_offset;

			{
				globals::logger.EnableChannel(tools::kChannelGeneral, true);
				std::stringstream message{ "" };
				message << std::endl << "NEW ORIGIN";
				LOG_INFO(tools::kChannelGeneral, message.str());
				globals::logger.EnableChannel(tools::kChannelGeneral, false);
			}
			for (int x_step = 0; x_step < target_subdiv.x; ++x_step)
			{
				maths::Decimal const x_target =
					-target_abs_bounds_cache.x +
					target_step_size_cache.x * x_step +
					target_step_size_cache.x * .5_d;

				for (int y_step = 0; y_step < target_subdiv.y; ++y_step)
				{
					maths::Decimal const y_target =
						-target_abs_bounds_cache.y +
						target_step_size_cache.y * y_step +
						target_step_size_cache.y * .5_d;
					maths::Point3f const target{
						maths::Point3f{ x_target, y_target, 0._d } + target_offset
					};
					maths::Vec3f const direction{ maths::Normalized(target - origin) };
#ifdef TESTROUTINE_LOGGING
					std::cout << origin.e << std::endl;
					std::cout << target.e << std::endl;
					std::cout << direction.e << std::endl;
#endif


					{
						globals::logger.EnableChannel(tools::kChannelGeneral, true);
						std::stringstream message{ "" };
						message << std::endl << "NEW RAY";
						LOG_INFO(tools::kChannelGeneral, message.str());
						globals::logger.EnableChannel(tools::kChannelGeneral, false);
					}

					maths::Ray const primary_ray{ origin, direction, ktMax, kTime };
					results_container.emplace_back(TestRoutine(_shape, primary_ray));

#ifdef TESTROUTINE_LOGGING
					std::cout << std::endl;
#endif
				}
			}
		}
	}

	return results_container;
}

template <typename ShapeType>
TestResults
TestContext::TestRoutine(ShapeType const &_shape, maths::Ray const &_primary_ray)
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
		result.NdotL = maths::Dot(result.primary.hit.geometry.normal(), -_primary_ray.direction);

		maths::Vec3f const primary_offset_w{ result.primary.hit.geometry.normal() };

		maths::Point3f const secondary_origin =
			result.primary.hit.OffsetOriginFromErrorBounds(primary_offset_w);

		maths::Ray const secondary_inward_ray{ secondary_origin,
											   _primary_ray.direction,
											   ktMax, kTime };
		maths::Ray const secondary_outward_ray{ secondary_origin,
												-_primary_ray.direction,
												ktMax, kTime };

		globals::logger.EnableChannel(tools::kChannelGeneral, true);
		result.inward.intersect = _shape.Intersect(secondary_inward_ray,
												   result.inward.t,
												   result.inward.hit);
		globals::logger.EnableChannel(tools::kChannelGeneral, false);
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

void
TestContext::PrintResults(ResultsContainer_t const &_results_container) const
{
	int64_t const primary_fail_count = std::count_if(
		_results_container.cbegin(), _results_container.cend(),
		[](TestResults const &_results)
		{
			return !_results.primary.intersect;
		});
	int64_t const inward_fail_count = std::count_if(
		_results_container.cbegin(), _results_container.cend(),
		[](TestResults const &_results)
		{
			return !_results.inward.intersect;
		});
	int64_t const outward_fail_count = std::count_if(
		_results_container.cbegin(), _results_container.cend(),
		[](TestResults const &_results)
		{
			return _results.outward.intersect;
		});

	maths::Vec3f const average_position = std::accumulate(
		_results_container.cbegin(), _results_container.cend(), maths::Vec3f( 0._d ),
		[_count = _results_container.size()](maths::Vec3f _acc, TestResults const &_results)
		{
			return _acc +
				static_cast<maths::Vec3f>(_results.primary.hit.position) /
				static_cast<maths::Decimal>(_count);
		});

	maths::Vector3<maths::DecimalBits> ulp_worst_error_bounds{ 0u, 0u, 0u };
	std::for_each(
		_results_container.cbegin(), _results_container.cend(),
		[&ulp_worst_error_bounds, &average_position](TestResults const &_results)
		{
			using Vec3DecimalToBits = maths::Vector3<maths::DecimalBitsMapper>;
			maths::Vec3f const position = static_cast<maths::Vec3f>(
				maths::Abs(_results.primary.hit.position));
			maths::Vec3f const position_plus_bounds =
				position + maths::Abs(_results.primary.hit.position_error);
			Vec3DecimalToBits const position_bits{ position.x, position.y, position.z };
			Vec3DecimalToBits const offset_position_bits{ position_plus_bounds.x,
														  position_plus_bounds.y,
														  position_plus_bounds.z };
			ulp_worst_error_bounds.x = std::max(ulp_worst_error_bounds.x,
												offset_position_bits.x.bits - position_bits.x.bits);
			ulp_worst_error_bounds.y = std::max(ulp_worst_error_bounds.y,
												offset_position_bits.y.bits - position_bits.y.bits);
			ulp_worst_error_bounds.z = std::max(ulp_worst_error_bounds.z,
												offset_position_bits.z.bits - position_bits.z.bits);
		});


	maths::Decimal const primary_fail_ratio =
		static_cast<maths::Decimal>(primary_fail_count) /
		static_cast<maths::Decimal>(exec_count());
	maths::Decimal const inward_fail_ratio =
		static_cast<maths::Decimal>(inward_fail_count) /
		static_cast<maths::Decimal>(exec_count());
	maths::Decimal const outward_fail_ratio =
		static_cast<maths::Decimal>(outward_fail_count) /
		static_cast<maths::Decimal>(exec_count());

	if (!short_output)
	{
		std::cout << "SUMMARY" << std::endl;
		std::cout << _results_container.size() << std::endl;
		std::cout << "TESTED " << exec_count() << " CASES" << std::endl;
		std::cout << "PRIMARY FAILS .. " << primary_fail_count << " .. "
				  << primary_fail_ratio << std::endl;
		std::cout << "INWARD FAILS .. " << inward_fail_count << " .. "
				  << inward_fail_ratio << std::endl;
		std::cout << "OUTWARD FAILS .. " << outward_fail_count << " .. "
				  << outward_fail_ratio << std::endl;

		std::cout << "LARGEST ERROR BOUNDS .. " << ulp_worst_error_bounds << std::endl;
		std::cout << "AVERAGE PRIMARY POSITION .. " << average_position << std::endl;
	}
	else
	{
		std::cout << exec_count() << " "
				  << primary_fail_count << " "
				  << inward_fail_count << " "
				  << outward_fail_count << " "
				  << ulp_worst_error_bounds
			//<< " "
			//<< std::endl;
			;
	}
#if 0
	std::cout << "INWARD FAILS LIST .. " << std::endl;
	std::for_each(_results_container.cbegin(), _results_container.cend(),
				  [](TestResults const &_results)
				  {
					  if (!_results.inward.intersect)
					  {
						  std::cout << "NdotL .. " << _results.NdotL << std::endl;
						  std::cout << _results.primary << std::endl << std::endl;
					  }
				  });
	std::cout << "INWARD SUCCESS SAMPLE .. " << std::endl;
	std::find_if(std::next(_results_container.cbegin(), _results_container.size() / 2), _results_container.cend(),
				 [](TestResults const &_results)
				 {
					 if (_results.inward.intersect)
					 {
						 std::cout << "NdotL .. " << _results.NdotL << std::endl;
						 std::cout << _results.primary << std::endl << std::endl;
						 return true;
					 }
					 return false;
				 });
#endif
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


namespace watertight_intersect {

template <typename InternalFloat>
bool Triangle<InternalFloat>::Intersect(maths::Ray const &_ray,
										maths::Decimal &_tHit,
										raytracer::SurfaceInteraction &_hit_info) const
{
	uint32_t const ray_direction_largest_dim = maths::MaximumDimension(maths::Abs(_ray.direction));
	uint32_t reordered_x = (ray_direction_largest_dim + 1u) % 3u;
	uint32_t reordered_y = (ray_direction_largest_dim + 2u) % 3u;
	if (_ray.direction[ray_direction_largest_dim] < 0._d)
	{
		std::swap(reordered_x, reordered_y);
	}

	InternalVec3 const reordered_direction{ maths::Swizzle(_ray.direction,
														   reordered_x,
														   reordered_y,
														   ray_direction_largest_dim) };
	InternalVec3 const shear_constants{ -reordered_direction.x / reordered_direction.z,
										-reordered_direction.y / reordered_direction.z,
										static_cast<InternalFloat>(1) / reordered_direction.z };

	std::array<InternalVec3, 3u> const local_vertices{
		maths::Swizzle(vertices[0] - static_cast<InternalVec3>(_ray.origin),
					   reordered_x,
					   reordered_y,
					   ray_direction_largest_dim),
		maths::Swizzle(vertices[1] - static_cast<InternalVec3>(_ray.origin),
					   reordered_x,
					   reordered_y,
					   ray_direction_largest_dim),
		maths::Swizzle(vertices[2] - static_cast<InternalVec3>(_ray.origin),
					   reordered_x,
					   reordered_y,
					   ray_direction_largest_dim)
	};

	maths::Vector<InternalVec3, 3> const transformed_vertices{
		local_vertices[0] * InternalVec3{ 1._d, 1._d, 0._d } + local_vertices[0].z * shear_constants,
		local_vertices[1] * InternalVec3{ 1._d, 1._d, 0._d } + local_vertices[1].z * shear_constants,
		local_vertices[2] * InternalVec3{ 1._d, 1._d, 0._d } + local_vertices[2].z * shear_constants
	};

	InternalVec3 const uvw{
		transformed_vertices[2].x * transformed_vertices[1].y -
		transformed_vertices[2].y * transformed_vertices[1].x,

		transformed_vertices[0].x * transformed_vertices[2].y -
		transformed_vertices[0].y * transformed_vertices[2].x,

		transformed_vertices[1].x * transformed_vertices[0].y -
		transformed_vertices[1].y * transformed_vertices[0].x
	};

	if (uvw.x < 0._d || uvw.y < 0._d || uvw.z < 0._d)
		return false;

	InternalFloat const determinant{ maths::FoldSum(uvw) };
	if (determinant == 0._d)
		return false;

	InternalFloat const scaled_hit_distance =
		maths::Dot(uvw,
				   InternalVec3{ transformed_vertices[0].z,
								 transformed_vertices[1].z,
								 transformed_vertices[2].z });

	if (scaled_hit_distance > (determinant * _ray.tMax) ||
		scaled_hit_distance < static_cast<InternalFloat>(0._d))
		return false;

	{
		InternalVec3 const barycentric{ uvw / determinant };
		InternalFloat const t = scaled_hit_distance / determinant;
		_tHit = static_cast<maths::Decimal>(t);

		InternalVec3 const barycentric_hit_point =
			barycentric[0] * vertices[0]
			+ barycentric[1] * vertices[1]
			+ barycentric[2] * vertices[2];

		InternalVec3 const t_hit_point =
			InternalVec3(_ray.origin) + static_cast<InternalVec3>(_ray.direction) * t;
		InternalVec3 const &ref_hit = t_hit_point;
		maths::Point3f const hit_point{ static_cast<maths::Decimal>(ref_hit.x),
										static_cast<maths::Decimal>(ref_hit.y),
										static_cast<maths::Decimal>(ref_hit.z) };

		maths::Vec3f error_bounds{ 0._d, 0._d, 0._d };
#if 0
		{
			maths::Vec3f const max_bounds{ std::max(std::abs(ref_hit.x.value - ref_hit.x.high_bound),
													std::abs(ref_hit.x.value - ref_hit.x.low_bound)),
										   std::max(std::abs(ref_hit.y.value - ref_hit.y.high_bound),
													std::abs(ref_hit.y.value - ref_hit.y.low_bound)),
										   std::max(std::abs(ref_hit.z.value - ref_hit.z.high_bound),
													std::abs(ref_hit.z.value - ref_hit.z.low_bound))
			};
			error_bounds = max_bounds;
		}
#endif
		{
			maths::Vec3f const gamma_bounds{ maths::Vec3f{ static_cast<maths::Decimal>(ref_hit.x),
														   static_cast<maths::Decimal>(ref_hit.y),
														   static_cast<maths::Decimal>(ref_hit.z) }
											 * maths::gamma(11u) };
			error_bounds = gamma_bounds;
		}

		maths::Point3f p0{ static_cast<maths::Decimal>(vertices[0].x),
						   static_cast<maths::Decimal>(vertices[0].y),
						   static_cast<maths::Decimal>(vertices[0].z) };
		maths::Point3f p1{ static_cast<maths::Decimal>(vertices[1].x),
						   static_cast<maths::Decimal>(vertices[1].y),
						   static_cast<maths::Decimal>(vertices[1].z) };
		maths::Point3f p2{ static_cast<maths::Decimal>(vertices[2].x),
						   static_cast<maths::Decimal>(vertices[2].y),
						   static_cast<maths::Decimal>(vertices[2].z) };
		maths::Norm3f const geometry_normal{
			static_cast<maths::Norm3f>(maths::Normalized(maths::Cross(p1 - p0, p2 - p0)))
		};
		raytracer::SurfaceInteraction::GeometryProperties const geometry{
			geometry_normal, maths::Vec3f(0._d), maths::Vec3f(0._d), maths::Norm3f(0._d), maths::Norm3f(0._d)
		};
		_hit_info = raytracer::SurfaceInteraction{
			hit_point, error_bounds, _tHit, -_ray.direction, nullptr, maths::Point2f(0._d), geometry, geometry
		};
	}

	return true;
}

template <typename InternalFloat>
bool
TriangleList<InternalFloat>:: Intersect(maths::Ray const &_ray,
										maths::Decimal &_tHit,
										raytracer::SurfaceInteraction &_hit_info) const
{
	bool hit = false;
	std::for_each(triangles.cbegin(), triangles.cend(), [&](Triangle_t const &_triangle) {
		hit = _triangle.Intersect(_ray, _tHit, _hit_info) || hit;
	});
	return hit;
}


} // namespace watertight_intersect

