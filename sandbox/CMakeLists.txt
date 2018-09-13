
option(WATERTIGHT_PRINT_ROUND_COUNTS "Activates round counts printing in ray-bounce" OFF)

# ==============================================================================

project(sampler_gnuplot)
add_executable(sampler_gnuplot sampler_gnuplot.cc)
target_link_libraries( sampler_gnuplot
					   PRIVATE
					   raytracer
					   )

# ==============================================================================

project(ray-bounce)

if (WATERTIGHT_PRINT_ROUND_COUNTS)
add_definitions(-DWATERTIGHT_PRINT_ROUND_COUNTS)
endif(WATERTIGHT_PRINT_ROUND_COUNTS)

add_executable(ray-bounce ray-bounce/main_test.cc)
target_link_libraries( ray-bounce
					   PRIVATE
					   raytracer
					   )

