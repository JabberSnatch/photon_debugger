
option(WATERTIGHT_PRINT_ROUND_COUNTS "Activates round counts printing in test_ray-bounce" OFF)

# ==============================================================================

project(sampler_gnuplot)
add_executable(sampler_gnuplot sampler_gnuplot.cc)
target_link_libraries( sampler_gnuplot
					   PRIVATE
					   raytracer
					   )

# ==============================================================================

project(test_ray-bounce)

if (WATERTIGHT_PRINT_ROUND_COUNTS)
add_definitions(-DWATERTIGHT_PRINT_ROUND_COUNTS)
endif(WATERTIGHT_PRINT_ROUND_COUNTS)

add_executable(test_ray-bounce test_ray-bounce.cc)
target_link_libraries( test_ray-bounce
					   PRIVATE
					   raytracer
					   )

