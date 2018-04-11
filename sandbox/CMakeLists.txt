
project(sampler_gnuplot)
add_executable(sampler_gnuplot sampler_gnuplot.cc)
target_link_libraries( sampler_gnuplot
					   PRIVATE
					   raytracer
					   )

project(test_ray-bounce)
add_executable(test_ray-bounce test_ray-bounce.cc)
target_link_libraries( test_ray-bounce
					   PRIVATE
					   raytracer
					   )

