
project(sampler_gnuplot)
add_executable(sampler_gnuplot sampler_gnuplot.cc)
target_link_libraries( sampler_gnuplot
					   PRIVATE
					   raytracer
					   )
