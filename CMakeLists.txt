cmake_minimum_required(VERSION 3.0)

set(SOURCE_DIR sources/src)
set(INCLUDE_DIR sources/inc)

set(LIB_DIR lib/)
set( LIB_INCLUDE_DIR
	 ${LIB_DIR}/boost/include/
	 ${LIB_DIR}/assimp/include/
	 ${LIB_DIR}/stb/include/
	 )

set(API_SOURCE_DIR ${SOURCE_DIR}/api)
set( API_SOURCES
	 ${API_SOURCE_DIR}/factory_functions.cc
	 ${API_SOURCE_DIR}/input_processor.cc
	 ${API_SOURCE_DIR}/param_set.cc
	 ${API_SOURCE_DIR}/render_context.cc
	 ${API_SOURCE_DIR}/resource_context.cc
	 ${API_SOURCE_DIR}/transform_cache.cc
	 ${API_SOURCE_DIR}/translation_state.cc
	 )

set(CORE_SOURCE_DIR ${SOURCE_DIR}/core)
set( CORE_SOURCES
	 ${CORE_SOURCE_DIR}/alloc.cc
	 ${CORE_SOURCE_DIR}/logger.cc
	 ${CORE_SOURCE_DIR}/memory_region.cc
	 ${CORE_SOURCE_DIR}/rng.cc
	 ${CORE_SOURCE_DIR}/win32_timer.cc
	 ${SOURCE_DIR}/globals.cc
	 )

set(MATHS_SOURCE_DIR ${SOURCE_DIR}/maths)
set( MATHS_SOURCES
	 ${MATHS_SOURCE_DIR}/maths.cc
	 ${MATHS_SOURCE_DIR}/quaternion.cc
	 ${MATHS_SOURCE_DIR}/redecimal.cc
	 ${MATHS_SOURCE_DIR}/transform.cc
	 )

set(RAYTRACER_SOURCE_DIR ${SOURCE_DIR}/raytracer)
set( RAYTRACER_SOURCES
	 ${RAYTRACER_SOURCE_DIR}/bvh_accelerator.cc
	 ${RAYTRACER_SOURCE_DIR}/camera.cc
	 ${RAYTRACER_SOURCE_DIR}/film.cc
	 ${RAYTRACER_SOURCE_DIR}/integrator.cc
	 ${RAYTRACER_SOURCE_DIR}/light.cc
	 ${RAYTRACER_SOURCE_DIR}/primitive.cc
	 ${RAYTRACER_SOURCE_DIR}/sampler.cc
	 ${RAYTRACER_SOURCE_DIR}/shape.cc
	 ${RAYTRACER_SOURCE_DIR}/surface_interaction.cc
	 ${RAYTRACER_SOURCE_DIR}/triangle_mesh_data.cc
	 ${RAYTRACER_SOURCE_DIR}/integrators/direct_lighting_integrator.cc
	 ${RAYTRACER_SOURCE_DIR}/samplers/halton_sampler.cc
	 ${RAYTRACER_SOURCE_DIR}/samplers/random_sampler.cc
	 ${RAYTRACER_SOURCE_DIR}/shapes/sphere.cc
	 ${RAYTRACER_SOURCE_DIR}/shapes/triangle.cc
	 ${RAYTRACER_SOURCE_DIR}/shapes/triangle_mesh.cc
	 )

set( PHDE_SOURCES
	 ${SOURCE_DIR}/main.cc
	 )



# PROJECTS DECLARATIONS
project(core)
project(maths)
project(raytracer)
project(api)
project(PhotonDebugger)


# COMPILER CONFIGURATION
include_directories(${INCLUDE_DIR} ${LIB_INCLUDE_DIR})
link_directories(${LIB_DIR}/boost/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")

add_definitions(-DNOMINMAX -D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS)

# boost compile error with msvc 2017 workaround
add_definitions(-D_HAS_AUTO_PTR_ETC)

# EXTERNAL LIBRARIES
find_library(ASSIMP_LIB assimp PATHS ${LIB_DIR}/assimp/lib64/)


# PROJECTS DEFINITIONS
# CORE
add_library(core STATIC ${CORE_SOURCES})

# MATHS
add_library(maths STATIC ${MATHS_SOURCES})

# RAYTRACER
add_library(raytracer STATIC ${RAYTRACER_SOURCES})
target_link_libraries( raytracer
					   PUBLIC
					   core
					   maths
					   )

# API
add_library(api STATIC ${API_SOURCES})
target_link_libraries( api
					   PUBLIC
					   raytracer
					   PRIVATE
					   "${ASSIMP_LIB}"
					   )

# PHOTONDEBUGGER
add_executable(PhotonDebugger ${PHDE_SOURCES})
target_link_libraries( PhotonDebugger
					   PRIVATE
					   api
					   )



# SUBDIRECTORIES
add_subdirectory(sandbox)
