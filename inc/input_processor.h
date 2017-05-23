#pragma once
#ifndef __YS_INPUT_PARSER_HPP__
#define __YS_INPUT_PARSER_HPP__

#include <string>

namespace api {

// Camera, Film, Shape
// Scene : list of shapes
// Shape -> Transform + flip_normals + individual properties
// For instance : Sphere{Transform, flip_normals, radius, min_z, max_z, max_phi}
// Each parameter can have default value.
// Specifying that we just want a sphere would be fine.
// Values should be assigned explicitly
//
// scope_id [ value ] et ref[ id ] last
//
// Film scope_id [ default_film ] { resolution [ 1920 1080 ] }
// Camera {
// position float [ 0 -5 5 ] lookat float [ 0 1 .5 ]
// up float [ 0 0 1 ]
// fov float [ 60 ]
// film ref[ default_film ]
// } 
// Translate 10 0 0
// Rotate 30 0 0
// Shape sphere {
// radius float [ 5 ]
// min_z float [ -10 ]
// max_z float [ 10 ]
// max_phi float [ 360 ]
// flip_normals bool [ false ]
// Translate 5 0 0
// Rotate 25 0 25
// Scale 3 3 1
// }
// Shape sphere {
// radius float [ 5 ]
// min_z float [ -3 ]
// max_z float [ 4 ]
// max_phi float [ 360 ]
// }

bool ProcessInputFile(std::string const &_path);

} // namespace api


#endif // __YS_INPUT_PARSER_HPP__
