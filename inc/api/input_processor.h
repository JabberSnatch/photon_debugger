#pragma once
#ifndef __YS_INPUT_PARSER_HPP__
#define __YS_INPUT_PARSER_HPP__

#include <string>


namespace api {


bool	ProcessInputFile(std::string const &_path);


enum TokenId
{
	kNone,

	kOutput,
	kFilm,
	kCamera,
	kShape,
	kScopeBegin,
	kScopeEnd,
	kParamBegin,
	kParamEnd,
	kType,
	kTransformIdentity,
	kTranslate,
	kRotate,
	kScale,
	kString,
	kNumber,

	kValueGroup,
	kDefinitionGroup,
	kParamGroup,
	kPropertiesGroup,
	kAttributeGroup,
	kTranslateGroup,
	kRotateGroup,
	kScaleGroup,
	kFilmGroup,
	kCameraGroup,
	kShapeGroup,
	kOutputGroup,
	kSceneGroup,

	kEnd,
	kDefault,
};


struct Token
{
	TokenId		id;
	std::string	text;
};


} // namespace api


#endif // __YS_INPUT_PARSER_HPP__
