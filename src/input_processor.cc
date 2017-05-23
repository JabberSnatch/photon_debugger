#include "input_processor.h"

/*
param_group ->		kString (kType kParamBegin kString* kParamEnd)+
translate_group ->	kTranslate kString kString kString
rotate_group ->		kRotate kString kString kString
scale_group ->		kScale kString kString kString
film_group ->		kFilm kScopeBegin param_group* kScopeEnd
camera_group ->		kCamera kScopeBegin param_group* kScopeEnd
shape_group ->		kShape kString kScopeBegin (param_group|translate_group|rotate_group|scale_group)* kScopeEnd
scene_group ->		(translate_group|rotate_group|scale_group|film_group|camera_group|shape_group)*
*/

/*
value_group ->		kString
value_group ->		kString kString
value_group ->		kString kString kString
value_group ->		kString kString kString kString

definition_group ->	kType kParamBegin value_group kParamEnd
definition_group ->	{}

param_group ->		kString definition_group

properties_group ->	param_group properties_group
properties_group -> {}

attribute_group ->	param_group attribute_group
attribute_group ->	translate_group attribute_group
attribute_group ->	rotate_group attribute_group
attribute_group ->	scale_group attribute_group
attribute_group ->	{}

translate_group ->	kTranslate kString kString kString
rotate_group ->		kRotate kString kString kString
scale_group ->		kScale kString kString kString
film_group ->		kFilm kScopeBegin properties_group kScopeEnd
camera_group ->		kCamera kScopeBegin properties_group kScopeEnd
shape_group ->		kShape kString kScopeBegin attribute_group kScopeEnd

scene_group ->		translate_group scene_group
scene_group ->		rotate_group scene_group
scene_group ->		scale_group scene_group
scene_group ->		film_group scene_group
scene_group ->		camera_group scene_group
scene_group ->		shape_group scene_group
scene_group ->		{}
*/

#include <vector>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <map>


namespace api {

enum Token
{
	kFilm,
	kCamera,
	kShape,
	kScopeBegin,
	kScopeEnd,
	kParamBegin,
	kParamEnd,
	kType,
	kTranslate,
	kRotate,
	kScale,
	kString,
	
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
	kSceneGroup,
};

std::unordered_map<std::string, Token> const		token_table
{
	{ "Film", kFilm },
	{ "Camera", kCamera },
	{ "Shape", kShape },
	{ "{", kScopeBegin },
	{ "}", kScopeEnd },
	{ "[", kParamBegin }, 
	{ "]", kParamEnd }, 
	{ "float", kType }, 
	{ "int", kType }, 
	{ "uint", kType },
	{ "Translate", kTranslate },
	{ "Rotate", kRotate },
	{ "Scale", kScale },
};

std::multimap<Token, std::vector<Token>> const		production_rules
{
	{ kValueGroup, { kString } },
	{ kValueGroup, { kString, kString } },
	{ kValueGroup, { kString, kString, kString } },
	{ kValueGroup, { kString, kString, kString, kString } },
	{ kDefinitionGroup, { kType, kParamBegin, kValueGroup, kParamEnd } },
	{ kDefinitionGroup, {} },
	{ kParamGroup, { kString, kDefinitionGroup } },
	{ kPropertiesGroup, { kParamGroup, kPropertiesGroup } },
	{ kPropertiesGroup, {} },
	{ kAttributeGroup, { kParamGroup, kAttributeGroup } },
	{ kAttributeGroup, { kTranslateGroup, kAttributeGroup } },
	{ kAttributeGroup, { kRotateGroup, kAttributeGroup } },
	{ kAttributeGroup, { kScaleGroup, kAttributeGroup } },
	{ kAttributeGroup, {} },
	{ kTranslateGroup, { kTranslate, kString, kString, kString } },
	{ kRotateGroup, { kRotate, kString, kString, kString } },
	{ kScaleGroup, { kScale, kString, kString, kString } },
	{ kFilmGroup, { kFilm, kScopeBegin, kPropertiesGroup, kScopeEnd } },
	{ kCameraGroup, { kCamera, kScopeBegin, kPropertiesGroup, kScopeEnd } },
	{ kShapeGroup, { kShape, kString, kScopeBegin, kAttributeGroup, kScopeEnd } },
	{ kSceneGroup, { kTranslateGroup, kSceneGroup } },
	{ kSceneGroup, { kRotateGroup, kSceneGroup } },
	{ kSceneGroup, { kScaleGroup, kSceneGroup } },
	{ kSceneGroup, { kFilmGroup, kSceneGroup } },
	{ kSceneGroup, { kCameraGroup, kSceneGroup } },
	{ kSceneGroup, { kShapeGroup, kSceneGroup } },
	{ kSceneGroup, {} },
};

struct SyntaxNode
{
	Token						token;
	std::vector<SyntaxNode>		children;
};

std::vector<Token>	LexicalAnalysis(std::ifstream &_file_stream,
									std::vector<std::string> &_string_values);
void				SyntaxAnalysis(std::vector<Token> const &_input_string,
								   std::vector<std::string> const &_string_values);
Token				StringToToken(std::string const &_string);


bool
ProcessInputFile(std::string const &_path)
{
	std::ifstream file_stream(_path);
	std::vector<std::string>	string_values;
	std::vector<Token>			input_string = LexicalAnalysis(file_stream, string_values);

	return false;
}


std::vector<Token>
LexicalAnalysis(std::ifstream &_file_stream,
				std::vector<std::string> &_string_values)
{
	std::vector<Token>	tokens;
	std::string			delimiters{ '\t', ' ', '\n' };

	for (;;)
	{
		std::stringstream	token_stream;
		
		while (std::find(delimiters.begin(), delimiters.end(), _file_stream.get())
			   != delimiters.end());
		_file_stream.unget();

		char buf;
		while (std::find(delimiters.begin(), delimiters.end(), buf = _file_stream.get())
			   == delimiters.end() && _file_stream.good())
			token_stream << buf;

		if (_file_stream.good())
		{
			Token	token = StringToToken(token_stream.str());
			tokens.push_back(token);
			if (token == kString)
				_string_values.push_back(token_stream.str());
		}
		else
			break;
	}

	return tokens;
}

void
SyntaxAnalysis(std::vector<Token> const &_input_string, 
			   std::vector<std::string> const &_string_values)
{
	SyntaxNode			root{ kSceneGroup };
	std::vector<Token>	planar_view{ kSceneGroup };
	size_t				input_index = 0;

	while (!std::equal(planar_view.begin(), planar_view.end(),
					   _input_string.begin(), _input_string.end()))
	{
		Token	token = planar_view[0];
		for (auto it = production_rules.find(token);
			 it != production_rules.end() && it->first == token;
			 ++it)
		{

		}
	}
}

Token
StringToToken(std::string const &_string)
{
	auto	it = token_table.find(_string);
	if (it != token_table.end())
		return it->second;
	return kString;
}

} // namespace api

