#include "input_processor.h"

#include "common_macros.h"

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
value_group ->		kNumber
value_group ->		kNumber kNumber
value_group ->		kNumber kNumber kNumber
value_group ->		kNumber kNumber kNumber kNumber

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

translate_group ->	kTranslate kNumber kNumber kNumber 
rotate_group ->		kRotate kNumber kNumber kNumber 
scale_group ->		kScale kNumber kNumber kNumber 

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

Lookahead identifier deductions

kTranslate => translate_group
kRotate => rotate_group
kScale => scale_group
kFilm => film_group
kCamera => camera_group
kShape => shape_group

kString => param_group
kNumber => value_group

param_group => properties_group / attribute_group

translate_group => kTranslate
rotate_group => kRotate
scale_group => kScale
film_group => 
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
	kSceneGroup,
	
	kEnd,
	kDefault,
};

using TokenTable_t = std::unordered_map<std::string, Token>;
using Production_t = std::vector<Token>;
using TokenProductions_t = std::map<Token, Production_t>;
using ProductionRules_t = std::map<Token, TokenProductions_t>;

TokenTable_t const		token_table
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


ProductionRules_t const		production_rules
{
	{ kSceneGroup, {
		{ kTranslate, { kTranslateGroup, kSceneGroup } },
		{ kRotate, { kRotateGroup, kSceneGroup } },
		{ kScale, { kScaleGroup, kSceneGroup } },
		{ kFilm, { kFilmGroup, kSceneGroup } },
		{ kCamera, { kCameraGroup, kSceneGroup } },
		{ kShape, { kShapeGroup, kSceneGroup } },
		{ kDefault, { kEnd } }
	} },

	{ kTranslateGroup, {
		{ kDefault, { kTranslate, kNumber, kNumber, kNumber } },
	} },

	{ kRotateGroup, {
		{ kDefault, { kRotate, kNumber, kNumber, kNumber } },
	} },

	{ kScaleGroup, {
		{ kDefault, { kScale, kNumber, kNumber, kNumber } },
	} },

	{ kFilmGroup, {
		{ kDefault, { kFilm, kScopeBegin, kPropertiesGroup, kScopeEnd } },
	} },

	{ kCameraGroup, {
		{ kDefault, { kCamera, kScopeBegin, kPropertiesGroup, kScopeEnd } },
	} },

	{ kShapeGroup, {
		{ kDefault, { kShape, kString, kScopeBegin, kAttributeGroup, kScopeEnd } },
	} },

	{ kPropertiesGroup, {
		{ kString, { kParamGroup, kPropertiesGroup } },
		{ kDefault, {} }
	} },

	{ kParamGroup, {
		{ kDefault, { kString, kDefinitionGroup } },
	} },

	{ kDefinitionGroup, {
		{ kType, { kType, kParamBegin, kValueGroup, kParamEnd } },
		{ kDefault, {} },
	} },

	{ kValueGroup, {
		{ kNumber, { kNumber, kValueGroup } },
		{ kDefault, {} } ,
	} },
	
	{ kAttributeGroup, {
		{ kString, { kParamGroup, kAttributeGroup } },
		{ kTranslate, { kTranslateGroup, kAttributeGroup } },
		{ kRotate, { kRotateGroup, kAttributeGroup } },
		{ kScale, { kScaleGroup, kAttributeGroup } },
		{ kDefault, {} },
	} },
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
	SyntaxAnalysis(input_string, string_values);

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
			if (token == kString || token == kNumber)
				_string_values.push_back(token_stream.str());
		}
		else
			break;
	}

	tokens.push_back(kEnd);
	return tokens;
}

void
SyntaxAnalysis(std::vector<Token> const &_input_string, 
			   std::vector<std::string> const &_string_values)
{
	YS_ASSERT(_input_string.size() > 0);

	std::vector<Token>	planar_view{ kSceneGroup };
	size_t				input_index = 0;
	size_t				value_index = 0;
	size_t				view_index = 0;

	while (input_index < _input_string.size())
	{
		Token	planar_token = planar_view[view_index];
		Token	lookahead = _input_string[input_index];

		auto	token_productions_it = production_rules.find(planar_token);
		if (token_productions_it != production_rules.end())
		{
			TokenProductions_t const &available_productions = token_productions_it->second;
			auto	production_it = available_productions.find(lookahead);
			if (production_it == available_productions.end())
				production_it = available_productions.find(kDefault);

			if (production_it != available_productions.end())
			{
				Production_t const &production = production_it->second;
				//input_index += production.size();
				planar_view.insert(planar_view.cbegin() + view_index + 1,
								   production.begin(), production.end());
				planar_view.erase(planar_view.cbegin() + view_index);
			}
			else
			{
				// Production not found for current lookahead => syntax error, unexpected token
				int i = 0;
			}
		}
		else
		{
			// No production rules found for current token => terminal token
			if (lookahead == planar_token)
			{
				view_index++;
				input_index++;
			}
			else
			{
				// Token mismatch
				int i = 0;
			}
		}
	}

	int k = 0;
}

Token
StringToToken(std::string const &_string)
{
	YS_ASSERT(_string.size() > 0);

	auto	it = token_table.find(_string);
	if (it != token_table.end())
		return it->second;
	else
	{
		if ((_string[0] >= '0' && _string[0] <= '9') ||
			(_string[0] == '-' || _string[0] == '+' || _string[0] == '.'))
			return kNumber;
	}
	return kString;
}

} // namespace api

