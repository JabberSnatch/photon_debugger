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
rotate_group ->		kRotate kNumber kNumber kNumber kNumber
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
#include <functional>

#include "common_macros.h"
#include "vector.h"
#include "transform.h"
#include "quaternion.h"
#include "camera.h"
#include "memory_region.h"


namespace api {


struct ProductionMetadata
{
	TokenId		parent;
	uint32_t	unscanned_count;
	std::vector<Token>::const_iterator	begin;
};

using TokenTable_t = std::unordered_map<std::string, TokenId>;
using Production_t = std::vector<TokenId>;
using TokenProductions_t = std::map<TokenId, Production_t>;
using ProductionRules_t = std::map<TokenId, TokenProductions_t>;
using ProductionStack_t = std::vector<ProductionMetadata>;

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
	{ "bool", kType },
	{ "false", kNumber },
	{ "true", kNumber },
	{ "Identity", kTransformIdentity },
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
		{ kDefault, { kRotate, kNumber, kNumber, kNumber, kNumber } },
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

using TranslationCallback_t = std::function<
	void (TranslationState&, std::vector<Token>::const_iterator, std::vector<Token>::const_iterator)
>;
using TranslationTable_t = std::unordered_map<TokenId, TranslationCallback_t>;

void	ParamGroup(TranslationState &_state, 
				   std::vector<Token>::const_iterator _production_begin,
				   std::vector<Token>::const_iterator _production_end);
void	ShapeGroup(TranslationState &_state, 
				   std::vector<Token>::const_iterator _production_begin,
				   std::vector<Token>::const_iterator _production_end);
void	FilmGroup(TranslationState &_state,
				  std::vector<Token>::const_iterator _production_begin,
				  std::vector<Token>::const_iterator _production_end);
void	CameraGroup(TranslationState &_state,
					std::vector<Token>::const_iterator _production_begin,
					std::vector<Token>::const_iterator _production_end);
void	TranslateGroup(TranslationState &_state,
					   std::vector<Token>::const_iterator _production_begin,
					   std::vector<Token>::const_iterator _production_end);
void	RotateGroup(TranslationState &_state,
					std::vector<Token>::const_iterator _production_begin,
					std::vector<Token>::const_iterator _production_end);
void	ScaleGroup(TranslationState &_state,
				   std::vector<Token>::const_iterator _production_begin,
				   std::vector<Token>::const_iterator _production_end);

void	IdentityTerminal(TranslationState &_state,
						 std::vector<Token>::const_iterator _production_begin,
						 std::vector<Token>::const_iterator _production_end);
void	ScopeBeginTerminal(TranslationState &_state,
						   std::vector<Token>::const_iterator _production_begin,
						   std::vector<Token>::const_iterator _production_end);
void	ScopeEndTerminal(TranslationState &_state,
						 std::vector<Token>::const_iterator _production_begin,
						 std::vector<Token>::const_iterator _production_end);

TranslationTable_t	semantic_actions = {
	{ kParamGroup, &api::ParamGroup },
	{ kShapeGroup, &api::ShapeGroup },
	{ kFilmGroup, &api::FilmGroup },
	{ kCameraGroup, &api::CameraGroup },
	{ kTranslateGroup, &api::TranslateGroup },
	{ kRotateGroup, &api::RotateGroup },
	{ kScaleGroup, &api::ScaleGroup },
	{ kTransformIdentity, &api::IdentityTerminal },
	{ kScopeBegin, &api::ScopeBeginTerminal },
	{ kScopeEnd, &api::ScopeEndTerminal },
};


std::vector<Token>	LexicalAnalysis(std::ifstream &_file_stream);
bool				SyntaxAnalysis(std::vector<Token> const &_input_string);
Token				StringToToken(std::string const &_string);



bool
ProcessInputFile(std::string const &_path)
{
	std::ifstream file_stream(_path);
	std::vector<Token>			input_string = LexicalAnalysis(file_stream);
	if (SyntaxAnalysis(input_string))
	{
		int i = 0;
	}

	return false;
}


std::vector<Token>
LexicalAnalysis(std::ifstream &_file_stream)
{
	std::vector<Token>	tokens;
	std::string const	delimiters{ '\t', ' ' , '\n' };

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
		}
		else
			break;
	}

	tokens.push_back(Token{ kEnd, "" });
	return tokens;
}

bool
SyntaxAnalysis(std::vector<Token> const &_input_string)
{
	YS_ASSERT(_input_string.size() > 0);
	bool				input_is_valid = true;

	TranslationState	state{};
	ProductionStack_t	production_stack;

	std::vector<TokenId>	planar_view{ kSceneGroup };
	size_t					input_index = 0;
	size_t					value_index = 0;
	size_t					view_index = 0;


	while (input_index < _input_string.size())
	{
		TokenId const	planar_token = planar_view[view_index];
		Token const		lookahead = _input_string[input_index];

		auto const		token_productions_it = production_rules.find(planar_token);
		if (token_productions_it != production_rules.end())
		{
			TokenProductions_t const &available_productions = token_productions_it->second;

			auto	production_it = available_productions.find(lookahead.id);
			if (production_it == available_productions.end())
				production_it = available_productions.find(kDefault);

			if (production_it != available_productions.end())
			{
				Production_t const &production = production_it->second;
				planar_view.insert(planar_view.cbegin() + view_index + 1,
								   production.begin(), production.end());
				planar_view.erase(planar_view.cbegin() + view_index);

				uint32_t const		production_size = static_cast<uint32_t>(production.size());
				std::vector<Token>::const_iterator const	production_begin =
					_input_string.begin() + input_index;

				production_stack.push_back({ planar_token, production_size, production_begin });
			}
			else
			{
				// Production not found for current lookahead => syntax error, unexpected token
				input_is_valid = false;
			}
		}
		else
		{
			// No production rules found for current token => terminal token
			if (lookahead.id == planar_token)
			{
				if (production_stack.back().unscanned_count != 0)
				{
					view_index++;
					input_index++;

					auto	action_it = semantic_actions.find(lookahead.id);
					if (action_it != semantic_actions.end())
						action_it->second(state,
										  _input_string.begin() + input_index - 1,
										  _input_string.begin() + input_index);
				}
				else
				{
					// An empty production should be ignored and there can't be more than one in a row.
					production_stack.pop_back();
				}
				production_stack.back().unscanned_count--;
			}
			else
			{
				// Token mismatch
				input_is_valid = false;
			}
		}

		while (production_stack.back().unscanned_count == 0)
		{
			auto	action_it = semantic_actions.find(production_stack.back().parent);
			if (action_it != semantic_actions.end())
				action_it->second(state, 
								  production_stack.back().begin, 
								  _input_string.begin() + input_index);

			production_stack.pop_back();
			if (production_stack.size() == 0)
				break;

			production_stack.back().unscanned_count--;
		}
	}

	return input_is_valid;
}

Token
StringToToken(std::string const &_string)
{
	YS_ASSERT(_string.size() > 0);

	TokenId		id = TokenId::kNone;
	auto	it = token_table.find(_string);
	if (it != token_table.end())
		id = it->second;
	else
	{
		if ((_string[0] >= '0' && _string[0] <= '9') ||
			(_string[0] == '-' || _string[0] == '+' || _string[0] == '.'))
			id = kNumber;
		else
			id = kString;
	}

	return Token{ id, _string };
}


void
ParamGroup(TranslationState &_state, 
		   std::vector<Token>::const_iterator _production_begin, 
		   std::vector<Token>::const_iterator _production_end)
{
	YS_ASSERT(_production_begin->id == api::kString);

	std::string const	identifier = _production_begin->text;

	if (_production_end[-1].id != api::kString)
	{
		uint32_t	value_count = static_cast<uint32_t>(
			std::count_if(_production_begin, _production_end,
						  [](Token const &_t) -> bool { return _t.id == kNumber; }));

		auto		cursor = 
			std::find_if(_production_begin, _production_end,
						 [](Token const &_t) -> bool { return _t.id == kParamBegin; });

		std::string	type_string = cursor[-1].text;
		cursor++;

		if (type_string == "float")
		{
			maths::Decimal	*values = new maths::Decimal[value_count];
			for (uint32_t i = 0; i < value_count; ++i)
				values[i] = maths::stof(cursor[i].text);
			_state.param_set().PushFloat(identifier, values, value_count);
		}
		else if (type_string == "int")
		{
			int32_t			*values = new int32_t[value_count];
			for (uint32_t i = 0; i < value_count; ++i)
				values[i] = std::stoi(cursor[i].text);
			_state.param_set().PushInt(identifier, values, value_count);
		}
		else if (type_string == "bool")
		{
			bool			*values = new bool[value_count];
			for (uint32_t i = 0; i < value_count; ++i)
			{
				values[i] = (cursor[i].text == "true") ? true : false;
				if (cursor[i].text != "true" && cursor[i].text != "false")
					LOG_WARNING(tools::kChannelGeneral, "Illegal expression for boolean value. Expected true or false.");
			}
			_state.param_set().PushBool(identifier, values, value_count);
		}
	}
	else
	{
		// In case we find a lone identifier, we treat it as a bool set to true
		_state.param_set().PushBool(identifier, true);
	}
}

void
ShapeGroup(TranslationState &_state,
		   std::vector<Token>::const_iterator _production_begin,
		   std::vector<Token>::const_iterator _production_end)
{}
void
FilmGroup(TranslationState &_state,
		  std::vector<Token>::const_iterator _production_begin,
		  std::vector<Token>::const_iterator _production_end)
{}
void
CameraGroup(TranslationState &_state,
			std::vector<Token>::const_iterator _production_begin,
			std::vector<Token>::const_iterator _production_end)
{}
void
TranslateGroup(TranslationState &_state,
			   std::vector<Token>::const_iterator _production_begin,
			   std::vector<Token>::const_iterator _production_end)
{
	std::vector<Token>::const_iterator const	it = _production_begin + 1;
	_state.Translate({ maths::stof(it->text),
					   maths::stof((it + 1)->text),
					   maths::stof((it + 2)->text) });
}
void
RotateGroup(TranslationState &_state,
			std::vector<Token>::const_iterator _production_begin,
			std::vector<Token>::const_iterator _production_end)
{
	std::vector<Token>::const_iterator const	it = _production_begin + 1;
	_state.Rotate(maths::stof(it->text),
				  { maths::stof((it + 1)->text),
					maths::stof((it + 2)->text),
					maths::stof((it + 3)->text) });
}
void
ScaleGroup(TranslationState &_state,
		   std::vector<Token>::const_iterator _production_begin,
		   std::vector<Token>::const_iterator _production_end)
{
	std::vector<Token>::const_iterator const	it = _production_begin + 1;
	_state.Scale(maths::stof(it->text),
				 maths::stof((it + 1)->text),
				 maths::stof((it + 2)->text));
}
void
IdentityTerminal(TranslationState &_state,
				 std::vector<Token>::const_iterator _production_begin,
				 std::vector<Token>::const_iterator _production_end)
{
	_state.Identity();
}
void
ScopeBeginTerminal(TranslationState &_state,
				   std::vector<Token>::const_iterator _production_begin,
				   std::vector<Token>::const_iterator _production_end)
{
	_state.ScopeBegin();
}
void
ScopeEndTerminal(TranslationState &_state,
				 std::vector<Token>::const_iterator _production_begin,
				 std::vector<Token>::const_iterator _production_end)
{
	_state.ScopeEnd();
}


TranslationState::TranslationState() :
	render_context_{}, parameters_{},
	scope_depth_{ 1 }, transform_stack_{ maths::Transform{} }
{}
void
TranslationState::Identity()
{
	transform_stack_.back() = maths::Transform{};
}
void
TranslationState::Translate(maths::Vec3f const &_t)
{
	transform_stack_.back() = transform_stack_.back() * maths::Translate(_t);
}
void
TranslationState::Rotate(maths::Decimal _angle, maths::Vec3f const &_axis)
{
	transform_stack_.back() = transform_stack_.back() * maths::Rotate(_angle, _axis);
}
void
TranslationState::Scale(maths::Decimal _x, maths::Decimal _y, maths::Decimal _z)
{
	transform_stack_.back() = transform_stack_.back() * maths::Scale(_x, _y, _z);
}
void
TranslationState::ScopeBegin()
{
	transform_stack_.push_back(transform_stack_.back());
	scope_depth_++;
}
void
TranslationState::ScopeEnd()
{
	scope_depth_--;
	transform_stack_.pop_back();
	parameters_.Clear();
}

} // namespace api

