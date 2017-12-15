#include "api/input_processor.h"

#include <vector>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <functional>

#include "boost/lexical_cast.hpp"
#include "boost/filesystem.hpp"

#include "api/factory_functions.h"
#include "api/translation_state.h"
#include "common_macros.h"
#include "core/logger.h"
#include "core/memory_region.h"
#include "maths/quaternion.h"
#include "maths/transform.h"
#include "maths/vector.h"
#include "raytracer/camera.h"
#include "raytracer/film.h"
#include "raytracer/primitive.h"
#include "raytracer/shapes/sphere.h"


namespace api {



enum TokenId
{
	kNone,

	kOutput,
	kObjectId,
	kFilm,
	kCamera,
	kShape,
	kLight,
	kSampler,
	kIntegrator,
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
	kObjectIdGroup,
	kFilmGroup,
	kCameraGroup,
	kShapeGroup,
	kLightGroup,
	kSamplerGroup,
	kIntegratorGroup,
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
	{ "Output", kOutput },
	{ "ID", kObjectId },
	{ "Film", kFilm },
	{ "Camera", kCamera },
	{ "Shape", kShape },
	{ "Light", kLight },
	{ "Sampler", kSampler },
	{ "Integrator", kIntegrator },
	{ "{", kScopeBegin },
	{ "}", kScopeEnd },
	{ "[", kParamBegin }, 
	{ "]", kParamEnd }, 
	{ "float", kType }, 
	{ "int", kType },
	{ "uint", kType },
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
		{ kTransformIdentity, { kTransformIdentity, kSceneGroup } },
		{ kOutput, { kOutputGroup, kSceneGroup } },

		{ kFilm, { kFilmGroup, kSceneGroup } },
		{ kCamera, { kCameraGroup, kSceneGroup } },
		{ kShape, { kShapeGroup, kSceneGroup } },
		{ kLight, { kLightGroup, kSceneGroup } },
		{ kSampler, { kSamplerGroup, kSceneGroup } },
		{ kIntegrator, { kIntegratorGroup, kSceneGroup } },
		{ kDefault, { kEnd } }
	} },

	{ kOutputGroup, {
		{ kDefault, { kOutput, kString } },
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

	{ kObjectIdGroup, {
		{ kDefault, { kObjectId, kString } },
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
	{ kLightGroup, {
		{ kDefault, { kLight, kString, kScopeBegin, kAttributeGroup, kScopeEnd } },
	} },
	{ kSamplerGroup, {
		{ kDefault, { kSampler, kString, kScopeBegin, kAttributeGroup, kScopeEnd } },
	} },
	{ kIntegratorGroup,{
		{ kDefault, { kIntegrator, kString, kScopeBegin, kAttributeGroup, kScopeEnd } },
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
		{ kString, { kString } },
		{ kDefault, {} },
	} },

	{ kValueGroup, {
		{ kNumber, { kNumber, kValueGroup } },
		{ kDefault, {} } ,
	} },
	
	{ kAttributeGroup, {
		{ kTranslate, { kTranslateGroup, kAttributeGroup } },
		{ kRotate, { kRotateGroup, kAttributeGroup } },
		{ kScale, { kScaleGroup, kAttributeGroup } },
		{ kTransformIdentity, { kTransformIdentity, kAttributeGroup } },
		{ kObjectId, { kObjectIdGroup, kAttributeGroup } },

		{ kString, { kParamGroup, kAttributeGroup } },
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
void	ObjectIdGroup(TranslationState &_state,
					  std::vector<Token>::const_iterator _production_begin,
					  std::vector<Token>::const_iterator _production_end);
void	ShapeGroup(TranslationState &_state, 
				   std::vector<Token>::const_iterator _production_begin,
				   std::vector<Token>::const_iterator _production_end);
void	LightGroup(TranslationState &_state,
				   std::vector<Token>::const_iterator _production_begin,
				   std::vector<Token>::const_iterator _production_end);
void	SamplerGroup(TranslationState &_state,
					 std::vector<Token>::const_iterator _production_begin,
					 std::vector<Token>::const_iterator _production_end);
void	IntegratorGroup(TranslationState &_state,
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
void	OutputGroup(TranslationState &_state,
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
	{ kObjectIdGroup, &api::ObjectIdGroup },
	{ kShapeGroup, &api::ShapeGroup },
	{ kLightGroup, &api::LightGroup },
	{ kSamplerGroup, &api::SamplerGroup },
	{ kIntegratorGroup, &api::IntegratorGroup },
	{ kFilmGroup, &api::FilmGroup },
	{ kCameraGroup, &api::CameraGroup },
	{ kTranslateGroup, &api::TranslateGroup },
	{ kRotateGroup, &api::RotateGroup },
	{ kScaleGroup, &api::ScaleGroup },
	{ kOutputGroup, &api::OutputGroup },

	{ kTransformIdentity, &api::IdentityTerminal },
	{ kScopeBegin, &api::ScopeBeginTerminal },
	{ kScopeEnd, &api::ScopeEndTerminal },
};


std::vector<Token>	LexicalAnalysis(std::ifstream &_file_stream);
bool				SyntaxAnalysis(std::vector<Token> const &_input_string,
								   api::TranslationState &_state);
Token				StringToToken(std::string const &_string);



bool
ProcessInputFile(std::string const &_path, TranslationState &_state)
{
	std::ifstream file_stream(_path);
	//
	std::vector<Token> const input_string = LexicalAnalysis(file_stream);
	//
	boost::filesystem::path const input_path{ boost::filesystem::absolute(_path) };
	_state.Workdir(input_path.parent_path().string());
	_state.Output(input_path.stem().string() + ".png");
	//
	bool const result = SyntaxAnalysis(input_string, _state);
	return result;
}


std::vector<Token>
LexicalAnalysis(std::ifstream &_file_stream)
{
	std::vector<Token>	tokens;
	std::string const	delimiters{ '\t', ' ' , '\n' };

	for (;;)
	{
		std::stringstream	token_stream;
		bool ignore_line = false;
		
		while (std::find(delimiters.begin(), delimiters.end(), _file_stream.get())
			   != delimiters.end());
		_file_stream.unget();

		char buf;
		while (std::find(delimiters.begin(), delimiters.end(), buf = _file_stream.get())
			   == delimiters.end() && _file_stream.good())
		{
			token_stream << buf;
			if (token_stream.str() == "//")
			{
				ignore_line = true;
				while (_file_stream.get() != '\n');
				break;
			}
		}

		if (ignore_line)
		{
		}
		else if (_file_stream.good())
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
SyntaxAnalysis(std::vector<Token> const &_input_string, TranslationState &_state)
{
	YS_ASSERT(_input_string.size() > 0);
	bool					input_is_valid = true;
	ProductionStack_t		production_stack;
	std::vector<TokenId>	planar_view{ kSceneGroup };
	size_t					input_index = 0;
	size_t					value_index = 0;
	size_t					view_index = 0;
	//
	_state.SceneBegin();
	while (input_index < _input_string.size())
	{
		TokenId const	planar_token = planar_view[view_index];
		Token const		lookahead = _input_string[input_index];
		auto const		token_productions_it = production_rules.find(planar_token);
		if (token_productions_it != production_rules.end())
		{
			// NOTE: Token expansion
			TokenProductions_t const &available_productions = token_productions_it->second;
			auto	production_it = available_productions.find(lookahead.id);
			if (production_it == available_productions.end())
				production_it = available_productions.find(kDefault);
			//
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
				//
				LOG_INFO(tools::kChannelParsing,
						 "Expanding token " + 
						 boost::lexical_cast<std::string>(planar_token) + 
						 " using lookahead " +
						 boost::lexical_cast<std::string>(lookahead.id)
				);
			}
			else
			{
				// Production not found for current lookahead => syntax error, unexpected token
				LOG_ERROR(tools::kChannelParsing,
						  "No production found for token " +
						  boost::lexical_cast<std::string>(planar_token) +
						  " matching lookahead " +
						  boost::lexical_cast<std::string>(lookahead.id)
				);
				input_is_valid = false;
			}
		}
		else
		{
			// NOTE: Syntax validation + Terminal semantic actions
			if (lookahead.id == planar_token)
			{
				LOG_INFO(tools::kChannelParsing,
						 "Successfully matched terminal token " +
						 boost::lexical_cast<std::string>(planar_token)
				);
				//
				if (production_stack.back().unscanned_count != 0)
				{
					TokenId	const	past_token = (input_index > 0) ? 
						_input_string[input_index - 1].id : kEnd;
					view_index++;
					input_index++;
					auto	action_it = semantic_actions.find(past_token);
					if (action_it != semantic_actions.end())
					{
						LOG_INFO(tools::kChannelParsing, "Found semantic action ");
						action_it->second(_state,
										  _input_string.begin() + input_index - 2,
										  _input_string.begin() + input_index - 1);
					}
					LOG_INFO(tools::kChannelParsing, "Stepped forward");
				}
				else
				{
					// An empty production should be ignored and there can't be more than one in a row.
					production_stack.pop_back();
					LOG_INFO(tools::kChannelParsing, "Popped empty production");
				}
				production_stack.back().unscanned_count--;
			}
			else
			{
				// Token mismatch
				LOG_INFO(tools::kChannelParsing,
						 "Token mismatch, expected " +
						 boost::lexical_cast<std::string>(planar_token) +
						 " but found " +
						 boost::lexical_cast<std::string>(lookahead.id) +
						 " instead"
				);
				input_index++;
				input_is_valid = false;
			}
		}

		// NOTE: Non-terminal semantic actions
		while (production_stack.back().unscanned_count == 0)
		{
			LOG_INFO(tools::kChannelParsing,
					 "Popping nonterminal " +
					 boost::lexical_cast<std::string>(production_stack.back().parent) +
					 " from production stack"
			);
			auto	action_it = semantic_actions.find(production_stack.back().parent);
			if (action_it != semantic_actions.end())
			{
				LOG_INFO(tools::kChannelParsing, "Found semantic action ");
				action_it->second(_state,
								  production_stack.back().begin, 
								  _input_string.begin() + input_index);
			}
			production_stack.pop_back();
			if (production_stack.size() == 0)
				break;
			production_stack.back().unscanned_count--;
		}
		//
		LOG_INFO(tools::kChannelParsing, "");
	}

	if (input_is_valid)
		_state.SceneEnd();

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
	LOG_INFO(tools::kChannelParsing, "Parameter group ended, applying semantic action..");

	std::string const	identifier = _production_begin->text;

	//if (_production_end[-1].id != api::kString)
	int64_t const production_size = 
		static_cast<int64_t>(std::distance(_production_begin, _production_end));
	if (production_size > 1)
	{
		if (_production_begin[1].id == kType)
		{
			uint32_t	value_count = static_cast<uint32_t>(
				std::count_if(_production_begin, _production_end,
							  [](Token const &_t) -> bool { return _t.id == kNumber; }));

			auto		cursor = 
				std::find_if(_production_begin, _production_end,
							 [](Token const &_t) -> bool { return _t.id == kParamBegin; });

			std::string	const type_string = cursor[-1].text;
			cursor++;

			if (type_string == "float")
			{
				maths::Decimal	*values = new maths::Decimal[value_count];
				for (uint32_t i = 0; i < value_count; ++i)
					values[i] = boost::lexical_cast<maths::Decimal>(cursor[i].text);
				_state.param_set().PushFloat(identifier, values, value_count);
			}
			else if (type_string == "int")
			{
				int64_t			*values = new int64_t[value_count];
				for (uint32_t i = 0; i < value_count; ++i)
					values[i] = boost::lexical_cast<int64_t>(cursor[i].text);
				_state.param_set().PushInt(identifier, values, value_count);
			}
			else if (type_string == "uint")
			{
				uint64_t		*values = new uint64_t[value_count];
				for (uint32_t i = 0; i < value_count; ++i)
					values[i] = boost::lexical_cast<uint64_t>(cursor[i].text);
				_state.param_set().PushUint(identifier, values, value_count);
			}
			else if (type_string == "bool")
			{
				bool			*values = new bool[value_count];
				for (uint32_t i = 0; i < value_count; ++i)
				{
					values[i] = (cursor[i].text == "true") ? true : false;
					if (cursor[i].text != "true" && cursor[i].text != "false")
					{
						LOG_WARNING(tools::kChannelGeneral, "Illegal expression for boolean value. Expected true or		false.");
					}
				}
				_state.param_set().PushBool(identifier, values, value_count);
			}
		} // if (_production_begin[1].id == kType)
		else
		{
			YS_ASSERT(_production_begin[1].id == kString);
			YS_ASSERT(production_size == 2);
			std::string const value = _production_begin[1].text;
			_state.param_set().PushString(identifier, value);
		}
	} // if (production_size > 1)
	else
	{
		// If we find a lone identifier, we treat it as a bool set to true
		_state.param_set().PushBool(identifier, true);
	}
}

void
ObjectIdGroup(TranslationState &_state,
			  std::vector<Token>::const_iterator _production_begin,
			  std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "ObjectId group ended, applying semantic action..");
	std::string const	object_id = std::next(_production_begin, 1)->text;
	_state.ObjectId(object_id);
}

void
ShapeGroup(TranslationState &_state,
		   std::vector<Token>::const_iterator _production_begin,
		   std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "Shape group ended, applying semantic action..");
	std::string const	shape_type = std::next(_production_begin, 1)->text;
	_state.Shape(shape_type);
}
void
LightGroup(TranslationState &_state,
		   std::vector<Token>::const_iterator _production_begin,
		   std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "Light group ended, applying semantic action..");
	std::string const	light_type = std::next(_production_begin, 1)->text;
	_state.Light(light_type);
}
void
SamplerGroup(TranslationState &_state,
			 std::vector<Token>::const_iterator _production_begin,
			 std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "Sampler group ended, applying semantic action..");
	std::string const	sampler_type = std::next(_production_begin, 1)->text;
	_state.Sampler(sampler_type);
}
void
IntegratorGroup(TranslationState &_state, 
				std::vector<Token>::const_iterator _production_begin,
				std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "Integrator group ended, applying semantic action..");
	std::string const	integrator_type = std::next(_production_begin, 1)->text;
	_state.Integrator(integrator_type);
}
void
FilmGroup(TranslationState &_state,
		  std::vector<Token>::const_iterator _production_begin,
		  std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "Film group ended, applying semantic action..");
	_state.Film();
}
void
CameraGroup(TranslationState &_state,
			std::vector<Token>::const_iterator _production_begin,
			std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "Camera group ended, applying semantic action..");
	_state.Camera();
}
void
TranslateGroup(TranslationState &_state,
			   std::vector<Token>::const_iterator _production_begin,
			   std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "Translate group ended, applying semantic action..");
	std::vector<Token>::const_iterator const	it = std::next(_production_begin, 1);
	_state.Translate({ maths::stof(it->text),
					   maths::stof(std::next(it, 1)->text),
					   maths::stof(std::next(it, 2)->text) });
}
void
RotateGroup(TranslationState &_state,
			std::vector<Token>::const_iterator _production_begin,
			std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "Rotate group ended, applying semantic action..");
	std::vector<Token>::const_iterator const	it = std::next(_production_begin, 1);
	_state.Rotate(maths::stof(it->text),
				  { maths::stof(std::next(it, 1)->text),
					maths::stof(std::next(it, 2)->text),
					maths::stof(std::next(it, 3)->text) });
}
void
ScaleGroup(TranslationState &_state,
		   std::vector<Token>::const_iterator _production_begin,
		   std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "Scale group ended, applying semantic action..");
	std::vector<Token>::const_iterator const	it = std::next(_production_begin, 1);
	_state.Scale(maths::stof(it->text),
				 maths::stof(std::next(it, 1)->text),
				 maths::stof(std::next(it, 2)->text));
}
void
OutputGroup(TranslationState &_state,
			std::vector<Token>::const_iterator _production_begin,
			std::vector<Token>::const_iterator _production_end)
{
	LOG_ERROR(tools::kChannelParsing, "Output group ended, NO SEMANTIC ACTION");
}

void
IdentityTerminal(TranslationState &_state,
				 std::vector<Token>::const_iterator _production_begin,
				 std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "Identity terminal found, applying semantic action..");
	_state.Identity();
}
void
ScopeBeginTerminal(TranslationState &_state,
				   std::vector<Token>::const_iterator _production_begin,
				   std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "ScopeBegin terminal found, applying semantic action..");
	_state.ScopeBegin();
}
void
ScopeEndTerminal(TranslationState &_state,
				 std::vector<Token>::const_iterator _production_begin,
				 std::vector<Token>::const_iterator _production_end)
{
	LOG_INFO(tools::kChannelParsing, "ScopeEnd terminal found, applying semantic action..");
	_state.ScopeEnd();
}


} // namespace api

