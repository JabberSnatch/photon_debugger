#pragma once
#ifndef __YS_INPUT_PARSER_HPP__
#define __YS_INPUT_PARSER_HPP__

#include <string>



namespace api {

class TranslationState;

bool ProcessInputFile(std::string const &_path, TranslationState &_state);

} // namespace api


#endif // __YS_INPUT_PARSER_HPP__
