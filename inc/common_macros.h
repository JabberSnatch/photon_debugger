#pragma once
#ifndef __YS_COMMON_MACROS_HPP__
#define __YS_COMMON_MACROS_HPP__

//#define YS_DECIMAL_IS_DOUBLE
//#define YS_NO_LOGS
//#define QUADRATIC_RUNNING_ERROR_ON_DISCRIMINANT


#ifdef _WIN32
#define YS_WINDOWS
#endif // _WIN32


#ifdef _DEBUG
#define YS_DEBUG

#ifndef YS_DECIMAL_IS_DOUBLE
#define YS_REDECIMAL_HAS_PRECISE
#endif // YS_DECIMAL_IS_DOUBLE

#endif // _DEBUG


#include <cassert>
#ifdef YS_WINDOWS
#ifdef YS_DEBUG
#define YS_ASSERT(statement) assert(statement)
#else
#define YS_ASSERT(statement) {}
#endif
#endif // _WIN32


#endif // __YS_COMMON_MACROS_HPP__
