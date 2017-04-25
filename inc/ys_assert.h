#ifndef __YS_ASSERT_HPP__
#define __YS_ASSERT_HPP__

#include <cassert>


#ifdef _WIN32

#ifdef _DEBUG
#define YS_ASSERT(statement) assert(statement)
#else
#define YS_ASSERT(statement) {}
#endif

#endif // _WIN32


#endif // __YS_ASSERT_HPP__
