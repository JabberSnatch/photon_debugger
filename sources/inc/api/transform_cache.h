#pragma once
#ifndef __YS_TRANSFORM_CACHE_HPP__
#define __YS_TRANSFORM_CACHE_HPP__

#include <map>

#include "core/nonmovable.h"
#include "core/noncopyable.h"
#include "core/memory_region.h"
#include "maths/transform.h"



namespace api {


class TransformCache final :
	private core::noncopyable,
	private core::nonmovable
{
private:
	using LookupTable_t = std::map<maths::Transform, maths::Transform*>;
public:
	maths::Transform const &Lookup(maths::Transform const &_t);
private:
	LookupTable_t		lookup_table_;
	core::MemoryRegion	mem_region_;
};


} // namespace api


#endif // __YS_TRANSFORM_CACHE_HPP__
