#include "api/transform_cache.h"

#include "maths/transform.h"


namespace api {


maths::Transform const &
TransformCache::Lookup(maths::Transform const &_t)
{
	LookupTable_t::const_iterator cit = lookup_table_.find(_t);
	if (cit == lookup_table_.cend())
	{
		maths::Transform *const instance = new (mem_region_) maths::Transform{ _t };
		lookup_table_.emplace(*instance, instance);
		return *instance;
	}
	else
		return *(cit->second);
}


} // namespace api
