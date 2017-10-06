#include "api/param_set.h"

#include "maths/vector.h"

namespace api {

void
ParamSet::PushFloat(std::string const &_id, maths::Decimal *_v, uint32_t _count)
{
	InputParameter<maths::Decimal> &param = AllocateCommon_(_id, _count, float_parameters_);
	std::copy(_v, _v + _count, param.values);
}
void
ParamSet::PushFloat(std::string const &_id, maths::Decimal _v)
{
	PushFloat(_id, &_v, 1);
}

void
ParamSet::PushInt(std::string const &_id, int64_t *_v, uint32_t _count)
{
	InputParameter<int64_t>	&param = AllocateCommon_(_id, _count, int_parameters_);
	std::copy(_v, _v + _count, param.values);
}
void
ParamSet::PushInt(std::string const &_id, int64_t _v)
{
	PushInt(_id, &_v, 1);
}

void
ParamSet::PushUint(std::string const &_id, uint64_t *_v, uint32_t _count)
{
	InputParameter<uint64_t>	&param = AllocateCommon_(_id, _count, uint_parameters_);
	std::copy(_v, _v + _count, param.values);
}
void
ParamSet::PushUint(std::string const &_id, uint64_t _v)
{
	PushUint(_id, &_v, 1);
}

void
ParamSet::PushBool(std::string const &_id, bool *_v, uint32_t _count)
{
	InputParameter<bool>	&param = AllocateCommon_(_id, _count, bool_parameters_);
	std::copy(_v, _v + _count, param.values);
}
void
ParamSet::PushBool(std::string const &_id, bool _v)
{
	PushBool(_id, &_v, 1);
}

void
ParamSet::PushString(std::string const &_id, std::string const &_v)
{
	using IteratorBoolPair_t = std::pair<StringMap_t::const_iterator, bool>;
	IteratorBoolPair_t result = string_parameters_.emplace(_id, _v);
	if (!result.second)
	{
		LOG_ERROR(tools::kChannelGeneral, "Tried to push a string on a pre-existing identifier");
		YS_ASSERT(false);
	}
}


void
ParamSet::PushTransform(std::string const &_id, maths::Transform const &_transform)
{
	using IteratorBoolPair_t = std::pair<TransformMap_t::const_iterator, bool>;
	IteratorBoolPair_t result = transforms_.emplace(_id, &_transform);
	if (!result.second)
	{
		LOG_ERROR(tools::kChannelGeneral, "Tried to push a Transform on a pre-existing identifier");
		YS_ASSERT(false);
	}
}

maths::Decimal
ParamSet::FindFloat(std::string const &_id, maths::Decimal _default) const
{
	auto	it = float_parameters_.find(_id);
	if (it != float_parameters_.end())
	{
		if (it->second.count == 1)
			return *(it->second.values);
	}

	return _default;
}

int64_t
ParamSet::FindInt(std::string const &_id, int64_t _default) const
{
	auto	it = int_parameters_.find(_id);
	if (it != int_parameters_.end())
	{
		if (it->second.count == 1)
			return *(it->second.values);
	}

	return _default;
}

uint64_t
ParamSet::FindUint(std::string const &_id, uint64_t _default) const
{
	auto	it = uint_parameters_.find(_id);
	if (it != uint_parameters_.end())
	{
		if (it->second.count == 1)
			return *(it->second.values);
	}

	return _default;
}

bool
ParamSet::FindBool(std::string const &_id, bool _default) const
{
	auto	it = bool_parameters_.find(_id);
	if (it != bool_parameters_.end())
	{
		if (it->second.count == 1)
			return *(it->second.values);
	}

	return _default;
}

std::string const &
ParamSet::FindString(std::string const &_id, std::string const &_default) const
{
	StringMap_t::const_iterator scit = string_parameters_.find(_id);
	if (scit != string_parameters_.cend())
	{
		return scit->second;
	}
	else
	{
		return _default;
	}
}


maths::Transform const &
ParamSet::FindTransform(std::string const &_id, maths::Transform const &_default) const
{
	TransformMap_t::const_iterator tcit = transforms_.find(_id);
	if (tcit != transforms_.cend())
	{
		return *tcit->second;
	}
	else
	{
		return _default;
	}
}

void
ParamSet::Clear()
{
	float_parameters_.clear();
	int_parameters_.clear();
	bool_parameters_.clear();
	string_parameters_.clear();
	transforms_.clear();
	region_.Clear();
}


} // namespace api
