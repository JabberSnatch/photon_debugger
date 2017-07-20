#include "param_set.h"

#include "vector.h"

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
ParamSet::PushInt(std::string const &_id, int32_t *_v, uint32_t _count)
{
	InputParameter<int32_t>	&param = AllocateCommon_(_id, _count, int_parameters_);
	std::copy(_v, _v + _count, param.values);
}
void
ParamSet::PushInt(std::string const &_id, int32_t _v)
{
	PushInt(_id, &_v, 1);
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

maths::Decimal const *
ParamSet::FindFloat(std::string const &_id, uint32_t &_count) const
{
	auto	it = float_parameters_.find(_id);
	if (it != float_parameters_.end())
	{
		_count = it->second.count;
		return it->second.values;
	}

	_count = 0;
	return nullptr;
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

int32_t const *
ParamSet::FindInt(std::string const &_id, uint32_t &_count) const
{
	auto	it = int_parameters_.find(_id);
	if (it != int_parameters_.end())
	{
		_count = it->second.count;
		return it->second.values;
	}

	_count = 0;
	return nullptr;
}
int32_t
ParamSet::FindInt(std::string const &_id, int32_t _default) const
{
	auto	it = int_parameters_.find(_id);
	if (it != int_parameters_.end())
	{
		if (it->second.count == 1)
			return *(it->second.values);
	}

	return _default;
}

bool const *
ParamSet::FindBool(std::string const &_id, uint32_t &_count) const
{
	auto	it = bool_parameters_.find(_id);
	if (it != bool_parameters_.end())
	{
		_count = it->second.count;
		return it->second.values;
	}

	_count = 0;
	return nullptr;
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

void
ParamSet::Clear()
{
	float_parameters_.clear();
	int_parameters_.clear();
	bool_parameters_.clear();
	region_.Clear();
}


} // namespace api
