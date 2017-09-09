#pragma once
#ifndef __YS_PARAM_SET_HPP__
#define __YS_PARAM_SET_HPP__

#include <cstdint>
#include <unordered_map>
#include <string>

#include "maths/maths.h"
#include "maths/vector.h"
#include "core/logger.h"
#include "globals.h"
#include "core/memory_region.h"


namespace api {

template <typename T>
struct InputParameter
{
	T				*const values;
	uint32_t const	count;
};


class ParamSet final
{
public:
	ParamSet() = default;
	~ParamSet() = default;
	ParamSet(ParamSet const &) = delete;
	ParamSet(ParamSet &&) = delete;
	ParamSet &operator=(ParamSet const &) = delete;
	ParamSet &operator=(ParamSet &&) = delete;

	void	PushFloat(std::string const &_id, maths::Decimal *_v, uint32_t _count);
	void	PushFloat(std::string const &_id, maths::Decimal _v);
	template <uint32_t size>
	void	PushFloat(std::string const &_id, maths::Vector<maths::Decimal, size> _v);

	void	PushInt(std::string const &_id, int32_t *_v, uint32_t _count);
	void	PushInt(std::string const &_id, int32_t _v);
	template <uint32_t size>
	void	PushInt(std::string const &_id, maths::Vector<int32_t, size> _v);

	void	PushBool(std::string const &_id, bool *_v, uint32_t _count);
	void	PushBool(std::string const &_id, bool _v);
	template <uint32_t size>
	void	PushBool(std::string const &_id, maths::Vector<bool, size> _v);

	maths::Decimal const *FindFloat(std::string const &_id, uint32_t &_count) const;
	maths::Decimal FindFloat(std::string const &_id, maths::Decimal _default) const;
	template <uint32_t size>
	maths::VectorF<size> FindFloat(std::string const &_id,
								   maths::VectorF<size> const &_default) const;

	int32_t const *FindInt(std::string const &_id, uint32_t &_count) const;
	int32_t FindInt(std::string const &_id, int32_t _default) const;
	template <uint32_t size>
	maths::VectorI<size> FindInt(std::string const &_id,
								 maths::VectorI<size> const &_default) const;

	bool const *FindBool(std::string const &_id, uint32_t &_count) const;
	bool FindBool(std::string const &_id, bool _default) const;
	template <uint32_t size>
	maths::VectorB<size> FindBool(std::string const &_id,
								  maths::VectorB<size> const &_default) const;

	void	Clear();

private:
	std::unordered_map<std::string, InputParameter<maths::Decimal>>	float_parameters_;
	std::unordered_map<std::string, InputParameter<int32_t>>		int_parameters_;
	std::unordered_map<std::string, InputParameter<bool>>			bool_parameters_;

	core::MemoryRegion<>	region_;


	template <typename T>
	inline InputParameter<T> &AllocateCommon_(
		std::string const &_id, uint32_t _count,
		std::unordered_map<std::string, InputParameter<T>> &_map);
};


template <uint32_t size>
void
ParamSet::PushFloat(std::string const &_id, maths::Vector<maths::Decimal, size> _v)
{
	PushFloat(_id, _v.e.data(), _v.size);
}

template <uint32_t size>
void
ParamSet::PushInt(std::string const &_id, maths::Vector<int32_t, size> _v)
{
	PushInt(_id, _v.e.data(), _v.size);
}

template <uint32_t size>
void
ParamSet::PushBool(std::string const &_id, maths::Vector<bool, size> _v)
{
	PushBool(_id, _v.e.data(), _v.size);
}

template <uint32_t size>
maths::VectorF<size>
ParamSet::FindFloat(std::string const &_id, maths::VectorF<size> const &_default) const
{
	auto	it = float_parameters_.find(_id);
	if (it != float_parameters_.end())
	{
		if (it->second.count == size)
		{
			maths::VectorF<size>	result{};
			std::copy(it->second.values, it->second.values + it->second.count, result.e.begin());
			return result;
		}
	}

	return _default;
}

template <uint32_t size>
maths::VectorI<size>
ParamSet::FindInt(std::string const &_id,
				  maths::VectorI<size> const &_default) const
{
	auto	it = int_parameters_.find(_id);
	if (it != int_parameters_.end())
	{
		if (it->second.count == size)
		{
			maths::VectorI<size>	result{};
			std::copy(it->second.values, it->second.values + it->second.count, result.e.begin());
			return result;
		}
	}

	return _default;
}

template <uint32_t size>
maths::VectorB<size>
ParamSet::FindBool(std::string const &_id,
				   maths::VectorB<size> const &_default) const
{
	auto	it = bool_parameters_.find(_id);
	if (it != bool_parameters_.end())
	{
		if (it->second.count == size)
		{
			maths::VectorB<size>	result{};
			std::copy(it->second.values, it->second.values + it->second.count, result.e.begin());
			return result;
		}
	}

	return _default;
}

template <typename T>
InputParameter<T> &
ParamSet::AllocateCommon_(std::string const &_id, uint32_t _count,
						  std::unordered_map<std::string, InputParameter<T>> &_map)
{
	YS_ASSERT(_count != 0);

	auto	it = _map.find(_id);
	if (it != _map.end())
	{
		LOG_WARNING(tools::kChannelGeneral, "Allocating on an existing field, value will be overwritten.");
		if (_count == it->second.count)
			return it->second;
		else
			_map.erase(it);
	}

	T *location = new (region_) T();
	InputParameter<T>	parameter{ location, _count };

	return _map.emplace(_id, parameter).first->second;
}


} // namespace api


#endif // __YS_PARAM_SET_HPP__
