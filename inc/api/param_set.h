#pragma once
#ifndef __YS_PARAM_SET_HPP__
#define __YS_PARAM_SET_HPP__

#include <cstdint>
#include <unordered_map>
#include <string>

#include "core/logger.h"
#include "core/memory_region.h"
#include "core/noncopyable.h"
#include "globals.h"
#include "maths/maths.h"
#include "maths/vector.h"


namespace api {

template <typename T>
struct InputParameter
{
	T				*const values;
	uint32_t const	count;
};


class ParamSet final :
	core::noncopyable
{
public:
	void	PushFloat(std::string const &_id, maths::Decimal *_v, uint32_t _count);
	void	PushFloat(std::string const &_id, maths::Decimal _v);
	template <uint32_t size>
	void	PushFloat(std::string const &_id, maths::VectorF<size> _v);

	void	PushInt(std::string const &_id, int64_t *_v, uint32_t _count);
	void	PushInt(std::string const &_id, int64_t _v);
	template <uint32_t size>
	void	PushInt(std::string const &_id, maths::VectorI<size> _v);

	void	PushUint(std::string const &_id, uint64_t *_v, uint32_t _count);
	void	PushUint(std::string const &_id, uint64_t _v);
	template <uint32_t size>
	void	PushUint(std::string const &_id, maths::VectorU<size> _v);

	void	PushBool(std::string const &_id, bool *_v, uint32_t _count);
	void	PushBool(std::string const &_id, bool _v);
	template <uint32_t size>
	void	PushBool(std::string const &_id, maths::VectorB<size> _v);

	void	PushString(std::string const &_id, std::string const &_v);

	void	PushTransform(std::string const &_id, maths::Transform const &_transform);

	maths::Decimal FindFloat(std::string const &_id, maths::Decimal _default) const;
	template <uint32_t size>
	maths::VectorF<size> FindFloat(std::string const &_id,
								   maths::VectorF<size> const &_default) const;

	int64_t FindInt(std::string const &_id, int64_t _default) const;
	template <uint32_t size>
	maths::VectorI<size> FindInt(std::string const &_id,
								 maths::VectorI<size> const &_default) const;

	uint64_t FindUint(std::string const &_id, uint64_t _default) const;
	template <uint32_t size>
	maths::VectorU<size> FindUint(std::string const &_id,
								  maths::VectorU<size> const &_default) const;

	bool FindBool(std::string const &_id, bool _default) const;
	template <uint32_t size>
	maths::VectorB<size> FindBool(std::string const &_id,
								  maths::VectorB<size> const &_default) const;

	std::string const &FindString(std::string const &_id, std::string const &_default) const;

	maths::Transform const &FindTransform(std::string const &_id,
										  maths::Transform const &_default) const;

	void	Clear();

private:
	template <typename T>
	using ParamMap_t = std::unordered_map<std::string, InputParameter<T>>;
	using FloatMap_t = ParamMap_t<maths::Decimal>;
	using IntMap_t = ParamMap_t<int64_t>;
	using UintMap_t = ParamMap_t<uint64_t>;
	using BoolMap_t = ParamMap_t<bool>;
	using StringMap_t = std::unordered_map<std::string, std::string>;
	using TransformMap_t = std::unordered_map<std::string, maths::Transform const *const>;
	FloatMap_t		float_parameters_;
	IntMap_t		int_parameters_;
	UintMap_t		uint_parameters_;
	BoolMap_t		bool_parameters_;
	StringMap_t		string_parameters_;
	TransformMap_t	transforms_;

	core::MemoryRegion region_;

	template <typename T>
	inline InputParameter<T> &AllocateCommon_(
		std::string const &_id, uint32_t _count,
		std::unordered_map<std::string, InputParameter<T>> &_map);
};


template <uint32_t size>
void
ParamSet::PushFloat(std::string const &_id, maths::VectorF<size> _v)
{
	PushFloat(_id, _v.e.data(), _v.size);
}

template <uint32_t size>
void
ParamSet::PushInt(std::string const &_id, maths::VectorI<size> _v)
{
	PushInt(_id, _v.e.data(), _v.size);
}

template <uint32_t size>
void
ParamSet::PushUint(std::string const &_id, maths::VectorU<size> _v)
{
	PushUint(_id, _v.e.data(), _v.size);
}

template <uint32_t size>
void
ParamSet::PushBool(std::string const &_id, maths::VectorB<size> _v)
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
maths::VectorU<size>
ParamSet::FindUint(std::string const &_id,
				   maths::VectorU<size> const &_default) const
{
	auto	it = uint_parameters_.find(_id);
	if (it != uint_parameters_.end())
	{
		if (it->second.count == size)
		{
			maths::VectorU<size>	result{};
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
						  ParamMap_t<T> &_map)
{
	YS_ASSERT(_count != 0);

	ParamMap_t<T>::iterator it = _map.find(_id);
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
