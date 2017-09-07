#pragma once
#ifndef __YS_MEMORY_REGION_HPP__
#define __YS_MEMORY_REGION_HPP__

#include <cstdint>
#include <list>
#include <bitset>

#include "core/core.h"
#include "common_macros.h"
#include "globals.h"
#include "core/profiler.h"
#include "algorithms.h"
#include "core/alloc.h"


namespace core {

class IMemoryRegion
{
public:
	virtual ~IMemoryRegion() = default;
	
	virtual void	ReserveBlocks(size_t _count = 1u) = 0;
	virtual void	*Alloc(size_t _size) = 0;
	//template <typename T> 
	//virtual T		*Alloc(size_t _count = 1u, bool _default_constructed = false) = 0;
};


template <size_t BS, uint64_t AAlign, uint64_t BAlign>
class MemoryRegion final : public IMemoryRegion
{
	static_assert(algo::IsPowerOfTwo(AAlign), "Alignment parameters have to be powers of two.");
	static_assert(algo::IsPowerOfTwo(BAlign), "Alignment parameters have to be powers of two.");
	static_assert(AAlign < BAlign, "Local alignment must be smaller than block alignment.");
public:
	// NOTE: <occupied_space, block_adress>
	using BlockDesc_t = std::pair<size_t, uint8_t *const>;
	using BlockList_t = std::list<BlockDesc_t>;
	static constexpr size_t kBlockSize{ BS };
	static constexpr uint64_t kAllocAlignment{ AAlign };
	static constexpr uint64_t kAllocAlignMask{ kAllocAlignment - 1u };
	static constexpr uint64_t kBlockAlignment{ BAlign };
	static constexpr uint64_t kBlockAlignMask{ kBlockAlignment - 1u };

	~MemoryRegion();

	void	ReserveBlocks(size_t _count = 1u) override;

	void	Clear();

	// TODO: override new operator and replace Alloc everywhere
	void	*Alloc(size_t _size) override;
	template <typename T> T *Alloc(size_t _count, bool _default_constructed);
	template <typename T> T *Alloc(size_t _count = 1u);

private:
	BlockList_t		blocks_;

};


} // namespace core


namespace core {

template <size_t BS, uint64_t AAlign, uint64_t BAlign>
MemoryRegion<BS, AAlign, BAlign>::~MemoryRegion()
{
	for (auto &&pair : blocks_)
		core::FreeAligned(pair.second);
}

template <size_t BS, uint64_t AAlign, uint64_t BAlign>
void
MemoryRegion<BS, AAlign, BAlign>::ReserveBlocks(size_t _count)
{
	TIMED_SCOPE(MemoryRegion_ReserveBlocks);

	for (size_t i = 0; i < _count; ++i)
	{
		uint8_t		*block = core::AllocAligned<uint8_t, kBlockAlignment>(kBlockSize);
		blocks_.emplace_front(0u, block);
	}
}

template <size_t BS, uint64_t AAlign, uint64_t BAlign>
void
MemoryRegion<BS, AAlign, BAlign>::Clear()
{
	for (BlockDesc_t block : blocks_)
		block.first = 0u;
}

template <size_t BS, uint64_t AAlign, uint64_t BAlign>
void *
MemoryRegion<BS, AAlign, BAlign>::Alloc(size_t _size)
{
	YS_ASSERT(_size <= kBlockSize);

	TIMED_SCOPE(MemoryRegion_Alloc);

	for (BlockList_t::iterator it = blocks_.begin(); it != blocks_.end(); ++it)
	{
		BlockDesc_t &pair = *it;

		if (pair.first == 0)
		{
			pair.first = _size;
			return reinterpret_cast<void*>(pair.second);
		}

		if (kBlockSize - pair.first >= _size)
		{
			size_t	base_offset = pair.first;
			size_t	aligned_down = base_offset & ~kAllocAlignMask;
			if (base_offset == aligned_down)
			{
				pair.first += _size;
				return reinterpret_cast<void*>(pair.second + base_offset);
			}
			else
			{
				size_t	aligned_up = aligned_down + kAllocAlignment;
				if (kBlockSize - aligned_up >= _size)
				{
					pair.first = aligned_up + _size;
					return reinterpret_cast<void*>(pair.second + aligned_up);
				}
			}
		}
	}

	ReserveBlocks(1);
	BlockDesc_t &pair = *(blocks_.begin());
	pair.first = _size;
	return reinterpret_cast<void*>(pair.second);
}

template <size_t BS, uint64_t AAlign, uint64_t BAlign>
template <typename T>
T *
MemoryRegion<BS, AAlign, BAlign>::Alloc(size_t _count, bool _default_constructed)
{
	if (_default_constructed)
	{
		T *result = reinterpret_cast<T*>(Alloc(sizeof(T) * _count));
			for (size_t i = 0; i < _count; ++i)
				new (result + i) T();
		return result;
	}
	else
	{
		return Alloc<T>(_count);
	}
}

template <size_t BS, uint64_t AAlign, uint64_t BAlign>
template <typename T>
T *
MemoryRegion<BS, AAlign, BAlign>::Alloc(size_t _count)
{
	return reinterpret_cast<T*>(Alloc(sizeof(T) * _count));
}

} // namespace core


#endif // __YS_MEMORY_REGION_HPP__
