#pragma once
#ifndef __YS_MEMORY_REGION_HPP__
#define __YS_MEMORY_REGION_HPP__

#include <cstdint>
#include <list>
#include <bitset>

#include "common_macros.h"
#include "globals.h"
#include "profiler.h"


namespace core {


template <size_t BS = 256 * 1024, uint64_t AAlign = 16, uint64_t BAlign = 64>
class MemoryRegion
{
public:
	using BlockDesc_t = std::pair<size_t, uint8_t*>;
	using BlockList_t = std::list<BlockDesc_t>;
	static constexpr size_t kBlockSize{ BS };
	static constexpr uint64_t kAllocAlignment{ AAlign };
	static constexpr uint64_t kAllocAlignMask{ kAllocAlignment - 1u };
	static constexpr uint64_t kBlockAlignment{ BAlign };
	static constexpr uint64_t kBlockAlignMask{ kBlockAlignment - 1u };


#ifdef YS_DEBUG
	MemoryRegion()
	{
		std::bitset<sizeof(kAllocAlignment)*8>	alloc_alignment_bitset(kAllocAlignment);
		YS_ASSERT(alloc_alignment_bitset.count() == 1);
		std::bitset<sizeof(kBlockAlignment)*8>	block_alignment_bitset(kBlockAlignment);
		YS_ASSERT(block_alignment_bitset.count() == 1);
		YS_ASSERT(kAllocAlignment < kBlockAlignment);
	}
#endif
	~MemoryRegion();

	void	ReserveBlocks(size_t _count = 1u);

	void	*Alloc(size_t _size);
	template <typename T> T *Alloc(size_t _count = 1u, bool _default_constructed = false);

private:
	BlockList_t		blocks_;

};


} // namespace core


namespace core {

template <size_t BS, uint64_t AAlign, uint64_t BAlign>
MemoryRegion<BS, AAlign, BAlign>::~MemoryRegion()
{
	for (auto &&pair : blocks_)
	{
		size_t		block_begin_address = *(reinterpret_cast<size_t*>(pair.second) - 1);
		uint8_t		*block_begin = reinterpret_cast<uint8_t*>(block_begin_address);
		delete[] block_begin;
	}
}

template <size_t BS, uint64_t AAlign, uint64_t BAlign>
void
MemoryRegion<BS, AAlign, BAlign>::ReserveBlocks(size_t _count)
{
	TIMED_SCOPE(MemoryRegion_ReserveBlocks);

	for (size_t i = 0; i < _count; ++i)
	{
		uint8_t		*block_begin = new uint8_t[kBlockSize + kBlockAlignment + sizeof(size_t)];
#ifdef YS_DEBUG
		memset(block_begin, 0xcau, kBlockSize + kBlockAlignment + sizeof(size_t));
#endif

		size_t		block_begin_address = reinterpret_cast<size_t>(block_begin);
		size_t		aligned_begin_address =
			(block_begin_address + kBlockAlignMask + sizeof(size_t)) & ~kBlockAlignMask;
		YS_ASSERT(block_begin_address <= aligned_begin_address);

		uint8_t		*block = reinterpret_cast<uint8_t*>(aligned_begin_address);
		*(reinterpret_cast<size_t*>(block) - 1) = block_begin_address;

		blocks_.emplace_front(0u, block);
	}
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
	T *result = reinterpret_cast<T*>(Alloc(sizeof(T) * _count));
	if (_default_constructed)
		for (size_t i = 0; i < _count; ++i)
			new (result + i) T();
	return result;
}

} // namespace core


#endif // __YS_MEMORY_REGION_HPP__
