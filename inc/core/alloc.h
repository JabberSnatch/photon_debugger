#pragma once
#ifndef __YS_ALLOC_HPP__
#define __YS_ALLOC_HPP__

#include <cstdint>
#include "common_macros.h"


namespace core {

template <size_t Alignment = 16>
void			*AllocAligned(size_t _size);
template <typename T, size_t Alignment = 16>
T				*AllocAligned(size_t _count = 1, bool _default_constructed = false);
inline void		FreeAligned(void *_address);


} // namespace core


namespace core {

template <size_t Alignment>
void *
AllocAligned(size_t _size)
{
	static_assert(algo::IsPowerOfTwo(Alignment), "Alignment must be a power of two.");
	constexpr size_t kAlignMask = Alignment - 1;

	uint8_t		*alloc_begin = new uint8_t[_size + Alignment + sizeof(size_t)];
#ifdef YS_DEBUG
	memset(alloc_begin, 0xcau, _size + Alignment + sizeof(size_t));
#endif
	
	size_t		alloc_begin_address = reinterpret_cast<size_t>(alloc_begin);
	size_t		aligned_begin_address =
		(alloc_begin_address + kAlignMask + sizeof(size_t)) & ~kAlignMask;
	YS_ASSERT(alloc_begin_address <= aligned_begin_address);

	uint8_t		*alloc = reinterpret_cast<uint8_t*>(aligned_begin_address);
	*(reinterpret_cast<size_t*>(alloc) - 1) = alloc_begin_address;

	return reinterpret_cast<void*>(alloc);
}

template <typename T, size_t Alignment>
T *
AllocAligned(size_t _count, bool _default_constructed)
{
	T *alloc = reinterpret_cast<T*>(AllocAligned<Alignment>(sizeof(T) * _count));
	if (_default_constructed)
		for (size_t i = 0u; i < _count; ++i)
			new (alloc + i) T();
	return alloc;
}

inline void
FreeAligned(void *_address)
{
	uint8_t		*alloc = reinterpret_cast<uint8_t*>(_address);
	size_t		alloc_begin_address = *(reinterpret_cast<size_t*>(alloc) - 1);
	uint8_t		*alloc_begin = reinterpret_cast<uint8_t*>(alloc_begin_address);
	delete[] alloc_begin;
}

} // namespace core


#endif // __YS_ALLOC_HPP__
