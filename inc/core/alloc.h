#pragma once
#ifndef __YS_ALLOC_HPP__
#define __YS_ALLOC_HPP__

#include <cstdint>


namespace core {

void *AllocAligned(size_t const _size, size_t const _alignment = 16);
void FreeAligned(void *const _address);
template <typename T>
T *AllocAligned(size_t const _count = 1, size_t const _alignment = 16);


} // namespace core


namespace core {


template <typename T>
T *
AllocAligned(size_t const _count, size_t const _alignment)
{
	T *const alloc = reinterpret_cast<T*>(AllocAligned(sizeof(T) * _count, _alignment));
	return alloc;
}


} // namespace core


#endif // __YS_ALLOC_HPP__
