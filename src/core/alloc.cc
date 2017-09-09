#include "core/alloc.h"

#include "common_macros.h"


namespace core {


void *
AllocAligned(size_t const _size, size_t const _alignment)
{
	YS_ASSERT(algo::IsPowerOfTwo(_alignement));
	size_t const align_mask = _alignment - 1;
	uint8_t *const alloc_begin = new uint8_t[_size + _alignment + sizeof(size_t)];
#ifdef YS_DEBUG
	memset(alloc_begin, 0xcau, _size + Alignment + sizeof(size_t));
#endif
	size_t const alloc_begin_address = reinterpret_cast<size_t>(alloc_begin);
	size_t const aligned_begin_address =
		(alloc_begin_address + align_mask + sizeof(size_t)) & ~align_mask;
	YS_ASSERT(alloc_begin_address <= aligned_begin_address);
	uint8_t *const alloc = reinterpret_cast<uint8_t*>(aligned_begin_address);
	*(reinterpret_cast<size_t*>(alloc) - 1) = alloc_begin_address;
	return reinterpret_cast<void*>(alloc);
}


void
FreeAligned(void *const _address)
{
	uint8_t *const alloc = reinterpret_cast<uint8_t*>(_address);
	size_t const alloc_begin_address = *(reinterpret_cast<size_t*>(alloc) - 1);
	uint8_t *const alloc_begin = reinterpret_cast<uint8_t*>(alloc_begin_address);
	delete[] alloc_begin;
}


} // namespace core
