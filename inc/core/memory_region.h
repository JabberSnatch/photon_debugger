#pragma once
#ifndef __YS_MEMORY_REGION_HPP__
#define __YS_MEMORY_REGION_HPP__

#include <cstdint>
#include <list>


namespace core {


class MemoryRegion final
{
public:
	// NOTE: <occupied_space, block_adress>
	using BlockDesc_t = std::pair<size_t, uint8_t *const>;
	using BlockList_t = std::list<BlockDesc_t>;
public:
	MemoryRegion();
	MemoryRegion(size_t const _block_size);
	MemoryRegion(size_t const _block_size, uint64_t const _alloc_align, uint64_t const _block_align);
	~MemoryRegion();
	MemoryRegion(MemoryRegion&&) = default;
	MemoryRegion& operator = (MemoryRegion&&) = default;
	MemoryRegion(const MemoryRegion&) = delete;
	MemoryRegion& operator = (const MemoryRegion&) = delete;
public:
	void ReserveBlocks(size_t const _count = 1u);
	void Clear();
	void *Alloc(size_t const _size);
private:
	BlockList_t		blocks_;
	const size_t	block_size_;
	const uint64_t	alloc_align_;
	const uint64_t	alloc_mask_;
	const uint64_t	block_align_;
	const uint64_t	block_mask_;
};


} // namespace core

void* operator new (size_t _count, core::MemoryRegion &_region);
void* operator new[](size_t _count, core::MemoryRegion &_region);
void operator delete (void* _ptr, core::MemoryRegion &_region);
void operator delete[](void* _ptr, core::MemoryRegion &_region);


#endif // __YS_MEMORY_REGION_HPP__
