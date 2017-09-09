#include "core/memory_region.h"

#include "algorithms.h"
#include "common_macros.h"
#include "globals.h"
#include "core/profiler.h"
#include "core/alloc.h"


namespace core {

static constexpr size_t kDefaultBlockSize = 256 * 1024;
static constexpr uint64_t kDefaultAllocAlign = 16u;
static constexpr uint64_t kDefaultBlockAlign = 64u;


MemoryRegion::MemoryRegion() :
	MemoryRegion(kDefaultBlockSize, kDefaultAllocAlign, kDefaultBlockAlign)
{}


MemoryRegion::MemoryRegion(size_t const _block_size) :
	MemoryRegion(_block_size, kDefaultAllocAlign, kDefaultBlockAlign)
{}


MemoryRegion::MemoryRegion(size_t const _block_size,
						   uint64_t const _alloc_align,
						   uint64_t const _block_align) :
	blocks_{},
	block_size_{ _block_size },
	alloc_align_{ _alloc_align },
	alloc_mask_{ alloc_align_ - 1u },
	block_align_{ _block_align },
	block_mask_{ block_align_ - 1u }
{
	YS_ASSERT(algo::IsPowerOfTwo(alloc_align_));
	YS_ASSERT(algo::IsPowerOfTwo(block_align_));
	YS_ASSERT(alloc_align_ < block_align_);
}

MemoryRegion::~MemoryRegion()
{
	for (auto &&pair : blocks_)
		core::FreeAligned(pair.second);
}

void
MemoryRegion::ReserveBlocks(size_t const _count)
{
	TIMED_SCOPE(MemoryRegion_ReserveBlocks);
	for (size_t block_index = 0; block_index < _count; ++block_index)
	{
		uint8_t *const block = core::AllocAligned<uint8_t>(block_size_, block_align_);
		blocks_.emplace_front(0u, block);
	}
}

void
MemoryRegion::Clear()
{
	for (BlockDesc_t &block : blocks_)
		block.first = 0u;
}

void*
MemoryRegion::Alloc(size_t const _size)
{
	YS_ASSERT(_size <= block_size_);
	TIMED_SCOPE(MemoryRegion_Alloc);
	void *result = nullptr;
	for (BlockList_t::iterator bit = blocks_.begin(); bit != blocks_.end(); ++bit)
	{
		BlockDesc_t &desc = *bit;
		if (desc.first == 0)
		{
			desc.first = _size;
			result = static_cast<void*>(desc.second);
		}
		else if (block_size_ - desc.first >= _size)
		{
			size_t const base_offset = desc.first;
			size_t const aligned_down = base_offset & ~alloc_mask_;
			if (base_offset == aligned_down)
			{
				desc.first += _size;
				result = static_cast<void*>(desc.second + base_offset);
			}
			else
			{
				size_t const aligned_up = aligned_down + alloc_align_;
				if (block_size_ - aligned_up >= _size)
				{
					desc.first = aligned_up + _size;
					return static_cast<void*>(desc.second + aligned_up);
				}
			}
		}
	}
	if (result == nullptr)
	{
		ReserveBlocks(1u);
		BlockDesc_t &desc = blocks_.front();
		desc.first = _size;
		result = static_cast<void*>(desc.second);
	}
	return result;
}

} // namespace core


void*
operator new (size_t _count, core::MemoryRegion &_region)
{
	return _region.Alloc(_count);
}
void*
operator new[] (size_t _count, core::MemoryRegion &_region)
{
	return _region.Alloc(_count);
}
void
operator delete (void* _ptr, core::MemoryRegion &_region)
{}
void
operator delete[](void* _ptr, core::MemoryRegion &_region)
{}