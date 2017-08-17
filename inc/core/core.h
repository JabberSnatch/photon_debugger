#pragma once
#ifndef __YS_CORE_HPP__
#define __YS_CORE_HPP__

#include <cstdint>


namespace core {

class AtomicSpinLock;

class IMemoryRegion;
template <size_t BS = 256 * 1024, uint64_t AAlign = 16, uint64_t BAlign = 64>
class MemoryRegion;

} // namespace core

#endif // __YS_CORE_HPP__
