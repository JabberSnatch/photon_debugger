#pragma once
#ifndef __YS_SPINLOCK_HPP__
#define __YS_SPINLOCK_HPP__

#include <atomic>

namespace tools {


class AtomicSpinLock final
{
public:
	inline void Acquire() { uint32_t value = 1u; while (value == 1u) value = atomic_.exchange(1u); }
	inline void Wait() const { uint32_t value = 1u; while (value == 1u) value = atomic_.load(); }
	inline void Release() { atomic_.exchange(0u); }
	inline bool is_locked()	const { return atomic_.load() > 0u; }
private:
	std::atomic_uint32_t	atomic_{ 0u };
};


} // namespace tools

#endif // __YS_SPINLOCK_HPP__
