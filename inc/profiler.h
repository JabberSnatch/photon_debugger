#ifndef __YS_PROFILER_HPP__
#define __YS_PROFILER_HPP__

#include "common_macros.h"

#ifdef YS_WINDOWS
#include "win32_timer.h"
#endif // YS_WINDOWS

#include <cstdint>
#include <unordered_map>

#include "spinlock.h"

namespace tools
{

class Profiler final
{
public:
	template<typename T>
	struct IdentityHash
	{
		constexpr T operator() (T const &_v) const { return _v; }
	};

	struct FNV_Hash
	{
		static constexpr uint64_t kPrime{ 1099511628211Ui64 };
		static constexpr uint64_t kOffset{ 0xcbf29ce484222325Ui64 };
		constexpr uint64_t operator() (char const *_v) const
		{
			uint64_t	hash{ kOffset };
			char const	*ite{ _v };
			while (*ite != '\0')
			{
				hash *= kPrime;
				hash ^= *ite;
				++ite;
			}
			return hash;
		}
	};

	using TimerTable_t = std::unordered_map<uint64_t, Timer, IdentityHash<uint64_t>>;

	// Merges _other into caller.
	// NOTE: The whole design is losely protected against race conditions. Keep in mind that
	//		 you should always start a merge when you are definitely sure that _other won't
	//		 change during the merge.
	void GrabTimers(Profiler &_other)
	{
		lock_.Acquire();

		for (auto &&pair : _other.timers_)
		{
			auto it = timers_.find(pair.first);
			if (it != timers_.end())
				it->second.Add(pair.second);
			else
				timers_.emplace(pair);
		}

		lock_.Release();
	}

	Timer &GetTimer(char const *_key)
	{
		uint64_t	hash{ FNV_Hash{}(_key) };
		auto		ite{ timers_.find(hash) };
		if (ite != timers_.end())
			return ite->second;
		return timers_.emplace(hash, Timer{ _key }).first->second;
	}

	void Clear()
	{
		timers_.clear();
	}

	TimerTable_t const &timers() const { return timers_; }

private:
	TimerTable_t	timers_{};
	AtomicSpinLock	lock_;

};

} // namespace tools


#endif // __YS_PROFILER_HPP__
