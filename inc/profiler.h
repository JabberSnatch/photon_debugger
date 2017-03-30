#ifndef __YS_PROFILER_HPP__
#define __YS_PROFILER_HPP__

#ifdef _WIN32
#include "win32_timer.h"
#endif // _WIN32

#include <cstdint>
#include <unordered_map>

namespace tools
{

class Profiler
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

	Timer&
	GetTimer(char const *_key)
	{
		uint64_t	hash{ FNV_Hash{}(_key) };
		auto		ite{ timers_.find(hash) };
		if (ite != timers_.end())
			return ite->second;
		return timers_.emplace(hash, Timer{ _key }).first->second;
	}

	TimerTable_t const &
	timers() const { return timers_; }

private:
	TimerTable_t	timers_{};

};

} // namespace tools

thread_local tools::Profiler	gProfiler{};

#endif // __YS_PROFILER_HPP__
