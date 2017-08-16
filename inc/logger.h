#pragma once
#ifndef __YS_LOGGER_HPP__
#define __YS_LOGGER_HPP__

#include <cstdint>
#include <string>
#include <array>
#include <chrono>
#include <map>
#include <set>
#include <vector>
#include <memory>
#include <algorithm>
#include <fstream>

#include "algorithms.h"
#include "common_macros.h"
#include "spinlock.h"
#include "globals.h"
#include "profiler.h"

namespace tools {

enum LogChannel { kChannelGeneral = 0, kChannelProfiling, kChannelParsing, kChannelCount, kChannelStdOut };
enum LogLevel { kLevelDebug = 0, kLevelInfo, kLevelWarning, kLevelError, kLevelCount };

// One instance of FileLogger manages an array of stringstreams.
// When a buffer exceeds a certain size, it is written back to disk and cleared.
// LocalLoggers give it messages, and it orders them before writing them back to disk.
// But really how much ordering do we want / can we actually get ?
// A thread that issues a lot of messages might cause an early flush, and we don't really
// want to sort messages once they are written to an output file.
// Each thread could maintain its own set of files in a binary format (?) and they would be merged
// after execution, when no log can be issued anymore.
// Another idea could be that a central flush would cause every thread local logger to flush into
// the central logger before writing to disk, much like a memory fence.
// Problem is if threads are busy doing something else, they cannot necessarily handle the request.
// Either we find a way to halt execution for logging, or loggers are allocated on shared memory,
// and each thread is given an index into the local table.

struct LogEntry final
{
	using Clock_t = std::chrono::high_resolution_clock;
	using TimeStamp_t = Clock_t::duration;

	LogEntry() : channel{ kChannelCount }, level{ kLevelCount }, message{}, date{} {}
	LogEntry(LogChannel _channel, LogLevel _level, std::string const &_message) :
		channel{ _channel }, level{ _level }, message{ _message },
		date{ Clock_t::now() - start_time }
	{}

	inline bool operator<(LogEntry const &_other) const { return date < _other.date; }


	LogChannel	channel;
	LogLevel	level;
	std::string	message;
	TimeStamp_t	date;

private:
	static Clock_t::time_point	start_time;
};

namespace logging {
using ChannelLocks_t = std::array<core::AtomicSpinLock, LogChannel::kChannelCount>;
} // namespace logging

template <size_t hs, size_t bs>
class Logger;

template <size_t hs, size_t bs>
class LogBuffer final
{
	friend class Logger<hs, bs>;
public:
	static constexpr size_t		kBufferSize = bs;

	LogBuffer(Logger<hs, bs> &_writer);

	LogBuffer(LogBuffer<hs, bs> const &) = delete;
	LogBuffer &operator=(LogBuffer<hs, bs> const &) = delete;

private:
	void	Log(LogChannel _channel, LogLevel _level, std::string const &_message);

	using EntryBuffer_t = std::array<LogEntry, bs>;
	using ChannelBuffers_t = std::array<EntryBuffer_t, LogChannel::kChannelCount>;
	using ChannelSizes_t = std::array<size_t, LogChannel::kChannelCount>;

	Logger<hs, bs>			&log_writer_;
	ChannelBuffers_t		entry_buffer_;
	ChannelSizes_t			entry_count_ = algo::fill<kChannelCount>::apply<size_t>(0);

	logging::ChannelLocks_t		channel_lock_;
};


template <size_t hs, size_t bs>
class Logger final
{
	using EntryBuffer_t = typename LogBuffer<hs, bs>::EntryBuffer_t;
public:
	static constexpr size_t		kHistorySize = hs;
	static constexpr size_t		kBufferSize = bs;
	static constexpr size_t		kInvalidThreadIndex = std::numeric_limits<size_t>::max();
	thread_local static size_t	thread_index;

	// Allocates necessary resources for parallel logging. This is called once by a main thread.
	// Remember to set thread_index on each client thread.
	void	AllowMultipleThreads(size_t _thread_count);

	void	Log(LogChannel _channel, LogLevel _level, char const *_message);
	void	Log(LogChannel _channel, LogLevel _level, std::string const &_message);
	void	Merge(LogChannel _channel, EntryBuffer_t &_buffer, size_t &_entry_count);
	void	Flush(LogChannel _channel, bool _release_merge_lock = false);
	void	FlushAll();

	void	BindPath(LogChannel _channel, std::string const &_path, bool _clear_file = true);
	void	UnbindPath(LogChannel _channel, std::string const &_path);

	void	WaitFlushEnd(LogChannel _channel) const;
	bool	IsFlushing(LogChannel _channel) const;

	void		EnableChannel(LogChannel _channel, bool _enabled);
#ifdef YS_DEBUG
	LogLevel	lowest_active_level{ LogLevel::kLevelDebug };
#else
	LogLevel	lowest_active_level{ LogLevel::kLevelInfo };
#endif

private:
	using ChannelOutputs_t = std::multimap<LogChannel, std::string>;
	using LogHistory_t = std::multiset<LogEntry>;
	using ChannelHistories_t = std::array<LogHistory_t, LogChannel::kChannelCount>;
	using pLogBuffer = std::unique_ptr<LogBuffer<hs, bs>>;

	std::multimap<LogChannel, std::string>		paths_;
	std::array<bool, LogChannel::kChannelCount>	is_enabled_{ algo::fill<kChannelCount>::apply(true) };

	size_t										thread_count_{ 0 };
	std::vector<pLogBuffer>						thread_buffer_;

	ChannelHistories_t							log_history_;

	logging::ChannelLocks_t						merge_lock_;
	logging::ChannelLocks_t						flush_lock_;
	logging::ChannelLocks_t						write_lock_;

	logging::ChannelLocks_t						debug_history_lock_;
};

} // namespace tools


// IMPLEMENTATION
namespace tools {

template <size_t hs, size_t bs> thread_local size_t Logger<hs, bs>::thread_index = Logger<hs, bs>::kInvalidThreadIndex;

inline std::ostream &operator<<(std::ostream &_stream, LogEntry const &_entry)
{
	TIMED_SCOPE(LogEntry_StreamOperator);

	std::string level_string;
	switch (_entry.level) {
	case kLevelDebug:
		level_string = "  DEBUG";
		break;
	case kLevelInfo:
		level_string = "   INFO";
		break;
	case kLevelError:
		level_string = "  ERROR";
		break;
	case kLevelWarning:
		level_string = "WARNING";
		break;
	}
	_stream << _entry.date.count() << " [" << level_string << "] : " << _entry.message;
	return _stream;
}


template <size_t hs, size_t bs>
LogBuffer<hs, bs>::LogBuffer(Logger<hs, bs> &_writer) :
	log_writer_{ _writer }
{}


template <size_t hs, size_t bs>
void
LogBuffer<hs, bs>::Log(LogChannel _channel, LogLevel _level, std::string const &_message)
{
	log_writer_.WaitFlushEnd(_channel);

	channel_lock_[_channel].Acquire();
	entry_buffer_[_channel][entry_count_[_channel]++] = LogEntry{ _channel, _level, _message };
	//size_t entry_count = ++entry_count_[_channel];
	channel_lock_[_channel].Release();

	if (entry_count_[_channel] >= kBufferSize && !log_writer_.IsFlushing(_channel))
	{
		//entry_count_[_channel] = 0;
		//EntryBuffer_t	entry_buffer = std::move(entry_buffer_[_channel]);
		//EntryBuffer_t	entry_buffer = entry_buffer_[_channel];
		log_writer_.Merge(_channel, entry_buffer_[_channel], entry_count_[_channel]);
	}
	//else
	//	channel_lock_[_channel].Release();
}



template <size_t hs, size_t bs>
void
Logger<hs, bs>::AllowMultipleThreads(size_t _thread_count)
{
	thread_count_ = _thread_count;
	thread_buffer_.clear();
	thread_buffer_.reserve(thread_count_);
	for (size_t i = 0; i < thread_count_; ++i)
		thread_buffer_.emplace_back(new LogBuffer<hs, bs>(*this));
}
template <size_t hs, size_t bs>
void
Logger<hs, bs>::Log(LogChannel _channel, LogLevel _level, char const *_message)
{
	Log(_channel, _level, std::string{ _message });
}
template <size_t hs, size_t bs>
void
Logger<hs, bs>::Log(LogChannel _channel, LogLevel _level, std::string const &_message)
{
	TIMED_SCOPE(Logger_Log);

	if (_level < lowest_active_level) return;
	if (!is_enabled_[_channel]) return;
	YS_ASSERT(paths_.find(_channel) != paths_.end());

	if (thread_index != kInvalidThreadIndex)
		YS_ASSERT(thread_count_ > 0);

	if (thread_count_ > 0 && thread_index != kInvalidThreadIndex)
	{
		YS_ASSERT(thread_index < thread_count_);
		thread_buffer_[thread_index]->Log(_channel, _level, _message);
	}
	else
	{
		log_history_[_channel].emplace(LogEntry{ _channel, _level, _message });
		if (log_history_[_channel].size() >= kHistorySize)
			Flush(_channel);
	}
}
template <size_t hs, size_t bs>
void
Logger<hs, bs>::Merge(LogChannel _channel, EntryBuffer_t &_buffer, size_t &_entry_count)
{
	TIMED_SCOPE(Logger_Merge);

	merge_lock_[_channel].Acquire();
	
	//debug_history_lock_[_channel].Acquire();
	//log_history_[_channel].emplace(LogEntry{ _channel, kLevelDebug, "Thread " + std::to_string(thread_index) + " started merging" });
	//debug_history_lock_[_channel].Release();

	for (size_t i = 0; i < _entry_count; ++i)
		log_history_[_channel].emplace(std::move(_buffer[i]));
	_entry_count = 0;

	if (log_history_[_channel].size() >= kHistorySize && !IsFlushing(_channel))
		Flush(_channel, true);
	else
		merge_lock_[_channel].Release();

	//debug_history_lock_[_channel].Acquire();
	//log_history_[_channel].emplace(LogEntry{ _channel, kLevelDebug, "Thread " + std::to_string(thread_index) + " finished merging" });
	//debug_history_lock_[_channel].Release();
}
template <size_t hs, size_t bs>
void
Logger<hs, bs>::Flush(LogChannel _channel, bool _release_merge_lock)
{
	TIMED_SCOPE(Logger_Flush);

	flush_lock_[_channel].Acquire();

	//debug_history_lock_[_channel].Acquire();
	//log_history_[_channel].emplace(LogEntry{ _channel, kLevelDebug, "Thread " + std::to_string(thread_index) + " started flushing" });
	//debug_history_lock_[_channel].Release();

	// NOTE: This is the only place where a thread will access another thread's entries
	for (size_t i = 0; i < thread_count_; ++i)
	{
		thread_buffer_[i]->channel_lock_[_channel].Acquire();
		size_t entry_count = thread_buffer_[i]->entry_count_[_channel];
		for (size_t j = 0; j < entry_count; ++j)
			log_history_[_channel].emplace(std::move(thread_buffer_[i]->entry_buffer_[_channel][j]));
		thread_buffer_[i]->entry_count_[_channel] = 0;
		thread_buffer_[i]->channel_lock_[_channel].Release();
	}

	LogHistory_t	history = std::move(log_history_[_channel]);
	//LogHistory_t	history_copy = log_history_[_channel];
	//log_history_[_channel].clear();
	
	flush_lock_[_channel].Release();
	write_lock_[_channel].Acquire();
	if (_release_merge_lock) merge_lock_[_channel].Release();

	for (auto it = paths_.find(_channel); it != paths_.end(); ++it)
	{
		if (it->first != _channel)
			break;
		std::ofstream	file_stream{ it->second, std::ios_base::app | std::ios_base::out };
		for (auto &&entry : history)
		{
			file_stream << entry << std::endl;
		}
		file_stream.close();
	}

	//debug_history_lock_[_channel].Acquire();
	//log_history_[_channel].emplace(LogEntry{ _channel, kLevelDebug, "Thread " + std::to_string(thread_index) + " finished flushing" });
	//debug_history_lock_[_channel].Release();

	write_lock_[_channel].Release();
}
template <size_t hs, size_t bs>
void
Logger<hs, bs>::FlushAll()
{
	for (int i = 0; i < kChannelCount; ++i)
		Flush(LogChannel(i));
}

template <size_t hs, size_t bs>
void
Logger<hs, bs>::BindPath(LogChannel _channel, std::string const &_path, bool _clear_file)
{
	for (auto it = paths_.find(_channel); it != paths_.end(); ++it)
	{
		if (it->first != _channel)
			break;
		if (it->second == _path)
			return;
	}

	paths_.emplace(_channel, _path);
	if (_clear_file)
		std::ofstream file_stream{ _path };
}
template <size_t hs, size_t bs>
void
Logger<hs, bs>::UnbindPath(LogChannel _channel, std::string const &_path)
{
	for (auto it = paths_.find(_channel); it != paths_.end(); ++it)
	{
		if (it->first != _channel)
			break;
		if (it->second == _path)
		{
			paths_.erase(it);
			return;
		}
	}
}

template <size_t hs, size_t bs>
void
Logger<hs, bs>::WaitFlushEnd(LogChannel _channel) const
{
	flush_lock_[_channel].Wait();
}
template <size_t hs, size_t bs>
bool
Logger<hs, bs>::IsFlushing(LogChannel _channel) const
{
	return flush_lock_[_channel].is_locked();
}

template <size_t hs, size_t bs>
void
Logger<hs, bs>::EnableChannel(LogChannel _channel, bool _enabled)
{
	is_enabled_[_channel] = _enabled;
}


} // namespace tools

#endif // __YS_LOGGER_HPP__
