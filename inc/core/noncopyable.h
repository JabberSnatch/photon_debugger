#pragma once
#ifndef __YS_NONCOPYABLE_HPP__
#define __YS_NONCOPYABLE_HPP__


namespace core {


class noncopyable
{
protected:
	noncopyable() = default;
	~noncopyable() = default;
	noncopyable(noncopyable const &) = delete;
	noncopyable &operator = (noncopyable const &) = delete;
};


} // namespace core


#endif // __YS_NONCOPYABLE_HPP__
