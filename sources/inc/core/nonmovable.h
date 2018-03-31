#pragma once
#ifndef __YS_NONMOVABLE_HPP__
#define __YS_NONMOVABLE_HPP__


namespace core {


class nonmovable
{
protected:
	nonmovable() = default;
	~nonmovable() = default;
	nonmovable(nonmovable&&) = delete;
	nonmovable &operator = (nonmovable&&) = delete;
};


} // namespace core 


#endif // __YS_NONMOVABLE_HPP__
