#ifndef __YS_ALGORITHMS_HPP__
#define __YS_ALGORITHMS_HPP__

#include <array>


namespace algo
{


template <size_t N>
struct fill
{
	template <typename T, typename ...Tn>
	static constexpr auto apply(T head, Tn ...tail)
	{
		return fill<N - 1>::apply(head, head, tail...);
	}
};

template <>
struct fill<1>
{
	template <typename T, typename ...Tn>
	static constexpr auto apply(T head, Tn ...tail)
	{
		return std::array<T, sizeof...(tail)+1>{head, tail...};
	}
};

template <typename T>
constexpr bool IsPowerOfTwo(T _value)
{
	return _value && ((_value & (_value - 1)) == 0);
}


} // namespace algo


#endif // __YS_ALGORITHMS_HPP__
