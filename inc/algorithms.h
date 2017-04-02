#ifndef __YS_ALGORITHMS_HPP__
#define __YS_ALGORITHMS_HPP__


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

} // namespace algo


#endif // __YS_ALGORITHMS_HPP__
