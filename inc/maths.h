#ifndef __YS_MATHS_HPP__
#define __YS_MATHS_HPP__


namespace maths
{


template <typename T>
constexpr T Min(T _lhs, T _rhs) { return (_lhs < _rhs) ? _lhs : _rhs; }
template <typename T>
constexpr T Max(T _lhs, T _rhs) { return (_lhs > _rhs) ? _lhs : _rhs; }
template <typename T>
constexpr T Clamp(T _v, T _min, T _max) { return Min(Max(_v, _min), _max); }

template <typename T>
constexpr T AlmostOne() { return T{ 1 } - std::numeric_limits<T>::epsilon(); }


} // namespace maths


#endif // __YS_MATHS_HPP__
