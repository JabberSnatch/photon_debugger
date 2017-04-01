#ifndef __YS_MATHS_HPP__
#define __YS_MATHS_HPP__


namespace maths
{

template <typename T> struct Zero { static constexpr T value = T::zero ; };
template <> struct Zero<float> { static constexpr float value = 0.f; };
template <> struct Zero<double> { static constexpr double value = 0.; };
template <> struct Zero<int> { static constexpr int value = 0; };
template <typename T> constexpr T zero = Zero<T>::value;

template <typename T> struct One { static constexpr T value = T::one; };
template <> struct One<float> { static constexpr float value = 1.f; };
template <> struct One<double> { static constexpr double value = 1.; };
template <> struct One<int> { static constexpr int value = 1; };
template <typename T> constexpr T one = One<T>::value;

/*
template <typename T> constexpr T zero = T::zero;
template <> constexpr float zero<float> = 0.f;
template <> constexpr double zero<double> = 0.0;
template <> constexpr int32_t zero<int32_t> = 0;
template <> constexpr int64_t zero<int64_t> = 0l;
template <> constexpr uint32_t zero<uint32_t> = 0u;
template <> constexpr uint64_t zero<uint64_t> = 0ul;
*/

template <typename T>
constexpr T AlmostOne() { return one<T> - std::numeric_limits<T>::epsilon(); }

template <typename T>
constexpr T Min(T _lhs, T _rhs) { return (_lhs < _rhs) ? _lhs : _rhs; }
template <typename T>
constexpr T Max(T _lhs, T _rhs) { return (_lhs > _rhs) ? _lhs : _rhs; }
template <typename T>
constexpr T Clamp(T _v, T _min, T _max) { return Min(Max(_v, _min), _max); }

template <typename T>
constexpr T Abs(T _v) { return (_v > zero<T>) ? _v : -_v; }

} // namespace maths


#endif // __YS_MATHS_HPP__
