#ifndef __YS_MATHS_HPP__
#define __YS_MATHS_HPP__


#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

namespace maths
{

template <typename T> static constexpr T LowestValue = std::numeric_limits<T>::lowest();
template <typename T> static constexpr T HighestValue = std::numeric_limits<T>::max();
template <typename T> static constexpr T Infinity = std::numeric_limits<T>::infinity();

template <typename T> static constexpr T Pi = T( 3.14159235658979323846 );

template <typename T>
constexpr T Radians(T _degrees) { return (Pi<T> / T( 180 )) * _degrees; }
template <typename T>
constexpr T Degrees(T _radians) { return (T( 180 ) / Pi<T>) * _radians; }


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


template <typename T>
constexpr T AlmostOne() { return one<T> - std::numeric_limits<T>::epsilon(); }

// NOTE: Maybe a Scalar<T> class could make these functions a little more specific.
template <typename T>
constexpr T Min(T _lhs, T _rhs) { return (_lhs < _rhs) ? _lhs : _rhs; }
template <typename T>
constexpr T Max(T _lhs, T _rhs) { return (_lhs > _rhs) ? _lhs : _rhs; }
template <typename T>
constexpr T Clamp(T _v, T _min, T _max) { return Min(Max(_v, _min), _max); }

template <typename T>
constexpr T Abs(T _v) { return (_v > zero<T>) ? _v : -_v; }


template <typename T>
constexpr T Lerp(T _a, T _b, float _t) { return _a*(1.f - _t) + _b*_t; }
template <typename T>
constexpr T Lerp(T _a, T _b, double _t) { return _a*(1. - _t) + _b*_t; }

} // namespace maths


#endif // __YS_MATHS_HPP__
