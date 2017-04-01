#ifndef __YS_QUATERNION_INL__
#define __YS_QUATERNION_INL__


namespace maths
{

template <typename T>
constexpr Quaternion<T>
&operator+=(Quaternion<T> &_lhs, Quaternion<T> const &_rhs)
{
	_lhs.v += _rhs.v;
	_lhs.w += _rhs.w;
	return *lhs;
}
template <typename T>
constexpr Quaternion<T>
&operator-=(Quaternion<T> &_lhs, Quaternion<T> const &_rhs)
{
	_lhs.v -= _rhs.v;
	_lhs.w -= _rhs.w;
	return *lhs;
}
template <typename T>
constexpr Quaternion<T>
&operator*=(Quaternion<T> &_lhs, Quaternion<T> const &_rhs)
{
	_lhs = _lhs * _rhs;
	return *lhs;
}

template <typename T>
constexpr Quaternion<T>
operator+(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs)
{
	return Quaternion<T>{ _lhs.v + _rhs.v, _lhs.w + _rhs.w };
}
template <typename T>
constexpr Quaternion<T>
operator-(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs)
{
	return Quaternion<T>{ _lhs.v - _rhs.v, _lhs.w - _rhs.w };
}
template <typename T>
constexpr Quaternion<T>
operator*(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs)
{
	return Quaternion<T>{
		vector::Cross(_lhs.v, _rhs.v) + _lhs.v * _rhs.w + _rhs.v * _lhs.w,
		_lhs.w * _rhs.w + vector::Dot(_lhs.v, _rhs.v)
	};
}

template <typename T>
constexpr Quaternion<T>
&operator*=(Quaternion<T> &_lhs, T _rhs)
{
	_lhs.v *= _rhs;
	_lhs.w *= _rhs;
	return *_lhs;
}
template <typename T>
constexpr Quaternion<T>
&operator/=(Quaternion<T> &_lhs, T _rhs)
{
	_lhs.v /= _rhs;
	_lhs.w /= _rhs;
	return *_lhs;
}

template <typename T>
constexpr Quaternion<T>
operator*(Quaternion<T> const &_lhs, T _rhs)
{
	return Quaternion<T>{ _lhs.v * _rhs, _lhs.w * _rhs };
}
template <typename T>
constexpr Quaternion<T>
operator*(T _lhs, Quaternion<T> const &_rhs)
{
	return _rhs * _lhs;
}
template <typename T>
constexpr Quaternion<T>
operator/(Quaternion<T> const &_lhs, T _rhs)
{
	return Quaternion<T>{ _lhs.v / _rhs, _lhs.w / _rhs };
}


namespace quaternion
{

template <typename T>
constexpr T
Dot(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs)
{
	return vector::Dot(_lhs.v, _rhs.v) + _lhs.w * _rhs.w;
}
template <typename T>
constexpr Quaternion<T>
Normalized(Quaternion<T> const &_v)
{
	return _v / std::sqrt(Dot(_v, _v));
}

} // namespace quaternion
} // namespace maths


#endif // __YS_QUATERNION_INL__
