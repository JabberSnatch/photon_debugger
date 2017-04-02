#ifndef __YS_QUATERNION_INL__
#define __YS_QUATERNION_INL__


namespace maths
{

template <typename T>
constexpr bool operator==(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs)
{
	return _lhs.v == _rhs.v && _lhs.w == _rhs.w;
}
template <typename T>
constexpr bool operator!=(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs)
{
	return _lhs.v != _rhs.v || _lhs.w != _rhs.w;
}

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
template <typename T>
Quaternion<T>
Slerp(Quaternion<T> const &_a, Quaternion<T> const &_b, float _t)
{
	// Source is PBR, Pharr
	float	cos_theta = Dot(_a, _b);
	if (cos_theta > .9995f)
		return Normalized(Lerp(_a, _b, _t));
	else
	{
		float	theta = std::acos(Clamp(cos_theta, -1.f, 1.f));
		float	thetap = theta * _t;
		Quaternion qperp = Normalized(_b - _a * cos_theta);
		return _a * std::cos(thetap) + qperp * std::sin(thetap);
	}
}

} // namespace quaternion
} // namespace maths


#endif // __YS_QUATERNION_INL__
