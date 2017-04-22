#ifndef __YS_QUATERNION_INL__
#define __YS_QUATERNION_INL__


namespace maths
{

template <typename T>
Quaternion<T>::Quaternion(Transform const &_v)
{
	// Variant of Shoemake, 1991, as implemented in pbrt (Pharr)
	Mat4x4f const &m = _v.m();
	float trace = m[0][0] + m[1][1] + m[2][2];
	if (trace > 0.f)
	{
		// Compute w from matrix trace, then xyz
		float s = std::sqrt(trace + 1.0f);
		w = s * 0.5f;
		s = 0.5f / s;
		v.x = (m[2][1] - m[1][2]) * s;
		v.y = (m[0][2] - m[2][0]) * s;
		v.z = (m[1][0] - m[0][1]) * s;
	}
	else
	{
		int const nxt[3] = { 1, 2, 0 };
		float q[3];
		int i = 0;
		if (m[1][1] > m[0][0]) i = 1;
		if (m[2][2] > m[i][i]) i = 2;
		int j = nxt[i];
		int k = nxt[j];
		float s = std::sqrt((m[i][i] - (m[j][j] + m[k][k])) + 1.f);
		q[i] = s * 0.5f;
		if (s != 0.f) s = 0.5f / s;
		w = (m[k][j] - m[j][k]) * s;
		q[j] = (m[j][i] + m[i][j]) * s;
		q[k] = (m[k][i] + m[i][k]) * s;
		v.x = q[0];
		v.y = q[1];
		v.z = q[2];
	}
}
template <typename T>
constexpr
Quaternion<T>::operator Transform() const
{
	T	xx{ v.x * v.x }, yy{ v.y * v.y }, zz{ v.z * v.z },
		xy{ v.x * v.y }, xz{ v.x * v.z }, xw{ v.x * w },
		yz{ v.y * v.z }, yw{ v.y * w }, zw{ v.z * w };
	Mat4x4f result{};

	result[0][0] = one<T> - 2 * (yy + zz);
	result[0][1] = 2 * (xy + zw);
	result[0][2] = 2 * (xz - yw);
	result[1][0] = 2 * (xy - zw);
	result[1][1] = one<T> - 2 * (xx + zz);
	result[1][2] = 2 * (yz + xw);
	result[2][0] = 2 * (xz + yw);
	result[2][1] = 2 * (yz - xw);
	result[2][2] = one<T> - 2 * (xx + yy);

	return Transform{ matrix::Transpose(result), result };
}


template <typename T>
constexpr bool
operator==(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs)
{
	return _lhs.v == _rhs.v && _lhs.w == _rhs.w;
}
template <typename T>
constexpr bool
operator!=(Quaternion<T> const &_lhs, Quaternion<T> const &_rhs)
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
constexpr Quaternion<T>
Conjugate(Quaternion<T> const &_v)
{
	return Quaternion<T>{-_v.v, _v.w};
}
template <typename T>
constexpr Quaternion<T>
Inverse(Quaternion<T> const &_v)
{
	return Conjugate(_v) / Norm(_v);
}
template <typename T>
constexpr T
Norm(Quaternion<T> const &_v)
{
	return Dot(_v, _v);
}
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
