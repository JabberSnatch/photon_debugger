#include "quaternion.h"

namespace maths
{

Quaternion::Quaternion(Transform const &_v)
{
	// Variant of Shoemake, 1991, as implemented in pbrt (Pharr)
	Mat4x4f const &m = _v.m();
	Decimal trace = m[0][0] + m[1][1] + m[2][2];
	if (trace > 0._d)
	{
		// Compute w from matrix trace, then xyz
		Decimal s = std::sqrt(trace + 1._d);
		w = s * 0.5_d;
		s = 0.5_d / s;
		v.x = (m[2][1] - m[1][2]) * s;
		v.y = (m[0][2] - m[2][0]) * s;
		v.z = (m[1][0] - m[0][1]) * s;
	}
	else
	{
		int const nxt[3] = { 1, 2, 0 };
		Decimal q[3];
		int i = 0;
		if (m[1][1] > m[0][0]) i = 1;
		if (m[2][2] > m[i][i]) i = 2;
		int j = nxt[i];
		int k = nxt[j];
		Decimal s = std::sqrt((m[i][i] - (m[j][j] + m[k][k])) + 1._d);
		q[i] = s * 0.5_d;
		if (s != 0._d) s = 0.5_d / s;
		w = (m[k][j] - m[j][k]) * s;
		q[j] = (m[j][i] + m[i][j]) * s;
		q[k] = (m[k][i] + m[i][k]) * s;
		v.x = q[0];
		v.y = q[1];
		v.z = q[2];
	}
}
Quaternion::operator Transform() const
{
	Decimal	xx{ v.x * v.x }, yy{ v.y * v.y }, zz{ v.z * v.z },
		xy{ v.x * v.y }, xz{ v.x * v.z }, xw{ v.x * w },
		yz{ v.y * v.z }, yw{ v.y * w }, zw{ v.z * w };
	Mat4x4f result{};

	result[0][0] = 1._d - 2 * (yy + zz);
	result[0][1] = 2 * (xy + zw);
	result[0][2] = 2 * (xz - yw);
	result[1][0] = 2 * (xy - zw);
	result[1][1] = 1._d - 2 * (xx + zz);
	result[1][2] = 2 * (yz + xw);
	result[2][0] = 2 * (xz + yw);
	result[2][1] = 2 * (yz - xw);
	result[2][2] = 1._d - 2 * (xx + yy);

	return Transform{ matrix::Transpose(result), result };
}


bool
operator==(Quaternion const &_lhs, Quaternion const &_rhs)
{
	return _lhs.v == _rhs.v && _lhs.w == _rhs.w;
}
bool
operator!=(Quaternion const &_lhs, Quaternion const &_rhs)
{
	return _lhs.v != _rhs.v || _lhs.w != _rhs.w;
}

Quaternion&
operator+=(Quaternion &_lhs, Quaternion const &_rhs)
{
	_lhs.v += _rhs.v;
	_lhs.w += _rhs.w;
	return _lhs;
}
Quaternion&
operator-=(Quaternion &_lhs, Quaternion const &_rhs)
{
	_lhs.v -= _rhs.v;
	_lhs.w -= _rhs.w;
	return _lhs;
}
Quaternion&
operator*=(Quaternion &_lhs, Quaternion const &_rhs)
{
	_lhs = _lhs * _rhs;
	return _lhs;
}

Quaternion
operator+(Quaternion const &_lhs, Quaternion const &_rhs)
{
	return Quaternion{ _lhs.v + _rhs.v, _lhs.w + _rhs.w };
}
Quaternion
operator-(Quaternion const &_lhs, Quaternion const &_rhs)
{
	return Quaternion{ _lhs.v - _rhs.v, _lhs.w - _rhs.w };
}
Quaternion
operator*(Quaternion const &_lhs, Quaternion const &_rhs)
{
	return Quaternion{
		vector::Cross(_lhs.v, _rhs.v) + _lhs.v * _rhs.w + _rhs.v * _lhs.w,
		_lhs.w * _rhs.w + vector::Dot(_lhs.v, _rhs.v)
	};
}

Quaternion&
operator*=(Quaternion &_lhs, Decimal _rhs)
{
	_lhs.v *= _rhs;
	_lhs.w *= _rhs;
	return _lhs;
}
Quaternion&
operator/=(Quaternion &_lhs, Decimal _rhs)
{
	_lhs.v /= _rhs;
	_lhs.w /= _rhs;
	return _lhs;
}

Quaternion
operator*(Quaternion const &_lhs, Decimal _rhs)
{
	return Quaternion{ _lhs.v * _rhs, _lhs.w * _rhs };
}
Quaternion
operator*(Decimal _lhs, Quaternion const &_rhs)
{
	return _rhs * _lhs;
}
Quaternion
operator/(Quaternion const &_lhs, Decimal _rhs)
{
	return Quaternion{ _lhs.v / _rhs, _lhs.w / _rhs };
}


namespace quaternion
{

Quaternion
Conjugate(Quaternion const &_v)
{
	return Quaternion{ -_v.v, _v.w };
}
Quaternion
Inverse(Quaternion const &_v)
{
	return Conjugate(_v) / Norm(_v);
}
Decimal
Norm(Quaternion const &_v)
{
	return Dot(_v, _v);
}
Decimal
Dot(Quaternion const &_lhs, Quaternion const &_rhs)
{
	return vector::Dot(_lhs.v, _rhs.v) + _lhs.w * _rhs.w;
}
Quaternion
Normalized(Quaternion const &_v)
{
	return _v / std::sqrt(Dot(_v, _v));
}
Quaternion
Slerp(Quaternion const &_a, Quaternion const &_b, Decimal _t)
{
	// Source is PBR, Pharr
	Decimal	cos_theta = Dot(_a, _b);
	if (cos_theta > .9995_d)
		return Normalized(Lerp(_a, _b, _t));
	else
	{
		Decimal	theta = std::acos(Clamp(cos_theta, -1._d, 1._d));
		Decimal	thetap = theta * _t;
		Quaternion qperp = Normalized(_b - _a * cos_theta);
		return _a * std::cos(thetap) + qperp * std::sin(thetap);
	}
}

} // namespace quaternion
} // namespace maths

