#ifndef __YS_BOUNDS_INL__
#define __YS_BOUNDS_INL__


namespace maths
{

template <typename T, uint32_t n>
Point<T, n> const&
Bounds<T, n>::operator[](uint32_t _i) const
{
	YS_ASSERT(_i < 2);
	if (_i == 0) return min;
	else return max;
}
template <typename T, uint32_t n>
Point<T, n>&
Bounds<T, n>::operator[] (uint32_t _i)
{
	YS_ASSERT(_i < 2);
	if (_i == 0) return min;
	else return max;
}

template <typename T, uint32_t n>
Point<T, n>
Bounds<T, n>::Corner(uint32_t _index) const
{
	YS_ASSERT(_index < (1 << n)); // Bounds have 2^n corners.
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = (*this)[(_index & (1 << i) ? 1 : 0)][i];
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
Bounds<T, n>::Diagonal() const
{
	return max - min;
}
template <typename T, uint32_t n>
T
Bounds<T, n>::SurfaceArea() const
{
	Vector<T, n> d{ Diagonal() };
	T result{ zero<T> };
	for (uint32_t i = 0; i < n; ++i)
		for (uint32_t j = i + 1; j < n; ++j)
			result += d[i] * d[j];
	return result * T{ 2 };
}
template <typename T, uint32_t n>
T
Bounds<T, n>::Volume() const
{
	return FoldProduct(Diagonal());
}
template <typename T, uint32_t n>
uint32_t
Bounds<T, n>::MaximumExtent() const
{
	return MaximumDimension(Diagonal());
}

template <typename T, uint32_t n>
Point<T, n>
Bounds<T, n>::Lerp(Point<float, n> const &_t) const
{
	Point<T, n> result{};
	for (uint32_t i = 0; i < n; ++i)
		result[i] = maths::Lerp(min[i], max[i], _t[i]);
	return result;
}
template <typename T, uint32_t n>
Vector<T, n>
Bounds<T, n>::Offset(Point<T, n> const &_p) const
{
	Vector<T, n> result = _p - min;
	Vector<T, n> d{ Diagonal() };
	for (uint32_t i = 0; i < n; ++i)
		if (d[i] > zero<T>)
			result[i] /= d[i];
	YS_ASSERT(!result.HasNaNs());
	return result;
}

template <typename T, uint32_t n>
void
Bounds<T, n>::BoundingSphere(Point<T, n> &_center, T &_radius) const
{
	_center = (min + max) / T{ 2 };
	_radius = Inside(*this, _center) ? Distance(_center, max) : zero<T>;
}


template <typename T, uint32_t n>
bool
operator==(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs)
{
	return _lhs.min == _rhs.min && _lhs.max == _rhs.max;
}
template <typename T, uint32_t n>
bool
operator!=(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs)
{
	return _lhs.min != _rhs.min || _lhs.max != _rhs.max;
}


template <typename T, uint32_t n>
Bounds<T, n>
Union(Bounds<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	Point<T, n> min, max;
	for (uint32_t i = 0; i < n; ++i)
	{
		min[i] = maths::Min(_lhs.min[i], _rhs[i]);
		max[i] = maths::Max(_lhs.max[i], _rhs[i]);
	}
	return Bounds<T, n>{min, max};
}
template <typename T, uint32_t n>
Bounds<T, n>
Union(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs)
{
	Point<T, n> min, max;
	for (uint32_t i = 0; i < n; ++i)
	{
		min[i] = maths::Min(_lhs.min[i], _rhs.min[i]);
		max[i] = maths::Max(_lhs.max[i], _rhs.max[i]);
	}
	return Bounds<T, n>{min, max};
}
template <typename T, uint32_t n>
Bounds<T, n>
Intersect(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs)
{
	Point<T, n> min, max;
	for (uint32_t i = 0; i < n; ++i)
	{
		min[i] = maths::Max(_lhs.min[i], _rhs.min[i]);
		max[i] = maths::Min(_lhs.max[i], _rhs.max[i]);
	}
	return Bounds<T, n>{min, max};
}
template <typename T, uint32_t n>
bool
Overlap(Bounds<T, n> const &_lhs, Bounds<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		if ((_lhs.max[i] < _rhs.min[i]) || (_lhs.min[i] > _rhs.max[i]))
			return false;
	return true;
}
template <typename T, uint32_t n>
bool
Inside(Bounds<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		if ((_lhs.max[i] < _rhs[i]) || (_lhs.min[i] > _rhs[i]))
			return false;
	return true;
}
template <typename T, uint32_t n>
bool
InsideExclusive(Bounds<T, n> const &_lhs, Point<T, n> const &_rhs)
{
	for (uint32_t i = 0; i < n; ++i)
		if ((_lhs.max[i] <= _rhs[i]) || (_lhs.min[i] > _rhs[i]))
			return false;
	return true;
}
template <typename T, uint32_t n>
Bounds<T, n>
Expand(Bounds<T, n> const &_v, T _delta)
{
	return Bounds<T, n>{_v.min - Vector<T, n>(_delta), _v.max + Vector<T, n>(_delta)};
}

} // namespace maths


#endif // __YS_BOUNDS_INL__
