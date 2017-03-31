#ifndef __YS_MATRIX_INL__
#define __YS_MATRIX_INL__


namespace maths
{

template <typename T, uint32_t R, uint32_t C, uint32_t common>
constexpr Matrix<T, R, C>
operator*(Matrix<T, R, common> const &_lhs, Matrix<T, common, C> const &_rhs)
{
	Matrix<T, R, C> result(T{ 0 });
	for (uint32_t i = 0; i < R; ++i)
		for (uint32_t j = 0; j < C; ++j)
			for (uint32_t k = 0; k < common; ++k)
				result[i][j] += _lhs[i][k] * _rhs[k][j];
	return result;
}
template <typename T, uint32_t R, uint32_t C>
constexpr Vector<T, R>
operator*(Matrix<T, R, C> const &_lhs, Vector<T, C> const &_rhs)
{
	Vector<T, C> result{};
	for (uint32_t i = 0; i < R; ++i)
		for (uint32_t j = 0; j < C; ++j)
			result[i] += _lhs[i][j] * _rhs[j];
	return result;
}

namespace matrix
{

template <typename T, uint32_t R, uint32_t C>
constexpr Matrix<T, C, R>
Transpose(Matrix<T, R, C> const &_v)
{
	Matrix<T, C, R> result(T{ 0 });
	for (uint32_t i = 0; i < R; ++i)
		for (uint32_t j = 0; j < C; ++j)
			result[j][i] = _v[i][j];
	return result;
}
template <typename T, uint32_t RC>
constexpr Matrix<T, RC, RC>
Inverse(Matrix<T, RC, RC> const &_v)
{

}

} // namespace matrix
} // namespace maths


#endif // __YS_MATRIX_INL__
