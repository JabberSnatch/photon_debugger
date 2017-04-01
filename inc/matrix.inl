#ifndef __YS_MATRIX_INL__
#define __YS_MATRIX_INL__

#include "ys_assert.h"


namespace maths
{

template <typename T, uint32_t R, uint32_t C, uint32_t common>
constexpr Matrix<T, R, C>
operator*(Matrix<T, R, common> const &_lhs, Matrix<T, common, C> const &_rhs)
{
	Matrix<T, R, C> result(zero<T>);
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
constexpr Matrix<T, R, C>
Identity()
{
	Matrix<T, R, C> result(zero<T>);
	for (uint32_t i = 0; i < Min(R, C); ++i)
		result.e[i * C + i] = one<T>;
	return result;
}

template <typename T, uint32_t R, uint32_t C>
constexpr Matrix<T, C, R>
Transpose(Matrix<T, R, C> const &_v)
{
	Matrix<T, C, R> result( zero<T> );
	for (uint32_t i = 0; i < R; ++i)
		for (uint32_t j = 0; j < C; ++j)
			result[j][i] = _v[i][j];
	return result;
}

template <typename T, uint32_t RC>
Matrix<T, RC, RC>
Inverse(Matrix<T, RC, RC> const &_v)
{
	// Source is pbrt. Implementation devised in "Numerical Recipes in C", 2.1
	// sample can be found at http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c2-1.pdf
	Matrix<T, RC, RC>	inv{ _v };

	Vector<uint32_t, RC>	pivot(0), indexc(0), indexr(0);

	for (uint32_t i = 0; i < RC; ++i)				// Loop over columns
	{
		uint32_t	icol{ 0 }, irow{ 0 };
		T			big = T{ 0 };
		for (uint32_t j = 0; j < RC; ++j)			// Outer pivot search loop
		{
			if (pivot[j] != 1)
			{
				for (uint32_t k = 0; k < RC; ++k)	// Inner pivot search loop
				{
					if (pivot[k] == 0)
					{
						T	abs_jk = Abs(inv[j][k]);
						if (abs_jk >= big)
						{
							big = abs_jk;
							irow = j; icol = k;
						}
					}
				}
			}
		}
		++pivot[icol];

		if (irow != icol)
			for (uint32_t k = 0; k < RC; ++k)
				std::swap(inv[irow][k], inv[icol][k]);

		indexr[i] = irow;
		indexc[i] = icol;
		
		// NOTE: This is the singular matrix case. Insert log call ASAP.
		YS_ASSERT(inv[icol][icol] != zero<T>);
		if (inv[icol][icol] == zero<T>) return Matrix<T, RC, RC>( zero<T> );

		T	pivot_inv = one<T> / inv[icol][icol];
		inv[icol][icol] = one<T>;
		for (uint32_t j = 0; j < RC; ++j) inv[icol][j] *= pivot_inv;

		for (uint32_t j = 0; j < RC; ++j)
		{
			if (j != icol)
			{
				T	backup = inv[j][icol];
				inv[j][icol] = zero<T>;
				for (uint32_t k = 0; k < RC; ++k) inv[j][k] -= inv[icol][k] * backup;
			}
		}
	}

	for (uint32_t i = 1; i <= RC; ++i)
	{
		uint32_t j = RC - i;
		if (indexr[j] != indexc[j])
			for (uint32_t k = 0; k < RC; ++k)
				std::swap(inv[k][indexr[j]], inv[k][indexc[j]]);
	}

	return inv;
}

} // namespace matrix
} // namespace maths


#endif // __YS_MATRIX_INL__
