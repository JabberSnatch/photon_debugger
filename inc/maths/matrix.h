#ifndef __YS_MATRIX_HPP__
#define __YS_MATRIX_HPP__

#include "algorithms.h"
#include "maths/vector.h"
#include "maths/maths.h"


namespace maths
{


// | 00,		01,		   .., 0(rows-1)		|
// | 10,		11,		   .., 1(rows-1)		|
// | ..,		..,		   ..,        ..		|
// | (cols-1)0, (cols-1)0, .., (cols-1)(rows-1) |
//					A
//					|
//					V
// | 0,				1,					.,	rows-1			|
// | rows,			rows+1,				.,	rows+rows-1		|
// | .,				.,					.,	.				|
// | (cols-1)*rows,	(cols-1)*rows+1,	.,	(cols)*rows-1	|

// m[i][j] <=> t[i * cols + j]

// Matrices can be cast from size to size. Any new rows and columns are initialized to
// identity.

template <typename T, uint32_t R, uint32_t C>
struct Matrix final
{
	Matrix();
	Matrix(std::initializer_list<T> _args);

	explicit constexpr Matrix(T _v);

	template <uint32_t sR, uint32_t sC>
	explicit Matrix(Matrix<T, sR, sC> const &_m);
	template <uint32_t sR, uint32_t sC>
	Matrix<T, R, C> &operator=(Matrix<T, sR, sC> const &_rhs);

	void SetRow(uint32_t _i, Vector<T, C> const &_v);
	void SetColumn(uint32_t _i, Vector<T, R> const &_v);

	union
	{
		std::array<T, R*C>				e;
		std::array<Vector<T, C>, R>		rows;
	};

	Vector<T, C> &operator[](uint32_t _i) { return rows[_i]; }
	Vector<T, C> const &operator[](uint32_t _i) const { return rows[_i]; }

	static Matrix<T, R, C> Identity();
};


using Mat3x3f = Matrix<maths::Decimal, 3, 3>;
using Mat4x4f = Matrix<maths::Decimal, 4, 4>;


// These cannot exist if Matrix's ctor is not constexpr.
//template <typename T, uint32_t R, uint32_t C> 
//struct Zero<Matrix<T, R, C>> { static constexpr Matrix<T, R, C> value = Matrix<T, R, C>(zero<T>); };
//template <typename T, uint32_t R, uint32_t C>
//struct One<Matrix<T, R, C>> { static constexpr Matrix<T, R, C> value = Matrix<T, R, C>(); };//matrix::Identity<T, R, C>(); };


template <typename T, uint32_t R, uint32_t C>
bool operator==(Matrix<T, R, C> const &_lhs, Matrix<T, R, C> const &_rhs);
template <typename T, uint32_t R, uint32_t C>
bool operator!=(Matrix<T, R, C> const &_lhs, Matrix<T, R, C> const &_rhs);

template <typename T, uint32_t R, uint32_t C, uint32_t common>
Matrix<T, R, C> operator*(Matrix<T, R, common> const &_lhs, Matrix<T, common, C> const &_rhs);
template <typename T, uint32_t R, uint32_t C>
Vector<T, R> operator*(Matrix<T, R, C> const &_lhs, Vector<T, C> const &_rhs);

template <typename T, uint32_t R, uint32_t C>
Matrix<T, C, R> FromRows(std::array<Vector<T, C>, R> const &_rows);
template <typename T, uint32_t R, uint32_t C>
Matrix<T, C, R> FromColumns(std::array<Vector<T, R>, C> const &_columns);

template <typename T, uint32_t R, uint32_t C>
Matrix<T, C, R> Transpose(Matrix<T, R, C> const &_v);
template <typename T, uint32_t RC>
Matrix<T, RC, RC> Inverse(Matrix<T, RC, RC> const &_v);

} // namespace maths


#include "maths/matrix.inl"

#endif // __YS_MATRIX_HPP__
