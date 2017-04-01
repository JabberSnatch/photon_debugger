#ifndef __YS_MATRIX_HPP__
#define __YS_MATRIX_HPP__

#include "vector.h"
#include "maths.h"


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


template <typename T, uint32_t R, uint32_t C>
struct Matrix
{
	// NOTE: std::array::fill isn't constexpr, so Matrix's ctor won't be either.
	constexpr Matrix() :
		Matrix(zero<T>)
	{
		for (uint32_t i = 0; i < Min(R, C); ++i)
			e[i * C + i] = one<T>;
	}
	// TODO: Profile a use case when replacing e.fill with a loop.
	explicit constexpr Matrix(T _v) { e.fill(_v); }
	constexpr Matrix(std::initializer_list<T> _args) {
		std::copy(_args.begin(), _args.end(), e.begin());
	}

	union
	{
		std::array<T, R*C>				e;
		std::array<Vector<T, C>, R>		rows;
	};

	constexpr Vector<T, C> &operator[](uint32_t _i) { return rows[_i]; }
	constexpr Vector<T, C> const &operator[](uint32_t _i) const { return rows[_i]; }
};


// These cannot exist if Matrix's ctor is not constexpr.
//template <typename T, uint32_t R, uint32_t C> 
//struct Zero<Matrix<T, R, C>> { static constexpr Matrix<T, R, C> value = Matrix<T, R, C>(zero<T>); };
//template <typename T, uint32_t R, uint32_t C>
//struct One<Matrix<T, R, C>> { static constexpr Matrix<T, R, C> value = matrix::Identity<T, R, C>(); };


template <typename T, uint32_t R, uint32_t C, uint32_t common>
constexpr Matrix<T, R, C> operator*(Matrix<T, R, common> const &_lhs, Matrix<T, common, C> const &_rhs);
template <typename T, uint32_t R, uint32_t C>
constexpr Vector<T, R> operator*(Matrix<T, R, C> const &_lhs, Vector<T, C> const &_rhs);


namespace matrix
{

// Not constexpr unless Matrix::ctor is.
template <typename T, uint32_t R, uint32_t C>
constexpr Matrix<T, R, C> Identity();

template <typename T, uint32_t R, uint32_t C>
constexpr Matrix<T, C, R> Transpose(Matrix<T, R, C> const &_v);
template <typename T, uint32_t RC>
Matrix<T, RC, RC> Inverse(Matrix<T, RC, RC> const &_v);

} // namespace matrix
} // namespace maths


#include "matrix.inl"

#endif // __YS_MATRIX_HPP__
