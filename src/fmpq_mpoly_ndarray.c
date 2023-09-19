/* This file is part of cagl.
 *
 * cagl is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * cagl is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cagl.  If not, see <https://www.gnu.org/licenses/>
 *
 * Authors:
 * Lorenzo Magherini (m4gh3) */

#include "include/fmpq_mpoly_ndarray.h"


void fmpq_mpoly_ndarray_init(fmpq_mpoly_ndarray_t A, size_t *shape, size_t ndim, size_t const fmpq_mpoly_ctx_t ctx )
{

	A->ndim = ndim;

	if( ndim < 7 )
	{
		A->shape = A->shape_;
		A->strides = A->strides_;
	}
	else
	{
		A->shape = (size_t *)malloc(sizeof(size_t)*ndim);
		A->strides = (size_t *)malloc(sizeof(size_t)*ndim);
	}

	size_t s = 1;
	for(size_t i=ndim; i > 0; i-- )
	{
		A->strides[i-1] = s;
		s *= shape[i-1];
		A->shape[i-1] = s;
	}

	A->nelem = s;

}

//multiplication of matrices: cannot be done in place !!!

int fmpq_mpoly_ndarray_mat_mul(fmpq_mpoly_ndarray_t A, fmpq_mpoly_ndarray_t B, fmpq_mpoly_ndarray_t C, const fmpq_mpoly_ctx_t ctx )
{

	mpq_mpoly_t temp;
	fmpq_mpoly_init(temp, ctx );

	if(B->ndim == 2)
	{
		if(C->ndim == 2)
		{
			for(size_t i_a=0,i_b=0; i_a < A->shape[0]; i_a += A->strides[0], i_b += B->strides[0]; )
			{
				for(size_t k_a=0, k_c=0; k_a < A->shape[1]; k_a += A->strides[1], k_c += C->strides[1]; )
				{
					fmpq_mpoly_zero(A->cfs[i_a+k_a], ctx );
					for(size_t j_b=0, j_c=0; j_b < B->shape[1]; j_b += B->strides[1], j_c += C->strides[0] )
					{
						fmpq_mpoly_mul(temp, B->cfs[b_i+b_j], C->cfs[c_j+c_k], ctx );
						fmpq_mpoly_add(A->cfs[i_a+k_a], A->cfs[i_a+k_a], temp, ctx );
					}
				}
			}
		}
		else if(C->ndim == 3)
		{
			for(size_t l_a=0,l_c=0; l_a < A->shape[0]; i_a += A->strides[0], c += C->strides[0] )
			{
				for(size_t i_a=0,i_b=0; i_a < A->shape[1]; i_a += A->strides[1], i_b += B->strides[0]; )
				{
					for(size_t k_a=0, k_c=0; k_a < A->shape[2]; k_a += A->strides[2], k_c += C->strides[2]; )
					{
						fmpq_mpoly_zero(A->cfs[l_a+i_a+k_a], ctx );
						for(size_t j_b=0, j_c=0; j_b < B->shape[1]; j_b += B->strides[1], j_c += C->strides[1] )
						{
							fmpq_mpoly_mul(temp, B->cfs[i_b+j_b], C->cfs[l_c+j_c+k_c], ctx );
							fmpq_mpoly_add(A->cfs[l_a+i_a+k_a], A->cfs[l_a+i_a+k_a], temp, ctx );
						}
					}
				}
			}
		}
	}
	else if(B->ndim == 3 )
	{
		if(C->ndim == 2)
		{
		}
		else if(C->ndim == 3)
		{
		}
	}

}

//elementwise product:

int fmpq_mpoly_ndarray_mul(fmpq_mpoly_ndarray_t A, fmpq_mpoly_ndarray_t B, fmpq_mpoly_ndarray_t C, const fmpq_mpoly_ctx_t ctx );

//matrix add:

int fmpq_mpoly_ndarray_add(fmpq_mpoly_ndarray_t A, fmpq_mpoly_ndarray_t B, fmpq_mpoly_ndarray_t C, const fmpq_mpoly_ctx_t ctx );

//matrix sub:

int fmpq_mpoly_ndarray_sub(fmpq_mpoly_ndarray_t A, fmpq_mpoly_ndarray_t B, fmpq_mpoly_ndarray_t C, const fmpq_mpoly_ctx_t ctx );

//squared frobenius norm:

void fmpq_mpoly_ndarray_sqfrob(fmpq_mpoly_t a, const fmpq_mpoly_ndarray_t B, const fmpq_mpoly_ctx_t ctx );

//pretty print to file

int fmpq_mpoly_matrix_fprint_pretty(FILE *file, const fmpq_mpoly_ndarray_t A, const char **x, const fmpq_mpoly_ctx_t ctx );

fmpq_mpoly_struct *fmpq_mpoly_coeff_at(const fmpq_mpoly_ndarray_t A, size_t i, size_t j );
//void fmpq_mpoly_coeff_set(const fmpq_mpoly_ndarray_t A, fmpq_mpoly_t b, size_t i, size_t j );
int fmpq_mpoly_matrix_coeff_set_str_pretty(const fmpq_mpoly_ndarray_t A, const char *str, const char **x, size_t i, size_t j, const fmpq_mpoly_ctx_t ctx );

#endif
