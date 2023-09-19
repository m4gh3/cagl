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

#pragma once

#ifndef _FMPQ_MPOLY_NDARRAY_H_
#define _FMPQ_MPOLY_NDARRAY_H_

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpq_mpoly.h>

typedef struct fmpq_mpoly_ndarray_struct
{
	size_t ndim;
	size_t nelem;
	size_t shape_[6];
	size_t strides_[6];
	size_t *shape;
	size_t *strides;
	fmpq_mpoly_t *cfs;
} fmpq_mpoly_ndarray_t[1];

void fmpq_mpoly_ndarray_init(fmpq_mpoly_ndarray_t A, size_t *shape, size_t *strides, const fmpq_mpoly_ctx_t ctx );

//multiplication of matrices: cannot be done in place !!!

int fmpq_mpoly_ndarray_mat_mul(fmpq_mpoly_ndarray_t A, fmpq_mpoly_ndarray_t B, fmpq_mpoly_ndarray_t C, const fmpq_mpoly_ctx_t ctx );

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
