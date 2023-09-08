#pragma once

#ifndef _FMPQ_MPOLY_MATRIX_H_
#define _FMPQ_MPOLY_MATRIX_H_

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpq_mpoly.h>

typedef struct fmpq_mpoly_matrix_struct
{
	size_t cols;
	size_t rows;
	fmpq_mpoly_t *cfs;
} fmpq_mpoly_matrix_t[1];

void fmpq_mpoly_matrix_init(fmpq_mpoly_matrix_t A, size_t cols, size_t rows, const fmpq_mpoly_ctx_t ctx );

//multiplication of matrices: cannot be done in place !!!

int fmpq_mpoly_matrix_mul(fmpq_mpoly_matrix_t A, fmpq_mpoly_matrix_t B, fmpq_mpoly_matrix_t C, const fmpq_mpoly_ctx_t ctx );

//pretty print to file

int fmpq_mpoly_matrix_fprint_pretty(FILE *file, const fmpq_mpoly_matrix_t A, const char **x, const fmpq_mpoly_ctx_t ctx );

fmpq_mpoly_struct *fmpq_mpoly_coeff_at(const fmpq_mpoly_matrix_t A, size_t i, size_t j );
//void fmpq_mpoly_coeff_set(const fmpq_mpoly_matrix_t A, fmpq_mpoly_t b, size_t i, size_t j );
int fmpq_mpoly_matrix_coeff_set_str_pretty(const fmpq_mpoly_matrix_t A, const char *str, const char **x, size_t i, size_t j, const fmpq_mpoly_ctx_t ctx );

#endif
