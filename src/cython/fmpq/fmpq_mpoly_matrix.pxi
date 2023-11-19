''' This file is part of cagl.
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
 * Lorenzo Magherini (m4gh3) '''

# cagl.pyx

include "fmpq_mpoly.pxi"

cdef extern from "../include/fmpq_mpoly_matrix.h":

    ctypedef struct fmpq_mpoly_matrix_struct:
        size_t cols
        size_t rows
        fmpq_mpoly_t *cfs

    ctypedef fmpq_mpoly_matrix_struct fmpq_mpoly_matrix_t[1]

    void fmpq_mpoly_matrix_init(fmpq_mpoly_matrix_t A, size_t rows, size_t cols, const fmpq_mpoly_ctx_t ctx )
    void fmpq_mpoly_matrix_clear(fmpq_mpoly_matrix_t A, const fmpq_mpoly_ctx_t ctx )
    int fmpq_mpoly_matrix_mul(fmpq_mpoly_matrix_t A, fmpq_mpoly_matrix_t B, fmpq_mpoly_matrix_t C, const fmpq_mpoly_ctx_t ctx )
    int fmpq_mpoly_matrix_hadamard(fmpq_mpoly_matrix_t A, fmpq_mpoly_matrix_t B, fmpq_mpoly_matrix_t C, const fmpq_mpoly_ctx_t ctx )
    int fmpq_mpoly_matrix_add(fmpq_mpoly_matrix_t A, fmpq_mpoly_matrix_t B, fmpq_mpoly_matrix_t C, const fmpq_mpoly_ctx_t ctx )
    int fmpq_mpoly_matrix_sub(fmpq_mpoly_matrix_t A, fmpq_mpoly_matrix_t B, fmpq_mpoly_matrix_t C, const fmpq_mpoly_ctx_t ctx )
    void fmpq_mpoly_matrix_squared_frobenius(fmpq_mpoly_t a, const fmpq_mpoly_matrix_t B, const fmpq_mpoly_ctx_t ctx )
    fmpq_mpoly_struct *fmpq_mpoly_coeff_at(const fmpq_mpoly_matrix_t A, size_t i, size_t j )
    int fmpq_mpoly_matrix_gens_fill(const fmpq_mpoly_matrix_t A, size_t g0, const fmpq_mpoly_ctx_t ctx )


cimport cython
import numpy
cimport numpy
cimport libc.math