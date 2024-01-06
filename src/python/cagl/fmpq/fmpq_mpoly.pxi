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

cdef extern from "flint/fmpq_mpoly.h":

    ctypedef enum ordering_t:
        ORD_LEX
        ORD_DEGLEX
        ORD_DEGREVLEX

    ctypedef struct fmpq_mpoly_ctx_t:
        pass

    void fmpq_mpoly_ctx_init(fmpq_mpoly_ctx_t ctx, signed long nvars, const ordering_t ord )
    signed long fmpq_mpoly_ctx_nvars(fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_ctx_clear(fmpq_mpoly_ctx_t ctx)

    ctypedef struct fmpq_mpoly_struct:
        pass

    ctypedef fmpq_mpoly_struct fmpq_mpoly_t[1]

    #mem manag

    void fmpq_mpoly_init(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx )
    void fmpq_mpoly_clear(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx )

    #I/O

    char *fmpq_mpoly_get_str_pretty(const fmpq_mpoly_t A, const char **x, const fmpq_mpoly_ctx_t ctx ) 
    int fmpq_mpoly_set_str_pretty(fmpq_mpoly_t A, const char *str, const char **x, const fmpq_mpoly_ctx_t ctx )

    #Basic manipulation

    void fmpq_mpoly_gen(fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx )
    void fmpq_mpoly_set(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx )

    #Constants
    void fmpq_mpoly_set_fmpq(fmpq_mpoly_t A, const fmpq_t c, const fmpq_mpoly_ctx_t ctx )

    #Add/Sub

    void fmpq_mpoly_add(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx )
    void fmpq_mpoly_sub(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx )

    #Mul

    void fmpq_mpoly_mul(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx )

    # Diff/Int
    void fmpq_mpoly_derivative(fmpq_mpoly_t A, const fmpq_mpoly_t B, signed long var, const fmpq_mpoly_ctx_t ctx )
    void fmpq_mpoly_integral(fmpq_mpoly_t A, const fmpq_mpoly_t B, signed long var, const fmpq_mpoly_ctx_t ctx )
