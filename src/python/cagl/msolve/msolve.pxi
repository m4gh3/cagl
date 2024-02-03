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

include "../fmpq/fmpq_mpoly.pxi"
cimport cython

#cdef extern from "../../msolve/src/fglm/data_fglm.h":
#    pass
cdef extern from "../../msolve/src/msolve/msolve-data.h":
    struct msolve_re_solutions_struct:
        long nb_real_roots
        void *real_roots
        void *real_pts
cdef extern from "../../msolve/src/msolve/msolve.h":
    ctypedef msolve_re_solutions_struct msolve_re_solutions_t[1]
    void msolve_from_fmpq_mpolys(
	msolve_re_solutions_t sols,
        fmpq_mpoly_t *polys,
        size_t n,
        char **vnames,
        fmpq_mpoly_ctx_t ctx
	)
    cdef void msolve_solutions_clear(msolve_re_solutions_t sols)
    cdef double *get_sols_buffer(msolve_re_solutions_t sols)

