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

include "msolve.pxi"

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def solve_from_gens(x):
    cdef unsigned long x_len
    cdef fmpq_mpoly_ctx ctx
    cdef fmpq_mpoly_t *x_
    cdef msolve_re_solutions_t sols
    cdef char **fake_vnames
    cdef double *sols_buffer
    if isinstance(x, collections.abc.Collection ):
        x_len = len(x)
        fake_vnames = <char **> malloc(x_len*sizeof(char *))
        for i in range(x_len):
            fake_vnames[i] = NULL
    else:
        raise ValueError("cannot solve_from_gens: arr needs to be a collection")
    if(x_len == 0):
        raise ValueError("the collection cannot be empty")
    if isinstance(x[0], fmpq_mpoly ):
        ctx = (<fmpq_mpoly>x[0]).ctx
    else:
        raise ValueError("the collection needs to be made of fmpq_mpoly s with the same ctx") 
    x_ = <fmpq_mpoly_t *> malloc(sizeof(fmpq_mpoly_t)*x_len)
    for i in range(x_len):
        if ctx != (<fmpq_mpoly>x[i]).ctx:
            raise ValueError("the collection needs to be made of fmpq_mpoly s with the same ctx")
         #x_[i] = (<fmpq_mpoly>x[i]).mpoly
        fmpq_mpoly_init(x_[i], ctx.ctx )
        fmpq_mpoly_set(x_[i], (<fmpq_mpoly>x[i]).mpoly, ctx.ctx )
        #c_string = fmpq_mpoly_get_str_pretty(x_[i], NULL, ctx.ctx )
        #print(c_string.decode('UTF-8'))
#
#    #print('sigh')
    msolve_from_fmpq_mpolys(sols, x_, x_len, fake_vnames, ctx.ctx )
    sols_buffer = get_sols_buffer(sols)
    result = numpy.copy(numpy.asarray(<double [:sols.nb_real_roots*ctx.nvars()]> sols_buffer).reshape(sols.nb_real_roots,ctx.nvars()))
    #done !!! now clear everything (except the final result)
    msolve_solutions_clear(sols)
    free(sols_buffer)
    free(fake_vnames)
    #you also need to free the copies of the generators
    for i in range(x_len):
        fmpq_mpoly_clear(x_[i], ctx.ctx )
    free(x_)
    return result
