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

from libc.stdlib cimport malloc, free

cdef extern from "Python.h":
    const char* PyUnicode_AsUTF8(object unicode)

cdef extern from "flint/flint.h":
    ctypedef signed long slong

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


    ctypedef struct fmpq_mpoly_t:
        pass

    #I/O

    char *fmpq_mpoly_get_str_pretty(const fmpq_mpoly_t A, const char **x, const fmpq_mpoly_ctx_t ctx ) 
    int fmpq_mpoly_set_str_pretty(fmpq_mpoly_t A, const char *str, const char **x, const fmpq_mpoly_ctx_t ctx )

    #Basic manipulation

    void fmpq_mpoly_gen(fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx )

    #Add/Sub

    void fmpq_mpoly_add(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx )
    void fmpq_mpoly_sub(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx )

    #Mul

    void fmpq_mpoly_mul(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx )

    # Diff/Int
    void fmpq_mpoly_derivative(fmpq_mpoly_t A, const fmpq_mpoly_t B, signed long var, const fmpq_mpoly_ctx_t ctx )
    void fmpq_mpoly_integral(fmpq_mpoly_t A, const fmpq_mpoly_t B, signed long var, const fmpq_mpoly_ctx_t ctx )


cdef class fmpq_mpoly_ctx:

    cdef fmpq_mpoly_ctx_t ctx

    def __cinit__(self, signed long nvars ):
        fmpq_mpoly_ctx_init(self.ctx, nvars, ORD_DEGREVLEX )

    def nvars(self):
        return fmpq_mpoly_ctx_nvars(self.ctx)

    def __dealloc__(self):
        fmpq_mpoly_ctx_clear(self.ctx)


cdef class fmpq_mpoly:

    cdef fmpq_mpoly_t mpoly
    cdef fmpq_mpoly_ctx ctx

    def __cinit__(self, ctx : fmpq_mpoly_ctx ):
        self.ctx = ctx

    def set_str_pretty(self, str, x ):
        py_byte_string = str.encode('UTF-8')
        cdef char *c_string = py_byte_string
        cdef unsigned long x_len = len(x)
        cdef const char **x_ = <char **> malloc(x_len*sizeof(char *))
        cdef i
        for i in range(x_len):
            x_[i] = PyUnicode_AsUTF8(x[i])
        fmpq_mpoly_set_str_pretty(self.mpoly, c_string, x_, self.ctx.ctx )
        free(x_)

    def get_str_pretty(self, x ):
        cdef const char *c_string
        cdef unsigned long x_len = len(x)
        cdef const char **x_ = <char **> malloc(x_len*sizeof(char *))
        cdef i
        for i in range(x_len):
            x_[i] = PyUnicode_AsUTF8(x[i])
        c_string = fmpq_mpoly_get_str_pretty(self.mpoly, x_, self.ctx.ctx )
        string = c_string.decode('UTF-8')
        free(c_string)
        free(x_)
        return string

    def __add__(self, other : fmpq_mpoly ):
        if isinstance(other, fmpq_mpoly ):
            if( self.ctx != other.ctx ):
                raise ValueError("cannot __add__ : first and second mpoly have different context")
            else:
                ret = fmpq_mpoly(self.ctx)
                fmpq_mpoly_add(ret.mpoly, self.mpoly, other.mpoly, self.ctx.ctx )
                return ret
        else:
            raise TypeError("Unsupported operand type")

    def __sub__(self, other : fmpq_mpoly ):
        if isinstance(other, fmpq_mpoly ):
            if( self.ctx != other.ctx ):
                raise ValueError("cannot __sub__ : first and second mpoly have different context")
            else:
                ret = fmpq_mpoly(self.ctx)
                fmpq_mpoly_sub(ret.mpoly, self.mpoly, other.mpoly, self.ctx.ctx )
                return ret
        else:
            raise TypeError("Unsupported operand type")

    def __mul__(self, other : fmpq_mpoly ):
        if isinstance(other, fmpq_mpoly ):
            if( self.ctx != other.ctx ):
                raise ValueError("cannot __mul__ : first and second mpoly have different context")
            else:
                ret = fmpq_mpoly(self.ctx)
                fmpq_mpoly_mul(ret.mpoly, self.mpoly, other.mpoly, self.ctx.ctx )
                return ret
        else:
            raise TypeError("Unsupported operand type")

    def der(self, var : long ):
        ret = fmpq_mpoly(self.ctx)
        fmpq_mpoly_derivative(ret.mpoly, self.mpoly, var, self.ctx.ctx )
        return ret

    def integ(self, var : long ):
        ret = fmpq_mpoly(self.ctx)
        fmpq_mpoly_integral(ret.mpoly, self.mpoly, var, self.ctx.ctx )
        return ret

    def gen(self, var : long ):
        fmpq_mpoly_gen(self.mpoly, var, self.ctx.ctx )


#next we'll work on fmpq_mpoly_tensor ? (Madness)

