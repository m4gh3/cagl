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

include "fmpq_mpoly_matrix.pxi"


cdef class fmpq_mpoly_matrix:

    cdef fmpq_mpoly_matrix_t mpoly_mat
    cdef fmpq_mpoly_ctx ctx

    def __cinit__(self, rows, cols, ctx ):#: fmpq_mpoly_ctx ):
        self.ctx = ctx
        cdef size_t rows_
        cdef size_t cols_
        rows_ = PyLong_AsSize_t(rows)
        cols_ = PyLong_AsSize_t(cols)
        fmpq_mpoly_matrix_init(self.mpoly_mat, rows_, cols_, (<fmpq_mpoly_ctx>ctx).ctx )

    def __dealloc__(self):
        fmpq_mpoly_matrix_clear(self.mpoly_mat, self.ctx.ctx )

    def __add__(self, other : fmpq_mpoly_matrix ):
        if isinstance(other, fmpq_mpoly_matrix ):
            if( self.ctx != other.ctx ):
                raise ValueError("cannot __add__ : first and second mpoly_matrix have different context")
            elif((self.mpoly_mat.rows != other.mpoly_mat.rows) or (self.mpoly_mat.cols != other.mpoly_mat.cols)):
                raise ValueError("cannot __add__: first and second mpoly_matrix have different shapes")
            else:
                ret = fmpq_mpoly_matrix(self.mpoly_mat.rows, self.mpoly_mat.cols, self.ctx )
                fmpq_mpoly_matrix_add(ret.mpoly_mat, self.mpoly_mat, other.mpoly_mat, self.ctx.ctx )
                return ret
        else:
            raise TypeError("Unsupported operand type")

    def __sub__(self, other : fmpq_mpoly_matrix ):
        if isinstance(other, fmpq_mpoly_matrix ):
            if( self.ctx != other.ctx ):
                raise ValueError("cannot __sub__ : first and second mpoly_matrix have different context")
            elif((self.mpoly_mat.rows != other.mpoly_mat.rows) or (self.mpoly_mat.cols != other.mpoly_mat.cols)):
                raise ValueError("cannot __sub__: first and second mpoly_matrix have different shapes")
            else:
                ret = fmpq_mpoly_matrix(self.mpoly_mat.rows, self.mpoly_mat.cols, self.ctx )
                fmpq_mpoly_matrix_sub(ret.mpoly_mat, self.mpoly_mat, other.mpoly_mat, self.ctx.ctx )
                return ret
        else:
            raise TypeError("Unsupported operand type")

    def __mul__(self, other : fmpq_mpoly_matrix ):
        if isinstance(other, fmpq_mpoly_matrix ):
            if( self.ctx != other.ctx ):
                raise ValueError("cannot __mul__ : first and second mpoly_matrix have different context")
            elif((self.mpoly_mat.rows != other.mpoly_mat.rows) or (self.mpoly_mat.cols != other.mpoly_mat.cols)):
                raise ValueError("cannot __mul__: first and second mpoly_matrix have different shapes")
            else:
                ret = fmpq_mpoly_matrix(self.mpoly_mat.rows, self.mpoly_mat.cols, self.ctx )
                fmpq_mpoly_matrix_hadamard(ret.mpoly_mat, self.mpoly_mat, other.mpoly_mat, self.ctx.ctx )
                return ret
        else:
            raise TypeError("Unsupported operand type")

    def __matmul__(self, other : fmpq_mpoly_matrix ):
        if isinstance(other, fmpq_mpoly_matrix ):
            if( self.ctx != other.ctx ):
                raise ValueError("cannot __matmul__ : first and second mpoly_matrix have different context")
            elif((self.mpoly_mat.cols != other.mpoly_mat.rows)):
                raise ValueError("cannot __matmul__: first and second mpoly_matrix have incompatible shapes")
            else:
                ret = fmpq_mpoly_matrix(self.mpoly_mat.rows, other.mpoly_mat.cols, self.ctx )
                fmpq_mpoly_matrix_mul(ret.mpoly_mat, self.mpoly_mat, other.mpoly_mat, self.ctx.ctx )
                return ret
        else:
            raise TypeError("Unsupported operand type")

    #void fmpq_mpoly_matrix_squared_frobenius(fmpq_mpoly_t a, const fmpq_mpoly_matrix_t B, const fmpq_mpoly_ctx_t ctx )
    def sqfrob(self):
        ret = fmpq_mpoly(self.ctx)
        fmpq_mpoly_matrix_squared_frobenius(ret.mpoly, self.mpoly_mat, self.ctx.ctx )
        return ret

    def __getitem__(self, key ):

        cdef size_t i, j
        cdef fmpq_mpoly_struct *coeff

        i = PyLong_AsSize_t(key[0])
        j = PyLong_AsSize_t(key[1])
        if( (i >= self.mpoly_mat.rows) or (j >= self.mpoly_mat.cols) ):
            raise ValueError("cannot __getitem__: index out of range")
        else:
            coeff = fmpq_mpoly_coeff_at(self.mpoly_mat, i, j )
            ret = fmpq_mpoly(self.ctx)
            fmpq_mpoly_set(ret.mpoly, coeff, self.ctx.ctx )
            return ret

    def __setitem__(self, key, value : fmpq_mpoly ):

        cdef size_t i, j
        cdef fmpq_mpoly_struct *coeff

        i = PyLong_AsSize_t(key[0])
        j = PyLong_AsSize_t(key[1])

        if( self.ctx != value.ctx ):
            raise ValueError("cannot __setitem__: the mpoly_matrix and the mpoly do not have the same context")
        elif( (i >= self.mpoly_mat.rows) or (j >= self.mpoly_mat.cols) ):
            raise ValueError("cannot __setitem__: index out of range")
        else:
            coeff = fmpq_mpoly_coeff_at(self.mpoly_mat, i, j )
            fmpq_mpoly_set(coeff, value.mpoly, self.ctx.ctx )

    #int fmpq_mpoly_matrix_gens_fill(const fmpq_mpoly_matrix_t A, size_t g0, const fmpq_mpoly_ctx_t ctx )

    def gens_fill(self, g0 ):
        if fmpq_mpoly_matrix_gens_fill(self.mpoly_mat, g0, self.ctx.ctx ) < 0:
            raise ValueError(f"cannot gens_fill : there are not enough generators to fill the matrix starting from {g0}")

    def get_str_pretty(self, x=None ):
        cdef const char *c_string
        cdef unsigned long x_len
        cdef const char **x_
        cdef long k
        cdef long i
        cdef long j
        string = ""
        if isinstance(x, collections.abc.Iterable ):
            x_len = len(x)
            x_ = <char **> malloc(x_len*sizeof(char *))
            for k in range(x_len):
                if isinstance(x[k], str ):
                    x_[k] = PyUnicode_AsUTF8(x[k])
                else:
                    free(x_)
                    raise ValueError(f"cannot get_str_pretty : element of x at index={k} is not a string" )
            string += "[\n"
            for i in range(self.mpoly_mat.rows):
                string += "\t[ "
                for j in range(self.mpoly_mat.cols):
                    c_string = fmpq_mpoly_get_str_pretty(self.mpoly_mat.cfs[i*self.mpoly_mat.cols+j], x_, self.ctx.ctx )
                    string += (c_string.decode('UTF-8') + ", ")
                    free(c_string)
                string += "], \n"
            string += "]"
            free(x_)
        elif x == None:
            string += "[\n"
            for i in range(self.mpoly_mat.rows):
                string += "\t[ "
                for j in range(self.mpoly_mat.cols):
                    c_string = fmpq_mpoly_get_str_pretty(self.mpoly_mat.cfs[i*self.mpoly_mat.cols+j], NULL, self.ctx.ctx )
                    string += (c_string.decode('UTF-8') + ", ")
                    free(c_string)
                string += "], \n"
            string += "]"
        else:
            raise ValueError("cannot get_str_pretty : x needs to be either an iterable (with integer indexes) containing strings or None")

        return string

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def set_from_np(self, arr : float [:,:], ee : int ):
        cdef long i, j, ek
        cdef float v, b
        cdef fmpq_t res
        cdef fmpz_t p, q
        fmpq_init(res)
        fmpz_init(p)
        fmpz_init(q)
        if ( arr.shape[0] != self.mpoly_mat.rows ) or ( arr.shape[1] != self.mpoly_mat.cols ):
            raise ValueError("cannot set_from_np : fmpq_mpoly_matrix and array shapes do not match")
        if ( ee <= 0 ):
            raise ValueError("cannot set_from_np : ee <= 0 and needs to be > 0")
        
        for i in range(self.mpoly_mat.rows):
            for j in range(self.mpoly_mat.cols):
                v = arr[i,j]
                if v != 0 :
                    ek = ee-libc.math.floor(libc.math.log2(libc.math.fabs(v)))
                    ek = 0 if ek < 0 else ek
                    fmpz_set_d_2exp(p, v, ek )
                    #fmpz_ui_pow_ui(q, 2, ek ) #it doesn't link yet probably need to wait until flint3
                    fmpz_set_d_2exp(q, 1, ek )
                    fmpq_set_fmpz_frac(res, p, q )
                    fmpq_mpoly_set_fmpq(self.mpoly_mat.cfs[i*self.mpoly_mat.cols+j], res, self.ctx.ctx )
        fmpq_clear(res)
        fmpz_clear(p)
        fmpz_clear(q)

cdef extern from "../msolve/src/fglm/data_fglm.h":
    pass
cdef extern from "../msolve/src/msolve/msolve-data.h":
    pass
cdef extern from "../msolve/src/msolve/msolve.h":
    ctypedef msolve_re_solutions_t
    void msolve_from_fmpq_mpolys(
	msolve_re_solutions_t sols,
        fmpq_mpoly_t *polys,
        size_t n,
        char **vnames,
        fmpq_mpoly_ctx_t ctx
	)
    cdef void msolve_solutions_clear(msolve_re_solutions_t sols)
    cdef double *get_sols_buffer(msolve_re_solutions_t sols)

def solve_from_gens(arr):
    pass
