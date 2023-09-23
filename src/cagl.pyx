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
import collections.abc

cdef extern from "Python.h":
    const char* PyUnicode_AsUTF8(object unicode)
    size_t PyLong_AsSize_t(object pylong)

cdef extern from "flint/flint.h":
    ctypedef signed long slong
    ctypedef unsigned long ulong

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

    def __cinit__(self, ctx : fmpq_mpoly_ctx, string=None, x=None ):
        self.ctx = ctx
        fmpq_mpoly_init(self.mpoly, ctx.ctx )
        if isinstance(string, str ):
            self.set_str_pretty(string, x )

    def __dealloc__(self):
        fmpq_mpoly_clear(self.mpoly, self.ctx.ctx )

    def set_str_pretty(self, string : str, x=None ):
        cdef const char *c_string
        cdef unsigned long x_len
        cdef const char **x_
        cdef i
        c_string = PyUnicode_AsUTF8(string)
        if isinstance(x, collections.abc.Iterable ):
            x_len = len(x)
            x_ = <char **> malloc(x_len*sizeof(char *))
            for i in range(x_len):
                if isinstance(x[i], str ):
                    x_[i] = PyUnicode_AsUTF8(x[i])
                else:
                    free(x_)
                    raise ValueError(f"cannot set_str_pretty : element of x at index={i} is not a string" )
            fmpq_mpoly_set_str_pretty(self.mpoly, c_string, x_, self.ctx.ctx )
            free(x_)
        elif x == None:
            fmpq_mpoly_set_str_pretty(self.mpoly, c_string, NULL, self.ctx.ctx )
        else:
            raise ValueError("cannot set_str_pretty : x needs to be either an iterable (with integer indexes) containing strings or None")

        return string

    def get_str_pretty(self, x=None ):
        cdef const char *c_string
        cdef unsigned long x_len
        cdef const char **x_
        cdef i
        if isinstance(x, collections.abc.Iterable ):
            x_len = len(x)
            x_ = <char **> malloc(x_len*sizeof(char *))
            for i in range(x_len):
                if isinstance(x[i], str ):
                    x_[i] = PyUnicode_AsUTF8(x[i])
                else:
                    free(x_)
                    raise ValueError(f"cannot get_str_pretty : element of x at index={i} is not a string" )
            c_string = fmpq_mpoly_get_str_pretty(self.mpoly, x_, self.ctx.ctx )
            string = c_string.decode('UTF-8')
            free(c_string)
            free(x_)
        elif x == None:
            c_string = fmpq_mpoly_get_str_pretty(self.mpoly, NULL, self.ctx.ctx )
            string = c_string.decode('UTF-8')
            free(c_string)
        else:
            raise ValueError("cannot get_str_pretty : x needs to be either an iterable (with integer indexes) containing strings or None")

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

cdef extern from "./include/fmpq_mpoly_matrix.h":

    ctypedef struct fmpq_mpoly_matrix_struct:
        size_t cols;
        size_t rows;
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


cdef extern from "flint/fmpz.h":

    ctypedef struct fmpz_t:
        pass

    void fmpz_init(fmpz_t f)
    void fmpz_clear(fmpz_t f)
    void fmpz_set_d_2exp(fmpz_t f, double d, slong exp )
    void fmpz_ui_pow_ui(fmpz_t f, ulong g, ulong x )


cdef extern from "flint/fmpq.h":

    ctypedef struct fmpq_t:
        pass

    void fmpq_init(fmpq_t x)
    void fmpq_clear(fmpq_t x)
    void fmpq_set_fmpz_frac(fmpq_t res, const fmpz_t p, const fmpz_t q )



cimport cython
import numpy
cimport numpy
cimport libc.math

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
