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
