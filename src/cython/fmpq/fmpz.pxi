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

cdef extern from "flint/fmpz.h":

    ctypedef struct fmpz_t:
        pass

    void fmpz_init(fmpz_t f)
    void fmpz_clear(fmpz_t f)
    void fmpz_set_d_2exp(fmpz_t f, double d, slong exp )
    void fmpz_ui_pow_ui(fmpz_t f, ulong g, ulong x )

