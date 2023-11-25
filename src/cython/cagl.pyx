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

include "c_utils.pxi"

#flint basic includes

include "flint.pxi"
include "fmpq/fmpz.pxi"
include "fmpq/fmpq.pxi"

#then fmpq_poly and fmpq_mpoly_matrix
include "fmpq/fmpq_mpoly.pyx"
include "fmpq/fmpq_mpoly_matrix.pyx"

#include msolve funcionality
include "msolve/msolve.pyx"