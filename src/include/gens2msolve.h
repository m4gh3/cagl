/* This file is part of cagl.
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
 * Lorenzo Magherini (m4gh3) */

#include <stdio.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpq_mpoly.h>

#include "../msolve/src/fglm/data_fglm.h"
#include "../msolve/src/usolve/data_usolve.h"
#include "../msolve/src/msolve/msolve.h"


data_gens_ff_t *gens2msolve(fmpq_mpoly_t *polys, size_t n, fmpq_mpoly_ctx_t ctx );
void free_data_gens(data_gens_ff_t *gens);
