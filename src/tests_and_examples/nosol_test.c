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

/* Note: this code takes also some inspiration from msolve's source code */


#include <stdio.h>
#include "../msolve/src/fglm/data_fglm.h"
//#include "../msolve/src/usolve/data_usolve.h"
#include "../msolve/src/msolve/msolve.h"
#include "../include/fmpq_mpoly_matrix.h"

int main(int argc, char **argv)
{

	fmpq_mpoly_ctx_t ctx;
	
	fmpq_mpoly_ctx_init(ctx, 2, ORD_DEGREVLEX );
	
	const char *ring_gen_names[] = {"x", "y" };
			
	fmpq_mpoly_t f[2];

	fmpq_mpoly_init(f[0], ctx );
	fmpq_mpoly_init(f[1], ctx );

	fmpq_mpoly_set_str_pretty(f[0], "x*y - 1", ring_gen_names, ctx );
	fmpq_mpoly_set_str_pretty(f[1], "x^2 + y^2 + 1", ring_gen_names, ctx );

	msolve_re_solutions_t sols;
	msolve_from_fmpq_mpolys(sols, f, 2, ring_gen_names, ctx );

	//now printing solutions
	printf("number of solutions found: %d\n", sols->nb_real_roots );

	printf("printing from buffer:\n");
	
	double *sols_buf = get_sols_buffer(sols);
	size_t nvars = fmpq_mpoly_ctx_nvars(ctx);

	for(size_t i=0; i < nvars*sols->nb_real_roots; i+=nvars )
	{
		for(size_t j=0; j < nvars; j++ )
			printf("%f ,", sols_buf[i+j] );
		putchar('\n');
	}

	printf("---------------------");

	msolve_solutions_clear(sols);

}
