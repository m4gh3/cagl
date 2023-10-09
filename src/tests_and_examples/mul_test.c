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

double *get_sols_buffer(msolve_re_solutions_t sols )
{
	long nb_real_roots = sols->nb_real_roots;
	long nvars = sols->real_pts[0]->nvars;
	double *buffer = malloc(sols->nb_real_roots*nvars*sizeof(double));
	double *buf_ptr = buffer;
	//convert all the coords to float
	for(size_t i=0; i < nb_real_roots; i++ )
		for(size_t j=0; j < nvars; j++ )
			*(buf_ptr++) = mpz_get_d(sols->real_pts[i]->coords[j]->val_do)/pow(2, sols->real_pts[i]->coords[j]->k_do ); //it's a bit wrong to just use val_do and k_do
	return buffer;
}

int main(int argc, char **argv)
{

	fmpq_mpoly_ctx_t ctx;
	
	fmpq_mpoly_ctx_init(ctx, 6, ORD_DEGREVLEX );
	
	const char *ring_gen_names[] = {"w0", "w1", "w2", "w3", "l0", "l1" };
	
	/*fmpq_mpoly_t P[2];
	
	fmpq_mpoly_init(P[0], ctx );
	fmpq_mpoly_init(P[1], ctx );
	
	fmpq_mpoly_set_str_pretty(P[0], "w0^2+w1^2-1", ring_gen_names, ctx );
	fmpq_mpoly_set_str_pretty(P[1], "w0-w1", ring_gen_names, ctx );*/
		
	fmpq_mpoly_matrix_t A, B, x_train, y_train, y_hat0, y_hat1, y_hat, E_vec;
	fmpq_mpoly_t E, constraint, L_grad[6];

	fmpq_mpoly_matrix_init(A, 1, 2, ctx );
	fmpq_mpoly_matrix_init(B, 1, 2, ctx );
	fmpq_mpoly_matrix_init(x_train, 2, 5, ctx );
	fmpq_mpoly_matrix_init(y_train, 1, 5, ctx );
	fmpq_mpoly_matrix_init(y_hat0, 1, 5, ctx );
	fmpq_mpoly_matrix_init(y_hat1, 1, 5, ctx );
	fmpq_mpoly_matrix_init(y_hat, 1, 5, ctx );
	fmpq_mpoly_matrix_init(E_vec, 1, 5, ctx );
	fmpq_mpoly_init(E, ctx );
	fmpq_mpoly_init(constraint, ctx );
	for(int i=0; i < 6; i++ )
		fmpq_mpoly_init(L_grad[i], ctx );


	fmpq_mpoly_matrix_coeff_set_str_pretty(A, "w0", ring_gen_names, 0, 0, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(A, "w1", ring_gen_names, 0, 1, ctx );

	fmpq_mpoly_matrix_coeff_set_str_pretty(B, "w2", ring_gen_names, 0, 0, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(B, "w3", ring_gen_names, 0, 1, ctx );

	fmpq_mpoly_set_str_pretty(constraint, "l0*w0^2+l0*w1^2-l0+l1*w2^2+l1*w3^2-l1", ring_gen_names, ctx );

	fmpq_mpoly_matrix_coeff_set_str_pretty(x_train, "1", ring_gen_names, 0, 0, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(x_train, "1", ring_gen_names, 1, 0, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(y_train, "1", ring_gen_names, 0, 0, ctx );
	
	fmpq_mpoly_matrix_coeff_set_str_pretty(x_train, "1", ring_gen_names, 0, 1, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(x_train, "-1", ring_gen_names, 1, 1, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(y_train, "-1", ring_gen_names, 0, 1, ctx );
	
	fmpq_mpoly_matrix_coeff_set_str_pretty(x_train, "-1", ring_gen_names, 0, 2, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(x_train, "1", ring_gen_names, 1, 2, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(y_train, "-1", ring_gen_names, 0, 2, ctx );

	
	fmpq_mpoly_matrix_coeff_set_str_pretty(x_train, "-1", ring_gen_names, 0, 3, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(x_train, "-1", ring_gen_names, 1, 3, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(y_train, "1", ring_gen_names, 0, 3, ctx );

	/*fmpq_mpoly_matrix_coeff_set_str_pretty(x_train, "1/2", ring_gen_names, 0, 4, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(x_train, "1/2", ring_gen_names, 1, 4, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(y_train, "1/4", ring_gen_names, 0, 4, ctx );*/


	fmpq_mpoly_matrix_fprint_pretty(stdout, A, ring_gen_names, ctx );
	printf("------------\n");
	fmpq_mpoly_matrix_fprint_pretty(stdout, B, ring_gen_names, ctx );
	printf("------------\n");
	fmpq_mpoly_matrix_fprint_pretty(stdout, x_train, ring_gen_names, ctx );
	printf("------------\n");

	fmpq_mpoly_matrix_mul(y_hat0, A, x_train, ctx );
	printf("------------\n");
	fmpq_mpoly_matrix_fprint_pretty(stdout, y_hat0, ring_gen_names, ctx );
	fmpq_mpoly_matrix_mul(y_hat1, B, x_train, ctx );
	printf("------------\n");
	fmpq_mpoly_matrix_fprint_pretty(stdout, y_hat1, ring_gen_names, ctx );
	fmpq_mpoly_matrix_hadamard(y_hat, y_hat0, y_hat1, ctx );
	printf("------------\n");
	fmpq_mpoly_matrix_fprint_pretty(stdout, y_hat, ring_gen_names, ctx );
	fmpq_mpoly_matrix_sub(E_vec, y_train, y_hat, ctx );
	printf("------------\n");
	fmpq_mpoly_matrix_fprint_pretty(stdout, E_vec, ring_gen_names, ctx );
	fmpq_mpoly_matrix_squared_frobenius(E, E_vec, ctx );
	printf("------------\n");
	fmpq_mpoly_print_pretty(E, ring_gen_names, ctx );
	fmpq_mpoly_add(E, E, constraint, ctx );
	printf("\n------------\n");
	fmpq_mpoly_print_pretty(E, ring_gen_names, ctx );

	printf("\n-----------\nLagrangian gradient:\n");
	for(int i=0; i < 6; i++ )
	{
		fmpq_mpoly_derivative(L_grad[i], E, i, ctx );
		fmpq_mpoly_print_pretty(L_grad[i], ring_gen_names, ctx );
		putchar('\n');
	}

	msolve_re_solutions_t sols;

	//data_gens_ff_t *gens = gens2msolve(L_grad, 6, ctx );
	msolve_from_fmpq_mpolys(sols, L_grad, 6, ring_gen_names, ctx );

	//now printing solutions
	printf("number of solutions found: %d\n", sols->nb_real_roots );

	printf("printing from buffer:\n");;
	
	double *sols_buf = get_sols_buffer(sols);

	for(size_t i=0; i < 6*16; i+=6 )
	{
		for(size_t j=0; j < 6; j++ )
			printf("%f ,", sols_buf[i+j] );
		putchar('\n');
	}

}
