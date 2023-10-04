/*
int main()
{
	fmpq_mpoly_ctx_t ctx;

	fmpq_mpoly_ctx_init(ctx, 2, ORD_DEGREVLEX );

	const char *ring_gen_names[] = {"x", "y" };

	fmpq_mpoly_matrix_t A, B;

	fmpq_mpoly_matrix_init(A, 2, 2, ctx );
	fmpq_mpoly_matrix_init(B, 2, 2, ctx );

	fmpq_mpoly_matrix_coeff_set_str_pretty(B, "x",   ring_gen_names, 0, 0, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(B, "y",   ring_gen_names, 0, 1, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(B, "x^2", ring_gen_names, 1, 0, ctx );
	fmpq_mpoly_matrix_coeff_set_str_pretty(B, "y^2", ring_gen_names, 1, 1, ctx );

	printf("A=\n");

	fmpq_mpoly_matrix_fprint_pretty(stdout, B, ring_gen_names, ctx );

	fmpq_mpoly_matrix_mul(A, B, B, ctx );

	printf("----------\nA^2=\n");

	fmpq_mpoly_matrix_fprint_pretty(stdout, A, ring_gen_names, ctx );

	fmpq_mpoly_t frob;
	fmpq_mpoly_init(frob, ctx );

	printf("----------\n");

	fmpq_mpoly_matrix_squared_frobenius(frob, A, ctx );
	fmpq_mpoly_print_pretty(frob, ring_gen_names, ctx );

	return 0;

}*/
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
#include "../msolve/src/usolve/data_usolve.h"
#include "../msolve/src/msolve/msolve.h"
#include "../include/fmpq_mpoly_matrix.h"

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

	//data_gens_ff_t *gens = gens2msolve(L_grad, 6, ctx );
	msolve_from_fmpq_mpolys(L_grad, 6, ring_gen_names, ctx );
	
}
