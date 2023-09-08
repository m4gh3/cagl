#include "../include/fmpq_mpoly_matrix.h"

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

	return 0;

}
