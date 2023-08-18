#include <stdio.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpq_mpoly.h>

int main()
{

    // Initialize the context: !!! void fmpq_mpoly_ctx_init(fmpq_mpoly_ctx_t ctx, slong nvars, const ordering_t ord) !!!
    fmpq_mpoly_ctx_t ctx;

    fmpq_mpoly_ctx_init(ctx, 2, ORD_DEGREVLEX );

    // Create array of string with the generators names
    
    const char *ring_gen_names[2] = {"x", "y" };

    // Initialize the variables
    fmpq_mpoly_t A, B, C;

    // Initialize the polynomials
    fmpq_mpoly_init(A, ctx );
    fmpq_mpoly_init(B, ctx );
    fmpq_mpoly_init(C, ctx );

    // Set coefficients for the first polynomial: A(x, y) = x^2 + 2xy + 1
    fmpq_t a1, a2, a3;
    fmpq_init(a1);
    fmpq_init(a2);
    fmpq_init(a3);

    fmpq_set_si(a1, 1, 1);
    fmpq_set_si(a2, 2, 1);
    fmpq_set_si(a3, 1, 1);

    // Create the array of ulongs for the exponents of the monomials of the first polynomial

    ulong e1[2] = {2,0};
    ulong e2[2] = {1,1};
    ulong e3[2] = {0,0};

    //void fmpq_mpoly_set_coeff_fmpq_ui(fmpq_mpoly_t A, const fmpq_t c, ulong const *exp, fmpq_mpoly_ctx_t ctx)

    fmpq_mpoly_set_coeff_fmpq_ui(A, a1, e1, ctx );
    fmpq_mpoly_set_coeff_fmpq_ui(A, a2, e2, ctx );
    fmpq_mpoly_set_coeff_fmpq_ui(A, a3, e3, ctx );

    // Set coefficients for the second polynomial: B(x, y) = 3xy^2 + 2
    fmpq_t b1, b2;
    fmpq_init(b1);
    fmpq_init(b2);

    fmpq_set_si(b1, 3, 1);
    fmpq_set_si(b2, 2, 1);

    // Create the array of ulongs for the exponents of the monomials of the second polynomial
    
    ulong f1[2] = {1,2};
    ulong f2[2] = {0,0};

    fmpq_mpoly_set_coeff_fmpq_ui(B, b1, f1, ctx );
    fmpq_mpoly_set_coeff_fmpq_ui(B, b2, f2, ctx );

    // Multiply the polynomials: C(x, y) = A(x, y) * B(x, y)
    fmpq_mpoly_mul(C, A, B, ctx );

    // Print the result
    fmpq_mpoly_print_pretty(C, ring_gen_names, ctx );
    putchar('\n'); 

    // get lenght of C
    printf("The polynomial C has length %lu\n", fmpq_mpoly_length(C, ctx ) );

    ulong C_len = fmpq_mpoly_length(C, ctx );
    
    //extract the coefficients:
   
    printf("extracting the coeffs:\n\t");

    mpq_t c; fmpq_t c_; mpq_init(c); fmpq_init(c_);
    for(ulong i=0; i < C_len; i++ )
    {
        fmpq_mpoly_get_term_coeff_fmpq(c_, C, i, ctx );
	//fmpq_print(c_); putchar('\n');
	fmpq_get_mpq(c, c_ );
	gmp_printf("%Qd,", c );
    }
    mpq_clear(c); fmpq_clear(c_);

    putchar('\n');

    //extract the exponents
   
    printf("extracting the exponents: {\n");

    ulong exps[2];
    for(ulong i=0; i < C_len; i++ )
    {
        //void fmpq_mpoly_get_term_exp_ui(ulong *exps, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)
	fmpq_mpoly_get_term_exp_ui(exps, C, i, ctx );
	printf("\t%lu : {", i );
	for(ulong j=0; j < 2; j++ )
	{
            printf("%lu, ", exps[j] );
	};
	printf("},\n");
    }
    printf("}\n");

    // Clear memory
    fmpq_clear(a1);
    fmpq_clear(a2);
    fmpq_clear(a3);
    fmpq_clear(b1);
    fmpq_clear(b2);
    fmpq_mpoly_clear(A, ctx );
    fmpq_mpoly_clear(B, ctx );
    fmpq_mpoly_clear(C, ctx );

    return 0;

}
