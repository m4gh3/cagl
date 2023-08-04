#include "mpoly.h"

#define FPRINT_MONOMIAL(x,s) \
	{ \
		gmp_fprintf(fp, s "%Qd", mpoly->coeffs[(x)] ); \
		for(int j=0; (j < explen)&&(exps[(x)*(explen<<1)+(j<<1)+1] != 0); j++ ) \
			printf("*%s^%d", ctx->varnames[exps[(x)*(explen<<1)+(j<<1)]], exps[(x)*(explen<<1)+(j<<1)+1] ); \
	}

void fprint_mpoly(FILE *fp, ring_ctx_t *ctx, mpq_mpoly_t *mpoly )
{

	int32_t explen = mpoly->explen;
	int32_t *exps  = mpoly->exps;

	FPRINT_MONOMIAL(0, "" )

	for(int i=1; i < mpoly->len; i++ )
		FPRINT_MONOMIAL(i, "+" )

}

static void _part_mpq_mpoly(mpq_mpoly_t *mpoly, size_t lo, size_t hi )
{
	int32_t p_exp[(mpoly->explen)<<1];
	memcpy(p_exp, mpoly );
}

static void _qsort_mpq_mpoly(mpq_mpoly_t *mpoly, size_t lo, size_t hi )
{

  	// Ensure indices are in correct order
  	if( (lo >= hi) || (lo < 0) )
		return;
    
	// Partition array and get the pivot index
	int p = partition(mpoly, lo, hi );
      
	// Sort the two partitions
	_qsort_mpq_mpoly(mpoly, lo, p - 1 ); // Left side of pivot
	_qsort_mpq_mpoly(mpoly, p + 1, hi ); // Right side of pivot

}


void qsort_mpq_mpoly(mpq_mpoly_t *mpoly )
{
	
}


int main()
{

	ring_ctx_t ring =
	{
		.numvars = 2,
		.varnames = (char *[2]) {"x", "y" }
	};

	mpq_t coeffs[2];
	mpq_inits(coeffs[0], coeffs[1], NULL );

	mpq_set_ui(coeffs[0], 1, 2 );
	mpq_set_ui(coeffs[1], 1, 2 );

	mpq_mpoly_t p =
	{
		.len = 2,
		.explen = 2,
		.exps = (int32_t [8] ){0,2,0,0, 1,2,0,0 },
		.coeffs = coeffs
	};
	
	fprint_mpoly(stdout, &ring, &p );

	return 0;

}
