#include "include/fmpq_mpoly_matrix.h"

void fmpq_mpoly_matrix_init(fmpq_mpoly_matrix_t A, size_t rows, size_t cols, const fmpq_mpoly_ctx_t ctx )
{

	A->cols = cols;
	A->rows = rows;
	A->cfs = malloc(cols*rows*sizeof(fmpq_mpoly_t));

	for(int i=0; i < cols*rows; i++ )
		fmpq_mpoly_init(A->cfs[i], ctx );

}

void fmpq_mpoly_matrix_clear(fmpq_mpoly_matrix_t A, const fmpq_mpoly_ctx_t ctx )
{

	for(int i=0; i < A->cols*A->rows; i++ )
		fmpq_mpoly_clear(A->cfs[i], ctx );

	free(A->cfs);
	
}

//Warning: matrix mul cannot be done in place !!!

int fmpq_mpoly_matrix_mul(fmpq_mpoly_matrix_t A, fmpq_mpoly_matrix_t B, fmpq_mpoly_matrix_t C, const fmpq_mpoly_ctx_t ctx )
{

	//Einsum notation
	//B[ij],C[jk] -> A[ik]

	if( (B->cols != C->rows)||(A->rows != B->rows)||(A->cols != C->cols) )
		return -1;

	fmpq_mpoly_t temp;
	fmpq_mpoly_init(temp, ctx );

	for(int i=0; i < B->rows; i++ )
	{
		for(int k=0; k < C->cols; k++ )
		{
			fmpq_mpoly_zero(A->cfs[A->cols*i+k], ctx );
			for(int j=0; j < B->cols; j++ )
			{
				fmpq_mpoly_mul(temp, B->cfs[B->cols*i+j], C->cfs[C->cols*j+k], ctx );
				fmpq_mpoly_add(A->cfs[A->cols*i+k], A->cfs[A->cols*i+k], temp, ctx );
			}
		}
	}

	return 0;

}

int fmpq_mpoly_matrix_hadamard(fmpq_mpoly_matrix_t A, fmpq_mpoly_matrix_t B, fmpq_mpoly_matrix_t C, const fmpq_mpoly_ctx_t ctx )
{
	
	if( (A->cols != B->cols)||(B->cols != C->cols)||(A->rows != B->rows)||(B->cols != C->cols) )
		return -1;

	for(int i=0; i < A->cols*A->rows; i++ )
		fmpq_mpoly_mul(A->cfs[i], B->cfs[i], C->cfs[i], ctx );

	return 0;

}

int fmpq_mpoly_matrix_add(fmpq_mpoly_matrix_t A, fmpq_mpoly_matrix_t B, fmpq_mpoly_matrix_t C, const fmpq_mpoly_ctx_t ctx )
{
	
	if( (A->cols != B->cols)||(B->cols != C->cols)||(A->rows != B->rows)||(B->cols != C->cols) )
		return -1;

	for(int i=0; i < A->cols*A->rows; i++ )
		fmpq_mpoly_add(A->cfs[i], B->cfs[i], C->cfs[i], ctx );

	return 0;

}

int fmpq_mpoly_matrix_sub(fmpq_mpoly_matrix_t A, fmpq_mpoly_matrix_t B, fmpq_mpoly_matrix_t C, const fmpq_mpoly_ctx_t ctx )
{
	
	if( (A->cols != B->cols)||(B->cols != C->cols)||(A->rows != B->rows)||(B->cols != C->cols) )
		return -1;

	for(int i=0; i < A->cols*A->rows; i++ )
		fmpq_mpoly_sub(A->cfs[i], B->cfs[i], C->cfs[i], ctx );

	return 0;

}


void fmpq_mpoly_matrix_squared_frobenius(fmpq_mpoly_t a, const fmpq_mpoly_matrix_t B, const fmpq_mpoly_ctx_t ctx )
{

	fmpq_mpoly_zero(a, ctx );
	fmpq_mpoly_t temp;
	fmpq_mpoly_init(temp, ctx );

	for(int i=0; i < B->cols*B->rows; i++ )
	{
		fmpq_mpoly_mul(temp, B->cfs[i], B->cfs[i], ctx );
		fmpq_mpoly_add(a, a, temp, ctx );
	}

}

int fmpq_mpoly_matrix_fprint_pretty(FILE *file, const fmpq_mpoly_matrix_t A, const char **x, const fmpq_mpoly_ctx_t ctx )
{
	for(int i=0; i < A->rows; i++ )
	{
		for(int j=0; j < A->cols; j++ )
		{
			if
			(
				( fmpq_mpoly_fprint_pretty(file, A->cfs[A->cols*i+j], x, ctx ) == -1 ) ||
				( putchar('\t') == EOF )
			) return -1;
		}
		if( putchar('\n') == EOF )
			return -1;
	}
	return 0;
}

fmpq_mpoly_struct *fmpq_mpoly_coeff_at(const fmpq_mpoly_matrix_t A, size_t i, size_t j )
{ return A->cfs[A->cols*i+j]; }

int fmpq_mpoly_matrix_coeff_set_str_pretty(const fmpq_mpoly_matrix_t A, const char *str, const char **x, size_t i, size_t j, const fmpq_mpoly_ctx_t ctx )
{ return fmpq_mpoly_set_str_pretty(A->cfs[A->cols*i+j], str, x, ctx ); }
