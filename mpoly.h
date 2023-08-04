#pragma once

#ifndef __MPOLY_SPARSE_H__
#define __MPOLY_SPARSE_H__

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <gmp.h>

//polynomial in QQ[..] block sparse rappresentation of exponents
typedef struct
{
	int32_t len;
	int32_t explen;
	int32_t *exps;
	mpq_t *coeffs;
} mpq_mpoly_t;

//QQ[..]
typedef struct
{
	int32_t numvars;
	char **varnames;
} ring_ctx_t;

void fprint_mpoly(FILE *fp, ring_ctx_t *ctx, mpq_mpoly_t *mpoly );
void qsort_mpq_mpoly(mpq_mpoly_t *mpoly );

#endif
