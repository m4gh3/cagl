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

#include "msolve/src/fglm/data_fglm.h"
#include "msolve/src/usolve/data_usolve.h"
#include "msolve/src/msolve/msolve.h"


static inline data_gens_ff_t *allocate_data_gens()
{

  data_gens_ff_t *gens = (data_gens_ff_t *)(malloc(sizeof(data_gens_ff_t)));
  gens->lens  = NULL;
  gens->exps  = NULL;
  gens->cfs   = NULL;
  gens->mpz_cfs = NULL;

  gens->elim = 0;
  return gens;

}

void gens2msolve(data_gens_ff_t *gens, fmpq_mpoly_t *polys, size_t n, fmpq_mpoly_ctx_t ctx )
{

	//slong nr_terms = 0;
	
	gens->nvars = fmpq_mpoly_ctx_nvars(ctx);
	gens->ngens = n;
	gens->nterms = 0;
	gens->field_char = 0;


	gens->lens = (int32_t *)malloc((unsigned long)(gens->ngens) * sizeof(int32_t));

	ulong exps[gens->nvars];

	for(size_t i=0; i < n; i++ )
	{
		gens->lens[i] = fmpq_mpoly_length(polys[i], ctx );
		gens->nterms +=  gens->lens[i];
	}

	gens->exps = (int32_t *) malloc(gens->nterms * gens->nvars * sizeof(int32_t));
	gens->mpz_cfs = (mpz_t **)(malloc(sizeof(mpz_t *) * 2 * gens->nterms));

	for(int i=0; i < 2*gens->nterms; i++ )
	{
		gens->mpz_cfs[i] = (mpz_t *)malloc(sizeof(mpz_t));
		mpz_init(*(gens->mpz_cfs[i]));
	}

	for(size_t i=0,p=0; i < n; i++ )
	{
		fmpq_t c;
		fmpq_init(c);
		for(size_t j=0; j < gens->lens[i]; j++, p++ )
		{

			fmpq_mpoly_get_term_coeff_fmpq(c, polys[i], j, ctx );
			fmpq_get_mpz_frac(gens->mpz_cfs[2*p][0], gens->mpz_cfs[2*p+1][0], c );
			fmpq_mpoly_get_term_exp_ui(exps/*gens->exps + gens->nvars*p*/, polys[i], j, ctx );

			for(size_t k=0; k < gens->nvars; k++ )
				(gens->exps + gens->nvars*p)[k] = (int32_t) exps[k];

		}
	}

}
