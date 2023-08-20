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
#include "../include/gens2msolve.h"


static void mpz_upoly_init(mpz_upoly_t poly, long alloc)
{
	mpz_t *tmp = NULL;
	if(alloc)
	{
		tmp = (mpz_t *)(malloc(alloc * sizeof(mpz_t)));
		if(tmp==NULL)
		{
			fprintf(stderr, "Unable to allocate in mpz_upoly_init\n");
			exit(1);
		}
		for(long i = 0; i < alloc; i++)
		{
			mpz_init(tmp[i]);
			mpz_set_ui(tmp[i], 0);
		}
	}
	poly->coeffs = tmp;
	poly->alloc = alloc;
	poly->length = -1;
}

void mpz_upoly_clear(mpz_upoly_t pol);

static inline void mpz_param_init(mpz_param_t param)
{
	param->nvars = 0;
	param->nsols = 0;
	mpz_upoly_init(param->elim, 0);
	mpz_upoly_init(param->denom, 0);
	param->coords=NULL;
	param->cfs=NULL;
}

static inline void mpz_param_clear(mpz_param_t param)
{
	mpz_upoly_clear(param->elim);
	mpz_upoly_clear(param->denom);
	if(param->coords != NULL)
	{
		for(long i = 0; i < param->nvars - 1; i++)
		{
			mpz_upoly_clear(param->coords[i]);
			mpz_clear(param->cfs[i]);
	  	}
	}
	free(param->coords);
	free(param->cfs);
	param->nvars  = 0;
	param->nsols  = 0;
}


void real_point_clear(real_point_t pt);


#define DEBUGGB 0
#define DEBUGBUILDMATRIX 0
#define IO_DEBUG 0

int main(int argc, char **argv)
{

	/* timinigs */
	double st0 = cputime();
	double rt0 = realtime();
	
	/**
	  We get values from the command line.
	 **/
	int32_t la_option             = 2; // by default
	int32_t use_signatures        = 0;
	int32_t nr_threads            = 1;
	int32_t info_level            = 0;
	int32_t initial_hts           = 17;
	int32_t max_pairs             = 0;
	int32_t elim_block_len        = 0;
	int32_t update_ht             = 0;
	int32_t generate_pbm          = 0;
	int32_t reduce_gb             = 1;
	int32_t print_gb              = 0;
	int32_t genericity_handling   = 2;
	int32_t saturate              = 0;
	int32_t colon                 = 0;
	int32_t normal_form           = 0;
	int32_t normal_form_matrix    = 0;
	int32_t is_gb                 = 0;
	int32_t get_param             = 0;
	int32_t precision             = 128;
	int32_t refine                = 0; /* not used at the moment */
	int32_t isolate               = 0; /* not used at the moment */
	
	files_gb *files = malloc(sizeof(files_gb));
	files->in_file = NULL;
	files->bin_file = NULL;
	files->out_file = NULL;
	files->bin_out_file = NULL;


	fmpq_mpoly_ctx_t ctx;
	
	fmpq_mpoly_ctx_init(ctx, 2, ORD_DEGREVLEX );
	
	const char *ring_gen_names[] = {"x", "y" };
	
	fmpq_mpoly_t P[2];
	
	fmpq_mpoly_init(P[0], ctx );
	fmpq_mpoly_init(P[1], ctx );
	
	fmpq_mpoly_set_str_pretty(P[0], "x^2+y^2-1", ring_gen_names, ctx );
	fmpq_mpoly_set_str_pretty(P[1], "x-y", ring_gen_names, ctx );
	
	data_gens_ff_t *gens = gens2msolve(P, 2, ctx );
	    
	/* data structures for parametrization */
	param_t *param  = NULL;
	mpz_param_t mpz_param;
	mpz_param_init(mpz_param);
	
	long nb_real_roots      = 0;
	interval *real_roots    = NULL;
	real_point_t *real_pts  = NULL;
	
	/* main msolve functionality */
	int ret = core_msolve(la_option, use_signatures, nr_threads, info_level,
	                      initial_hts, max_pairs, elim_block_len, update_ht,
	                      generate_pbm, reduce_gb, print_gb, get_param,
	                      genericity_handling, saturate, colon, normal_form,
	                      normal_form_matrix, is_gb, precision, 
	                      files, gens,
	        &param, &mpz_param, &nb_real_roots, &real_roots, &real_pts);
	
	/* free parametrization */
	//free(param);
	mpz_param_clear(mpz_param);
	
	
	if (nb_real_roots > 0) {
	    for(long i = 0; i < nb_real_roots; i++){
	      real_point_clear(real_pts[i]);
	      mpz_clear(real_roots[i].numer);
	    }
	    free(real_pts);
	}
	free(real_roots);
	
	/* timings */
	if (info_level > 0) {
	    double st1 = cputime();
	    double rt1 = realtime();
	    fprintf(stderr, "-------------------------------------------------\
	-------------------------------\n");
	    fprintf(stderr, "msolve overall time  %13.2f sec (elapsed) / %5.2f sec (cpu)\n",
	            rt1-rt0, st1-st0);
	    fprintf(stderr, "-------------------------------------------------\
	-------------------------------\n");
	}

	free_data_gens(gens);
	free(files);

	return ret;

}
