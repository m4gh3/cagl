/* This file is part of msolve.
 *
 * msolve is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * msolve is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with msolve.  If not, see <https://www.gnu.org/licenses/>
 *
 * Authors:
 * Jérémy Berthomieu
 * Christian Eder
 * Mohab Safey El Din */

//#include "libmsolve.c"
#include "msolve/src/fglm/data_fglm.h"
#include "msolve/src/usolve/data_usolve.h"
#include "msolve/src/msolve/msolve.h"

static void mpz_upoly_init(mpz_upoly_t poly, long alloc){
  mpz_t *tmp = NULL;
  if(alloc){
    tmp = (mpz_t *)(malloc(alloc * sizeof(mpz_t)));
    if(tmp==NULL){
      fprintf(stderr, "Unable to allocate in mpz_upoly_init\n");
      exit(1);
    }
    for(long i = 0; i < alloc; i++){
      mpz_init(tmp[i]);
      mpz_set_ui(tmp[i], 0);
    }
  }
  poly->coeffs = tmp;
  poly->alloc = alloc;
  poly->length = -1;
}

static void mpz_upoly_init2(mpz_upoly_t poly, long alloc, long nbits){
  mpz_t *tmp = NULL;
  if(alloc){
    tmp = (mpz_t *)(malloc(alloc * sizeof(mpz_t)));
    //    tmp = (mpz_t *)(calloc(alloc, sizeof(mpz_t)));
    if(tmp==NULL){
      fprintf(stderr, "Unable to allocate in mpz_upoly_init\n");
      exit(1);
    }
    for(long i = 0; i < alloc; i++){
      mpz_init2(tmp[i], nbits);
      mpz_set_ui(tmp[i], 0);
    }
  }
  poly->coeffs = tmp;
  poly->alloc = alloc;
  poly->length = -1;
}

void mpz_upoly_clear(mpz_upoly_t pol); /*{

  for(long i = 0; i < pol->alloc; i++){
    mpz_clear(pol->coeffs[i]);
  }
  free(pol->coeffs);
}*/

static inline void mpz_param_init(mpz_param_t param){
  param->nvars = 0;
  param->nsols = 0;
  mpz_upoly_init(param->elim, 0);
  mpz_upoly_init(param->denom, 0);
  param->coords=NULL;
  param->cfs=NULL;
}

static inline void mpz_param_clear(mpz_param_t param){
  mpz_upoly_clear(param->elim);
  mpz_upoly_clear(param->denom);
  if(param->coords != NULL){
    for(long i = 0; i < param->nvars - 1; i++){
      mpz_upoly_clear(param->coords[i]);
      mpz_clear(param->cfs[i]);
    }
  }
  free(param->coords);
  free(param->cfs);
  param->nvars  = 0;
  param->nsols  = 0;
}


void real_point_clear(real_point_t pt);/*{
  for(long i = 0; i < pt->nvars; i++){
    mpz_clear(pt->coords[i]->val_up);
    mpz_clear(pt->coords[i]->val_do);
  }
  free(pt->coords);
}*/

static inline data_gens_ff_t *allocate_data_gens(){
  data_gens_ff_t *gens = (data_gens_ff_t *)(malloc(sizeof(data_gens_ff_t)));
  gens->lens  = NULL;
  gens->exps  = NULL;
  gens->cfs   = NULL;
  gens->mpz_cfs = NULL;

  gens->elim = 0;
  return gens;
}

static inline void free_data_gens(data_gens_ff_t *gens){
  /*for(long i = 0; i < gens->nvars; i++){
    free(gens->vnames[i]);
  }
  free(gens->vnames);*/
  if (gens->field_char == 0) {
      for(long i = 0; i < 2*gens->nterms; i++){
          mpz_clear(*(gens->mpz_cfs[i]));
          free(gens->mpz_cfs[i]);
      }
  }
  free(gens->mpz_cfs);
  free(gens->lens);
  free(gens->cfs);
  free(gens->exps);
  free(gens->random_linear_form);
  free(gens);
}


#define DEBUGGB 0
#define DEBUGBUILDMATRIX 0
#define IO_DEBUG 0

int main(int argc, char **argv){

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
    
    FILE *fh  = fopen(files->in_file, "r");
    FILE *bfh  = fopen(files->bin_file, "r");

    /*if (fh == NULL && bfh == NULL) {
      fprintf(stderr, "Input file not found.\n");
      exit(1);
    }*/
    if(fh!=NULL){
      fclose(fh);
    }
    if(bfh != NULL){
      fclose(bfh);
    }
    fh =  NULL;
    bfh =  NULL;

    /* clear out_file if given */
    /*if(files->out_file != NULL){
      FILE *ofile = fopen(files->out_file, "w");
      if(ofile == NULL){
        fprintf(stderr, "Cannot open output file\n");
        exit(1);
      }
      fclose(ofile);
    }*/
    /**
       We get from files the requested data. 
    **/
    //  int32_t mon_order   = 0;
    int32_t nr_vars     = 2;
    int32_t nr_gens     = 2;
    int32_t nr_terms    = 5;
    int32_t field_char  = 0;

    data_gens_ff_t *gens = allocate_data_gens();
 
    /* get_data_from_file(files->in_file, &nr_vars, &field_char, &nr_gens, gens); */
  
    gens->nvars = nr_vars;
    gens->ngens = nr_gens;
    gens->nterms = nr_terms;
    gens->field_char = 0;
    gens->lens = (int32_t *) malloc((unsigned long)nr_gens * sizeof(int32_t));
    gens->exps = (int32_t *) malloc(nr_terms * nr_vars * sizeof(int32_t));
    gens->mpz_cfs = (mpz_t **)(malloc(sizeof(mpz_t *) * 2 * nr_terms));

    gens->lens[0] = 3; gens->lens[1] = 2;

    gens->exps[0] = 2; gens->exps[1] = 0;
    gens->exps[2] = 0; gens->exps[3] = 2;
    gens->exps[4] = 0; gens->exps[5] = 0;
    gens->exps[6] = 1; gens->exps[7] = 0;
    gens->exps[8] = 0; gens->exps[9] = 1;

    gens->vnames = (char *[2]){"x","y"};

    for(int i=0; i < 2*nr_terms; i++ )
    {
    	gens->mpz_cfs[i] = (mpz_t *)malloc(sizeof(mpz_t));
	mpz_init(*(gens->mpz_cfs[i]));
    }

    signed long int coeffs[2*5/*nr_terms*/] =
    {
	    1,1, 1,1, -1,1, 1,1, -1,1
    };

    for(int i=0; i < 2*nr_terms; i++ )
	mpz_set_si(*(gens->mpz_cfs[i]), coeffs[i] );

    /*gens->rand_linear           = 0;
    gens->random_linear_form = malloc(sizeof(int32_t)*(nr_vars));
    gens->elim = elim_block_len;*/
    
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
-----------------------------------\n");
        fprintf(stderr, "msolve overall time  %13.2f sec (elapsed) / %5.2f sec (cpu)\n",
                rt1-rt0, st1-st0);
        fprintf(stderr, "-------------------------------------------------\
-----------------------------------\n");
    }
    free_data_gens(gens);
    /* for(long i = 0; i < gens->nvars; i++){
        free(gens->vnames[i]);
    }
    free(gens->vnames);
    free(gens->lens);
    free(gens->cfs);
    free(gens->exps);
    free(gens->random_linear_form); */
    free(files);
    return ret;
}
