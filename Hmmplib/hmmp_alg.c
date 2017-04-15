/*  This file is part of Hmmplib.
*	Description: Hmmplib is a powerful multy-core library solution of Hidden Markov model in C.
*   Copyright (C) 2017  Iskren Tarkalanov
*
*   Hmmplib is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   any later version.
*
*   Hmmplib is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with Hmmplib.  If not, see <http://www.gnu.org/licenses/>.
*
*	Contact: isktark@yahoo.com Iskren Tarkalanov
*/
#include "hmmp_alg.h"
#include "hmmp_datatypes.h"
#include "hmmp_memop.h"
#include "hmmp_dataproc.h"
#include <math.h>
#include <omp.h>
int hmmp_forward_alg (	hmmp_Model model, hmmp_Sequence seq, dbl_matrix *o_alfa,
						dbl_array *o_alfa_scale )
{
	int i, j, t;
	double part_sum;
	dbl_matrix *transit_indx_ptr, *alfa_indx_ptr;

	if(!o_alfa_scale || !o_alfa )
		return E_PARAMETER;
	//initialization
	for ( i = 0 ; i < model.num_states ; ++i ){
		o_alfa[i] = model.initial[i] * model.emission[i*model.num_symbols+seq.sequence[0]]; //variant 2
	}
	o_alfa_scale[0] = hmmp_normalize_arr(o_alfa,model.num_states);

	alfa_indx_ptr = o_alfa;
	//induction
	for ( t = 1 ; t < seq.length ; ++t ){
		transit_indx_ptr = model.transition;
		for ( i = 0 ; i < model.num_states ; ++i ){
			part_sum = 0.0;
			for ( j = 0; j < model.num_states ; ++j ){
				part_sum += alfa_indx_ptr[j] * transit_indx_ptr[j*model.num_states];
			}
			alfa_indx_ptr[model.num_states+i]=part_sum * model.emission[i*model.num_symbols+seq.sequence[t]];
			++transit_indx_ptr;
		}
		alfa_indx_ptr += model.num_states;	// indexing value to the next time step
		o_alfa_scale[t] = hmmp_normalize_arr(alfa_indx_ptr,model.num_states);
	}
	return 0;
}

int hmmp_backward_alg ( hmmp_Model model, hmmp_Sequence seq, dbl_matrix *o_beta,
					    dbl_array *o_beta_scale )
{
	int i, j, t; //  i = 0,1,...,model.num_states-1,model.num_states(num_states) ; t = 0,1,...,seq.length-1( seq.length - sequence seq.sequence length )
	double part_sum;
	dbl_array *beta_helper, *indx_helper;

	if( !o_beta || !o_beta_scale )
		return E_PARAMETER;
	beta_helper = hmmp_create_dbl_array(model.num_states);
	if ( !beta_helper )
		return E_ALLOCATION;

	//initialization
	//scaling factors at last time step
	o_beta_scale[seq.length-1]=1;

	//using address arithmetic to do less operations
	indx_helper = o_beta +  ( seq.length - 1 ) * model.num_states ;
	for ( i = 0 ; i < model.num_states ; ++i ){
		indx_helper[i] = 1;
	}
	//induction
	for ( t = seq.length - 2 ; t >= 0 ; --t ){
		for ( i = 0 ; i < model.num_states ; ++i )
			//indx_helper[i] represents beta values at the next time step.
			beta_helper[i] = indx_helper[i] * 
							 model.emission[i*model.num_symbols+seq.sequence[t+1]];
		for ( i = 0 ; i < model.num_states ; ++i ){
			part_sum = 0.0;
			for ( j = 0; j < model.num_states ; ++j ){
				part_sum += model.transition[i*model.num_states + j] * beta_helper[j]; 
				// transision prob from each previus to current state indexed j->i
			}
			o_beta[t*model.num_states + i] = part_sum;
		}
		o_beta_scale[t]=hmmp_normalize_arr(o_beta+t*model.num_states,model.num_states);
		indx_helper-=model.num_states;
	}
	hmmp_delete_dbl_array(beta_helper);
	return E_SUCCESS;
}
//int hmm_viterbi_alg(hmmp_Model model, int_array *obs, int obs_len, int *backtrack, 
//					double *mu, int* bestpath, double *prob_bp )
int hmmp_backward_rescale( dbl_matrix *beta, int num_states, int seq_len,
						   dbl_array *scale_alfa, dbl_array *scale_beta )
{
	double rescale = 1.0;
	int t, i;
	if(!beta || !scale_alfa || !scale_beta )
		return E_PARAMETER;
	for ( t = seq_len-1 ; t >= 0 ; --t){
		rescale *= scale_alfa[t]/scale_beta[t];
		for ( i = 0 ; i < num_states ; ++i ){
			beta[t*num_states+i] *= rescale;
		}
	}
	return E_SUCCESS;
}

int hmmp_viterbi_alg( hmmp_Model log_model,	hmmp_Sequence seq, int_array *backtrack, 
					 dbl_matrix *mu,	hmmp_Sequence *o_state_seq, double *o_logP )
{
	int i, j, t;
	int backtrack_i = 0;
	double mu_max, swap_val;
	dbl_matrix *mu_old, *swap_ptr;

	if ( !backtrack || !mu || !o_state_seq || !o_logP )
		return E_PARAMETER;
	// Initialization
	for ( i = 0 ; i < log_model.num_states ; ++i )
		mu[i] = log_model.initial[i] + log_model.emission[i*log_model.num_symbols+seq.sequence[0]];

	mu_old = mu;
	mu += log_model.num_states;
	//Induction
	for ( t = 1 ; t < seq.length ; ++t ){


		for ( i = 0 ; i < log_model.num_states ; ++i ){
			mu_max = -HMMP_DBL_MAX;
			swap_val = 0;
			for ( j = 0; j < log_model.num_states ; ++j ){
				swap_val = mu_old[j] + log_model.transition[j*log_model.num_states+i];;
				if ( mu_max < swap_val ){
					mu_max = swap_val;
					backtrack_i = j;
				}
			}
			mu[i] = mu_max + log_model.emission[i*log_model.num_symbols+seq.sequence[t]];
			backtrack[t*log_model.num_states + i] = backtrack_i;
		}
		swap_ptr = mu;
		mu = mu_old;
		mu_old = swap_ptr;

	}
	//Termination
	mu_max = -HMMP_DBL_MAX;
	for ( i = 0 ; i < log_model.num_states ; ++i ){
		if ( mu_max < mu_old[i] ){
			mu_max = mu_old[i];
			backtrack_i = i;
		}
	}
	//Backtracking best path
	o_state_seq->sequence[seq.length-1] = backtrack_i;
	for ( t = seq.length-1 ; t > 0 ; --t ){
		backtrack_i = backtrack[t*log_model.num_states + backtrack_i];
		o_state_seq->sequence[t-1] = backtrack_i;	
	}
	o_state_seq->sequence[0] = backtrack_i;
	// Save OUT: param probability for best path
	*o_logP = mu_max;
	return E_SUCCESS;
}

int hmmp_bwa_gamma_alg ( dbl_matrix *o_gamma, dbl_matrix *alfa, dbl_matrix *beta, 
						dbl_array *alfa_scale, int num_states, int seq_length )
{
	int i, t;
	dbl_matrix *pa, *pb;
	if ( !o_gamma || !alfa || !beta || !alfa_scale )
		return E_PARAMETER;
#pragma omp for private(t,i,pa,pb) schedule(static)
	for ( t = 0 ; t < seq_length ; ++t ){
		pa = alfa + t*num_states;
		pb = beta + t*num_states;
		for ( i = 0 ; i < num_states ; ++i)
			o_gamma[i*seq_length+t] = pa[i] * pb[i] / alfa_scale[t];
	}
	return E_SUCCESS;
}
int hmmp_bwa_xi_alg ( dbl_matrix *o_xi, dbl_matrix *alfa, dbl_matrix *beta, 
					 hmmp_Model model, hmmp_Sequence seq )
{
	int i, j, t, NtT, N;
	dbl_matrix *pa, *pb;
	double emit;

	if ( !o_xi || !alfa || !beta )
		return E_PARAMETER;
	--seq.length;
	N = model.num_states;
	NtT = N*seq.length;
#pragma omp for private(t,i,j,pa,pb) schedule(static)
	for ( t = 0 ; t < seq.length ; ++t ){
		pa = alfa + t*N;
		pb = beta + (t+1)*N;
		for ( j = 0 ; j < N ; ++j ){
			emit = model.emission[j*model.num_symbols+seq.sequence[t+1]];
			for ( i = 0 ; i < N ; ++i )
				o_xi[i*NtT+j*seq.length+t]=pa[i]*model.transition[i*N+j]*emit*pb[j];
		}
	}
	return E_SUCCESS;
}

int hmmp_bwa_reest_alg ( hmmp_Model model,		hmmp_Sequence seq,
						 dbl_matrix *xi,		dbl_matrix *gamma,
						 dbl_matrix *o_a_num,	dbl_matrix *o_b_num,
						 dbl_array *o_a_denom,	dbl_array *o_b_denom	)
{
	int i, j ,t;
	double part_sum;
	dbl_matrix *xiindx, *gammaindx;
	--seq.length;

	if ( !xi || !gamma || !o_a_num || !o_b_num || !o_a_denom || !o_b_denom )
		return E_PARAMETER;
#pragma omp for private(i,j,t,xiindx) schedule(static)
	for ( i = 0 ; i < model.num_states; ++i ){
		for ( j = 0 ; j < model.num_states ; ++j ){
			// j loop for a ( transition ) numerators ksi(t)[i]->[j];
			xiindx = xi + i*model.num_states*seq.length + j*seq.length;
			for ( t = 0 ; t < seq.length ; ++t ){
				o_a_num[i*model.num_states + j] += xiindx[t];
			}
		}
	}
#pragma omp for private(i,j,t,gammaindx, part_sum) schedule(static)
	for ( i = 0 ; i < model.num_states; ++i ){
		gammaindx = gamma + i*(seq.length+1);
		part_sum = 0.0;
		for ( t = 0 ; t < seq.length ; ++t ){
			// denominator summation of gamma(t)[i] for each state over 't'
			part_sum += gammaindx[t];
			// adding gamma(t)[i] to b_ik when k depends on O(t)
			
			o_b_num[i*model.num_symbols + seq.sequence[t]] += gammaindx[t];
		}
		// denominator for: a_ij over t=1:T-1 | b_ik summation over t=1:T
		o_a_denom[i] += part_sum;
		o_b_denom[i] += part_sum + gammaindx[t];
		o_b_num[i*model.num_symbols + seq.sequence[t]] += gammaindx[t];

	}
	return E_SUCCESS;
}