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
#include "hmmp_general.h"
#include "hmmp_datatypes.h"
#include "hmmp_memop.h"
#include "hmmp_dataproc.h"
#include "hmmp_alg.h"
#include "hmmp_algwrap.h"
#include <malloc.h>
#include <omp.h>


int hmmp_decode(hmmp_Model model, hmmp_Sequence *observ_array, int num_obs, 
				hmmp_Sequence **o_state_array, dbl_array **o_logPS )
{
	hmmp_Model *log_model = 0;
	hmmp_Sequence *states_arr = 0;
	dbl_array *arr_prob_state_seq = 0;
	int_matrix *backtrack = 0; // T x N matrix
	dbl_matrix *mu = 0;	// 2 x N matrix
	int k , max_length = 0;
	int e_overflow;
	char fail_flag = 0;
	if ( !observ_array || !o_state_array || !o_logPS )
		return E_PARAMETER;
	for ( k = 0 ; k < num_obs ; ++k )
		if ( max_length < observ_array[k].length )
			max_length = observ_array[k].length;
	e_overflow = hmmp_memop_overflow(model.num_states, max_length);
	if (e_overflow == E_MEM_OVERFLOW_L1){
		*o_state_array = 0;
		*o_logPS = 0;
		return e_overflow;
	}

	states_arr = (hmmp_Sequence *) malloc ( sizeof(hmmp_Sequence)*num_obs);
	if ( ! states_arr ){ fail_flag = 1; goto SKIP_REST; }	
	for ( k = 0 ; k < num_obs ; ++k ){
		states_arr[k].sequence=hmmp_create_int_array(observ_array[k].length);
		if(!(states_arr[k].sequence)){
			--k;
			for (; k >= 0 ; --k )
				hmmp_delete_int_array( states_arr[k].sequence );
			free ( states_arr );
			states_arr = 0;
			fail_flag = 1;
			goto SKIP_REST;
		}
	}

	arr_prob_state_seq = hmmp_create_dbl_array ( num_obs );
	if ( ! arr_prob_state_seq ){ fail_flag = 1; goto SKIP_REST;	}

	log_model = hmmp_create_model_copy ( &model );
	if ( !log_model ){ fail_flag = 1; goto SKIP_REST; }

	hmmp_model_log_param ( log_model );

#pragma omp parallel firstprivate ( backtrack, mu ) default(shared) num_threads(HMMP_NUM_THREADS)
{
	#pragma omp critical
	{
	#pragma omp flush ( fail_flag )
		if (!fail_flag){
			backtrack = hmmp_create_int_matrix ( model.num_states * max_length );
			if ( !backtrack ) fail_flag = 1;
			else{
				mu = hmmp_create_dbl_matrix ( 2 * model.num_states );
				if (!mu ) fail_flag = 1;
			}
		}
	}//END OF CRITICAL SECTION
	#pragma omp barrier
	#pragma omp flush ( fail_flag )
	if(!fail_flag){
		#pragma omp for schedule(static)
		for ( k = 0 ; k < num_obs ; ++k )
			hmmp_viterbi_alg ( *log_model, observ_array[k], backtrack, 
								mu,states_arr+k,arr_prob_state_seq+k );
	}
	#pragma omp critical
	{
		if ( mu ) hmmp_delete_dbl_matrix ( mu );
		if ( backtrack ) hmmp_delete_int_array ( backtrack );
	}//end of critical section
}//end of paralell region
	*o_state_array = states_arr;
	*o_logPS = arr_prob_state_seq;
SKIP_REST:
	if ( log_model ) hmmp_delete_model ( log_model );

	if ( fail_flag && arr_prob_state_seq )
		hmmp_delete_dbl_array ( arr_prob_state_seq );
	if ( fail_flag && states_arr)
		hmmp_delete_arr_seq(states_arr,num_obs);
	if ( fail_flag ){
		*o_state_array = 0;
		*o_logPS = 0;
		return E_ALLOCATION;
	}
	return E_SUCCESS;
}

int hmmp_evaluate_models(hmmp_Model *arr_models, int num_models, 
						 hmmp_Sequence observ_seq, dbl_array **o_logP_arr )
{
	int k, max_num_states = 0;
	int e_overflow;
	dbl_array *scales = 0, *prob_arr = 0;
	dbl_matrix *alfa = 0;
	char fail_flag = 0;
	if ( !arr_models || !o_logP_arr )
		return E_PARAMETER;
	for ( k = 0 ; k < num_models ; ++k )
		if ( max_num_states < arr_models[k].num_states )
			max_num_states = arr_models[k].num_states;
	e_overflow = hmmp_memop_overflow(max_num_states, observ_seq.length);
	if (e_overflow == E_MEM_OVERFLOW_L1){
		*o_logP_arr = 0;
		return e_overflow;
	}
	prob_arr = hmmp_create_dbl_array ( num_models );
	if ( ! prob_arr ){
		return E_ALLOCATION;
	}
#pragma omp parallel private(k) firstprivate(alfa,scales) default(shared) num_threads(HMMP_NUM_THREADS)
{
	#pragma omp critical 
	{
		#pragma omp flush ( fail_flag )
		if(!fail_flag){
			alfa = hmmp_create_dbl_matrix(observ_seq.length*max_num_states);
			if( ! alfa )
				fail_flag = 1; 
			else{
				scales = hmmp_create_dbl_array ( observ_seq.length );
				if(!scales)
					fail_flag = 1; 
			}
		}
	}
	#pragma omp barrier
	#pragma omp flush ( fail_flag )
	if(!fail_flag){
		#pragma omp for schedule(static)
		for ( k = 0 ; k < num_models ; ++k )
		{
			hmmp_forward_alg(arr_models[k],observ_seq,alfa,scales);
			prob_arr[k] = hmmp_log_of_divisors(scales,observ_seq.length);
		}
	}
	#pragma omp critical
	{
		if ( alfa ) 
			hmmp_delete_dbl_matrix(alfa);
		if ( scales ) 
			hmmp_delete_dbl_array(scales);
	}

}
	if ( fail_flag ){
		*o_logP_arr = 0;
		hmmp_delete_dbl_array(prob_arr);
		return E_ALLOCATION;
	}
	*o_logP_arr = prob_arr;
	return E_SUCCESS;
}

int hmmp_evaluate_sequences(hmmp_Model model, hmmp_Sequence *observ_arr, int num_obs,
							double **o_logP_arr )
{
	int k, max_length = 0;
	int e_overflow;
	dbl_array *scales = 0, *prob_arr = 0;
	dbl_matrix *alfa = 0;
	char fail_flag = 0;
	if ( ! observ_arr || !o_logP_arr )
		return E_PARAMETER;
	for ( k = 0 ; k < num_obs ; ++k )
		if ( max_length < observ_arr[k].length )
			max_length = observ_arr[k].length;
	if (e_overflow = hmmp_memop_overflow(model.num_states, max_length)){
		*o_logP_arr = 0;
		return e_overflow;
	}
	prob_arr = hmmp_create_dbl_array ( num_obs );
	if ( ! prob_arr ){
		return E_ALLOCATION;
	}
#pragma omp parallel firstprivate(alfa,scales) default(shared) num_threads(HMMP_NUM_THREADS)
{
	#pragma omp critical 
	{
		#pragma omp flush ( fail_flag )
		if(!fail_flag){
			alfa = hmmp_create_dbl_matrix(model.num_states*max_length);
			if( ! alfa )
				fail_flag = 1; 
			else{
				scales = hmmp_create_dbl_array ( max_length );
				if(!scales)
					fail_flag = 1; 
			}
		}
	}//end of critical
	#pragma omp barrier
	#pragma omp flush ( fail_flag )
	if(!fail_flag){
		#pragma omp for private(k) schedule(static) nowait
		for ( k = 0 ; k < num_obs ; ++k ){
			hmmp_forward_alg(model,observ_arr[k],alfa,scales);
			prob_arr[k] = hmmp_log_of_divisors(scales,observ_arr[k].length);
		}
	}
	#pragma omp critical
	{
		if ( alfa ) 
			hmmp_delete_dbl_matrix(alfa);
		if ( scales ) 
			hmmp_delete_dbl_array(scales);
	}
}//end of parallel region
	if ( fail_flag ){
		o_logP_arr = 0;
		hmmp_delete_dbl_array(prob_arr);
		return E_ALLOCATION;
	}
	*o_logP_arr = prob_arr;
	return E_SUCCESS;
}

int hmmp_baum_welch ( hmmp_Model *model,hmmp_Sequence *seq_arr, int num_seq , int max_steps, double delta_p )
{
	dbl_matrix *alfa = 0, *beta = 0, *gamma, *xi = 0;
	dbl_array *scales_a = 0, *scales_b = 0;
	dbl_matrix *a_num = 0, *b_num = 0;
	dbl_array *a_denom = 0, *b_denom = 0, *pi_new = 0;

	int i, j, k, t, e_overflow, max_length = 0;
	double logP_current, pi_normalize = 0.0;
	char flag_failed = 0;

	if (!model || !seq_arr )
		return E_PARAMETER;
	for ( i = 0 ; i < num_seq ; ++i )
		if ( max_length < seq_arr[i].length )
			max_length = seq_arr[i].length;
	if (e_overflow = hmmp_memop_overflow(model->num_states, max_length))
		return e_overflow;
	pi_new = hmmp_create_dbl_array ( model->num_states );
	if ( !pi_new ) {flag_failed = 1; goto SKIP_REST; }
	alfa = hmmp_create_dbl_matrix ( max_length * model->num_states );
	if ( !alfa ) {flag_failed = 1; goto SKIP_REST; }
	beta = hmmp_create_dbl_matrix ( max_length * model->num_states );
	if ( !beta ) {flag_failed = 1; goto SKIP_REST; }
	gamma = hmmp_create_dbl_matrix ( max_length * model->num_states );
	if ( !gamma ) {flag_failed = 1; goto SKIP_REST; }
	i = (max_length-1) * model->num_states * model->num_states;
	xi = hmmp_create_dbl_matrix ( (max_length-1) * model->num_states * model->num_states);
	if ( !xi ) {flag_failed = 1; goto SKIP_REST; }
	
	scales_a = hmmp_create_dbl_array ( max_length );
	if ( !scales_a ) {flag_failed = 1; goto SKIP_REST; }
	scales_b = hmmp_create_dbl_array ( max_length );
	if ( !scales_b ) {flag_failed = 1; goto SKIP_REST; }
	
	a_num = hmmp_create_dbl_matrix(model->num_states*model->num_states);
	if ( !a_num ) {flag_failed = 1; goto SKIP_REST; }
	b_num = hmmp_create_dbl_matrix(model->num_states*model->num_symbols);
	if ( !b_num ) {flag_failed = 1; goto SKIP_REST; }
	a_denom = hmmp_create_dbl_array(model->num_states);
	if ( !a_denom ) {flag_failed = 1; goto SKIP_REST; }
	b_denom = hmmp_create_dbl_array(model->num_states);
	if ( !b_denom ) {flag_failed = 1; goto SKIP_REST; }
SKIP_REST:
	if(!flag_failed)
	for ( t = 0 ; t < max_steps ; ++t ){
		logP_current = 0.0;
		hmmp_init_dbl_dataset(pi_new,model->num_states,0.0);
		hmmp_init_dbl_dataset(a_num,model->num_states*model->num_states,0.0);
		hmmp_init_dbl_dataset(b_num,model->num_states*model->num_symbols,0.0);
		hmmp_init_dbl_dataset(a_denom,model->num_states,0.0);
		hmmp_init_dbl_dataset(b_denom,model->num_states,0.0);
		for ( k = 0 ; k < num_seq ; ++k ){
#pragma omp parallel num_threads(HMMP_NUM_THREADS) default(shared)
		{
		#pragma omp sections
			{
			#pragma omp section
				{
				hmmp_forward_alg(*model, seq_arr[k], alfa, scales_a );
				logP_current += hmmp_log_of_divisors( scales_a, seq_arr[k].length );
				}
			#pragma omp section
				{ hmmp_backward_alg(*model, seq_arr[k], beta,scales_b ); }
			}
			#pragma omp single
			{ hmmp_backward_rescale(beta,model->num_states,seq_arr[k].length,scales_a,scales_b);}
			hmmp_bwa_gamma_alg(gamma,alfa,beta,scales_a,model->num_states,seq_arr[k].length);
			hmmp_bwa_xi_alg(xi,alfa,beta,*model,seq_arr[k]);
			#pragma omp barrier
			hmmp_bwa_reest_alg(*model,seq_arr[k],xi,gamma,a_num,b_num,a_denom,b_denom);
		}
			for ( i = 0 ; i < model->num_states ; ++ i ){
				pi_new[i] += gamma[i*seq_arr[k].length];
			}
		}
		if ( hmmp_delta_logp ( model->prior, logP_current ) < delta_p ){
			model->prior = logP_current;
			break;
		}
		model->prior = logP_current;
		if ( logP_current > 0.0-HMMP_PRECISION ){
			model->prior = 0.0;
			break;
		}
		for ( i = 0 ; i < model->num_states ; ++i )
			model->initial[i] = pi_new[i];
		hmmp_normalize_arr( model->initial, model->num_states );

		for ( i = 0 ; i < model->num_states ; ++i ){
			for ( j = 0 ; j < model->num_states ; ++j ){
				model->transition[i*model->num_states+j] = a_num[i*model->num_states+j] / a_denom[i];
			}
			for ( j = 0 ; j < model->num_symbols ; ++j ){
				model->emission[i*model->num_symbols+j] = b_num[i*model->num_symbols+j] / b_denom[i];
			}
		}
	}
	if (b_denom) hmmp_delete_dbl_array(b_denom);
	if (a_denom) hmmp_delete_dbl_array(a_denom);
	if (b_num) hmmp_delete_dbl_matrix(b_num);
	if (a_num) hmmp_delete_dbl_matrix(a_num);
	if (scales_b) hmmp_delete_dbl_array(scales_b);
	if (scales_a) hmmp_delete_dbl_array(scales_a);
	if (xi) hmmp_delete_dbl_matrix(xi);
	if (gamma) hmmp_delete_dbl_matrix(gamma);
	if (beta) hmmp_delete_dbl_matrix(beta);
	if (alfa) hmmp_delete_dbl_matrix(alfa);
	if (pi_new) hmmp_delete_dbl_array(pi_new);
	if (flag_failed)
		return E_ALLOCATION;
	return t;
}