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
#include "hmmp_algwrap.h"
#include "hmmp_datatypes.h"
#include "hmmp_memop.h"
#include "hmmp_alg.h"
#include "hmmp_dataproc.h"
int hmmp_forward (	hmmp_Model *model, hmmp_Sequence seq, dbl_matrix **o_addr_alfa,
					dbl_array **o_scale_arr )
{
	dbl_matrix *alfa;
	dbl_array *scaling;
	int e_overflow;
	
	if(!model || !o_addr_alfa || !o_scale_arr )
		return E_PARAMETER;
	if (e_overflow = hmmp_memop_overflow(model->num_states, seq.length))
		return e_overflow;
	scaling = hmmp_create_dbl_array(seq.length);
	if(!scaling)
		return E_ALLOCATION;
	alfa = hmmp_create_dbl_matrix(seq.length*model->num_states);
	if(!alfa){
		hmmp_delete_dbl_array(scaling);
		return E_ALLOCATION;
	}
	hmmp_forward_alg(*model,seq,alfa,scaling);
	*o_scale_arr = scaling;
	*o_addr_alfa = alfa;
	return E_SUCCESS;
}

int hmmp_backward(	hmmp_Model *model, hmmp_Sequence seq, dbl_matrix **o_addr_beta,
					dbl_array **o_scale_arr)
{
	dbl_matrix *beta;
	dbl_array *scaling;
	int e_overflow;

	if(!model || !o_addr_beta || !o_scale_arr )
		return E_PARAMETER;
	if (e_overflow = hmmp_memop_overflow(model->num_states, seq.length))
		return e_overflow;
	scaling = hmmp_create_dbl_array(seq.length);
	if(!scaling)
		return E_ALLOCATION;
	beta = hmmp_create_dbl_matrix(seq.length*model->num_states);
	if(!beta){
		hmmp_delete_dbl_array(scaling);
		return E_ALLOCATION;
	}
	hmmp_backward_alg(*model,seq,beta,scaling);
	*o_scale_arr = scaling;
	*o_addr_beta = beta;
	return E_SUCCESS;
}

int hmmp_viterbi(	hmmp_Model log_model,  hmmp_Sequence obs_seq, hmmp_Sequence *o_state_seq,
				 dbl_array *o_logPS )
{
	int_matrix *backtrack; // T x N matrix
	dbl_matrix *mu;	// Only 2 x N matrix
	double logP_recv;
	int e_overflow;

	if(!o_state_seq || !o_logPS)
		return E_PARAMETER;
	if (e_overflow = hmmp_memop_overflow(log_model.num_states, obs_seq.length))
		return e_overflow;
	backtrack = hmmp_create_int_matrix ( log_model.num_states * obs_seq.length );
	if(!backtrack)
		return E_ALLOCATION;
	mu = hmmp_create_dbl_matrix ( 2 * log_model.num_states );
	if(!mu){
		hmmp_delete_int_matrix ( backtrack );
		return E_ALLOCATION;
	}
	o_state_seq->seq_id = obs_seq.seq_id;
	o_state_seq->length = obs_seq.length;
	o_state_seq->cardinality = log_model.num_states;

	hmmp_viterbi_alg ( log_model, obs_seq, backtrack, mu, o_state_seq, &logP_recv );
	*o_logPS = logP_recv;

	hmmp_delete_dbl_matrix ( mu );
	hmmp_delete_int_matrix ( backtrack );

	return E_SUCCESS;
}
