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
#include "hmmp_generate.h"
#include "hmmp_datatypes.h"
#include "hmmp_dataproc.h"
#include "hmmp_memop.h"
#include <stdlib.h>
#include <malloc.h>

hmmp_Model *hmmp_gen_from_seq ( hmmp_Sequence *seq )
{
	hmmp_Model *m;
	int i;
	if(!seq)
		return 0;
	m = hmmp_create_model(seq->length+1, seq->cardinality);
	if ( !m )
		return 0;
	for ( i = 0 ; i < m->num_states ; ++i )
		m->initial[i] = 0;
	for ( i = 0 ; i < m->num_states * m->num_states; ++i )
		m->transition[i] = 0;
	for ( i = 0 ; i < m->num_states * m->num_symbols ; ++i )
		m->emission[i] = 0;
	m->initial[0] = 1;
	for ( i = 0 ; i < m->num_states-1; ++i )
		m->transition[(i+1)*m->num_states + i ] = 1;
	for ( i = 0 ; i < m->num_symbols ; ++i )
		m->emission[i*m->num_symbols + seq->sequence[i]] = 1;
	m->prior = 1.0;
	return m;
}

hmmp_Model *hmmp_gen_random_models ( int count, int num_states, int num_symbols, int seed )
{
	int i, j;
	hmmp_Model *arr;
	srand(seed);
//allocate all models, and all matrices
	arr = hmmp_create_arr_models(count,num_states,num_symbols);
	if(!arr)
		return 0;
//initialize all models
	for ( i = 0 ; i < count ; ++i ){
		arr[i].model_id = i;
		arr[i].num_states = num_states;
		arr[i].num_symbols = num_symbols;
		arr[i].prior = 0.0;
		for ( j = 0 ; j < num_states ; ++j ){
			arr[i].initial[j] = 1/(double)(1+rand()%(num_states*num_symbols));
		}
		for ( j = 0 ; j < num_states*num_states ; ++j ){
			arr[i].transition[j] = 1/(double)(1+rand()%(num_states*num_symbols));
		}
		for ( j = 0 ; j < num_states*num_symbols ; ++j ){
			arr[i].emission[j] = 1/(double)(1+rand()%(num_states*num_symbols));
		}
		hmmp_normalize_arr(arr[i].initial,num_states);
		hmmp_normalize_rows(arr[i].transition,num_states,num_states);
		hmmp_normalize_rows(arr[i].emission,num_states,num_symbols);
	}
	return arr;
}

hmmp_Sequence *hmmp_gen_random_sequences ( int count, int cardinality, int length, int seed )
{
	int i, j;
	hmmp_Sequence *arr;
	srand(seed);
//allocate all sequences
	arr = hmmp_create_arr_seq(count, length, cardinality);
	if(!arr)
		return 0;
//initialize all sequences
	for ( i = 0 ; i < count ; ++i ){
		arr[i].seq_id = i;
		arr[i].length = length;
		arr[i].cardinality = cardinality;

		for ( j = 0 ; j < length ; ++j ){
			arr[i].sequence[j] = rand()%cardinality;
		}
	}
	return arr;
}
