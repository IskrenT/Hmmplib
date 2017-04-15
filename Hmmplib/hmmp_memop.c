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
#include "hmmp_memop.h"
#include "hmmp_datatypes.h"
#include "hmmp_dataproc.h"
#include <malloc.h>

hmmp_Model *hmmp_create_model ( int num_states, int num_symbols )
{
	hmmp_Model *model;
	int e_overflow = 0;
	unsigned int overflow[3];
	overflow[0] = num_states;
	overflow[1] = num_symbols;
	overflow[2] = sizeof(dbl_matrix);
	e_overflow |= hmmp_multiplication_overflow(overflow, 3);
	overflow[1] = num_states;
	e_overflow |= hmmp_multiplication_overflow(overflow, 3);
	if ( e_overflow )
		return 0;
	model = ( hmmp_Model*)malloc(sizeof(hmmp_Model));
	if(!model)
		return 0;
//init
	model->model_id = 0;
	model->num_states = num_states;
	model->num_symbols = num_symbols;
	if (!(model->initial = (dbl_array*) malloc ( sizeof(dbl_array) * num_states ))){
		free(model);
		return 0;
	}
	if (!(model->transition = ( dbl_matrix*) malloc ( sizeof(dbl_matrix) * num_states * num_states))){
		free(model->initial);
		free(model);
		return 0;
	}
	if (!(model->emission = (dbl_matrix*) malloc ( sizeof(dbl_matrix) * num_states * num_symbols ))){
		free(model->initial);
		free(model->transition);
		free(model);
		return 0;
	}
	return model;
}
hmmp_Model *hmmp_create_model_copy ( hmmp_Model *source )
{
	hmmp_Model *m;
	if (!source)
		return 0;
	m = hmmp_create_model(source->num_states,source->num_symbols);
	if(!m)
		return 0;
	hmmp_model_copy(m,source);
	return m;
}
int hmmp_delete_model ( hmmp_Model *model )
{
	if(!model)
		return E_PARAMETER;
	model->model_id=-1;
	model->num_states=0;
	model->num_symbols=0;
	model->prior=0.0;
	free(model->initial);
	free(model->transition);
	free(model->emission);
	free(model);
	return E_SUCCESS;
}
hmmp_Model *hmmp_create_arr_models ( int count, int num_states, int num_symbols )
{
	int i;
	hmmp_Model *arr;
	int e_overflow = 0;
	unsigned int overflow[3];
	overflow[0] = num_states;
	overflow[1] = num_symbols;
	overflow[2] = sizeof(dbl_matrix);
	e_overflow |= hmmp_multiplication_overflow(overflow, 3);
	overflow[1] = num_states;
	e_overflow |= hmmp_multiplication_overflow(overflow, 3);
	overflow[0] = count;
	overflow[1] = sizeof(hmmp_Model);
	e_overflow |= hmmp_multiplication_overflow(overflow, 2);
	if ( e_overflow )
		return 0;
//allocate all models, and all matrices
	arr = (hmmp_Model*)malloc(sizeof(hmmp_Model)*count);
	if(!arr)
		return 0;
	arr = (hmmp_Model*)malloc(sizeof(hmmp_Model)*count);
	if(!arr)
		return 0;
	for ( i = 0 ; i < count ; ++i ){
		arr[i].num_states = num_states;
		arr[i].num_symbols = num_symbols;
		arr[i].initial = (dbl_array*)malloc(sizeof(dbl_array)*num_states);
		if(!(arr[i].initial))
			break;
		arr[i].transition = (dbl_matrix*)malloc(sizeof(dbl_matrix)*num_states*num_states);
		if(!(arr[i].transition)){
			free(arr[i].initial);
			break;
		}
		arr[i].emission = (dbl_matrix*)malloc(sizeof(dbl_matrix)*num_states*num_symbols);
		if(!(arr[i].emission)){
			free(arr[i].initial);
			free(arr[i].transition);
			break;
		}
	}
	if(i!=count){
		while(i+1){
			free(arr[i].initial);
			free(arr[i].transition);
			free(arr[i].emission);
			--i;
		}
		free(arr);
		return 0;
	}
	return arr;
}
int hmmp_delete_arr_models(hmmp_Model *arr, int count)
{
	int i;
	if(!arr)
		return E_PARAMETER;
	for ( i = 0 ; i < count ; ++i ){
		arr[i].model_id=-1;
		free(arr[i].initial);
		free(arr[i].transition);
		free(arr[i].emission);
	}
	free(arr);
	return E_SUCCESS;
}

dbl_array *hmmp_create_dbl_array(unsigned num_elements)
{
	dbl_array *arr;
	arr = ( dbl_array * ) malloc ( sizeof(dbl_array) * num_elements );
	return arr;
}
int hmmp_delete_dbl_array(dbl_array* arr)
{
	free(arr);
	return E_SUCCESS;
}
dbl_matrix *hmmp_create_dbl_matrix(unsigned num_elements)
{
	dbl_matrix *mat;
	mat = ( dbl_matrix * ) malloc ( sizeof(dbl_matrix)*num_elements);
	return mat;
}
int hmmp_delete_dbl_matrix(dbl_matrix* mat)
{
	if(!mat)
		return E_PARAMETER;
	free(mat);
	return E_SUCCESS;
}
int_array *hmmp_create_int_array(unsigned num_elements)
{
	int_array *arr;
	arr = (int_array*) malloc ( sizeof(int_array) * num_elements );
	return arr;
}
int hmmp_delete_int_array(int_array* arr)
{
	if(!arr)
		return E_PARAMETER;
	free(arr);
	return E_SUCCESS;
}
int_matrix *hmmp_create_int_matrix(unsigned num_elements)
{
	int_matrix *arr;
	arr = (int_matrix*) malloc ( sizeof(int_matrix) * num_elements );
	return arr;
}
int hmmp_delete_int_matrix(int_matrix* mat)
{
	if(!mat)
		return E_PARAMETER;
	free(mat);
	return E_SUCCESS;
}
hmmp_Sequence *hmmp_create_sequence( int length )
{
	hmmp_Sequence *seq;
	seq = ( hmmp_Sequence* ) malloc ( sizeof(hmmp_Sequence) );
	seq->sequence = hmmp_create_int_array(length);
	if(!(seq->sequence)){
		free(seq);
		return 0;
	}
	seq->length = length;
	return seq;
}
int hmmp_delete_sequence(hmmp_Sequence *seq)
{
	if ( !seq )
		return E_PARAMETER;
	free(seq->sequence);
	seq->cardinality=0;
	seq->length=0;
	seq->seq_id=-1;
	free(seq);
	return E_SUCCESS;
}

hmmp_Sequence *hmmp_create_arr_seq(int count, int length, int num_symbols)
{
	int i;
	hmmp_Sequence *arr;
	int e_overflow = 0;
	unsigned int overflow[2];
	overflow[0] = count;
	overflow[1] = sizeof(hmmp_Sequence);
	e_overflow |= hmmp_multiplication_overflow(overflow, 2);
	overflow[0] = length;
	overflow[1] = sizeof(int_array);
	e_overflow |= hmmp_multiplication_overflow(overflow, 2);
	if ( e_overflow )
		return 0;
	//allocate all sequences
	arr = (hmmp_Sequence*)malloc(sizeof(hmmp_Sequence)*count);
	if(!arr)
		return 0;
	for ( i = 0 ; i < count ; ++i ){
		arr[i].sequence = (int_array*)malloc(sizeof(int_array)*length);
		if(!(arr[i].sequence))
			break;
		arr[i].seq_id = i;
		arr[i].length = length;
		arr[i].cardinality = num_symbols;
	}
	if(i!=count){
		while(i+1){
			free(arr[i].sequence);
			--i;
		}
		free(arr);
		return 0;
	}
	return arr;
}
int hmmp_delete_arr_seq(hmmp_Sequence *arr, int count)
{
	int i;
	if(!arr)
		return E_PARAMETER;
	for ( i = 0 ; i < count ; ++i ){
		arr[i].seq_id=-1;
		arr[i].length=0;
		free(arr[i].sequence);
	}
	free(arr);
	return E_SUCCESS;
}

int hmmp_memop_overflow(int num_states, int seq_length){
	unsigned int ovrfl[4];
	ovrfl[0] = num_states;
	ovrfl[1] = sizeof(double);
	ovrfl[2] = seq_length;
	ovrfl[3] = num_states;
	if(hmmp_multiplication_overflow(ovrfl, 3))
		return E_MEM_OVERFLOW_L1;
	if(hmmp_multiplication_overflow(ovrfl, 4))
		return E_MEM_OVERFLOW_L2;
	return E_SUCCESS;

}