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
#include "hmmp_dataproc.h"
#include "hmmp_datatypes.h"
#include "hmmp_memop.h"
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>

const double HMMP_DBL_MAX = 1.7976931348623157e+308;
const double HMMP_PRECISION = 2.2204460492503131e-016;

int hmmp_transpose_matrix ( dbl_matrix *mat, int num_rows, int num_cols )
{
	int i, j;
	dbl_matrix *mat2;
	if ( !mat )
		return E_PARAMETER;
	mat2 = hmmp_create_dbl_matrix(num_rows*num_cols);
	if(!mat2)
		return E_ALLOCATION;
	for ( i = 0 ; i < num_rows ; ++i )
		for ( j = 0 ; j < num_cols ; ++j )
			mat2[j*num_rows+i]=mat[i*num_cols+j];
	for ( i = 0 ; i < num_rows*num_cols ; ++i )
			mat[i]=mat2[i];
	hmmp_delete_dbl_matrix(mat2);
	return E_SUCCESS;
}

int hmmp_normalize_columns ( dbl_matrix *mat, int num_rows, int num_cols )
{
	int i, j;
	double partsum;
	if(!mat)
		return E_PARAMETER;
	for( i = 0 ; i < num_cols ; ++i ){
		partsum = 0.0;
		for( j = 0 ; j < num_rows ; ++j)
			partsum += mat[j*num_cols+i];
		partsum = 1 / partsum;
		for( j = 0 ; j < num_rows ; ++j)
			mat[j*num_cols+i]*=partsum;
	}
	return E_SUCCESS;
}

int hmmp_normalize_rows ( dbl_matrix *mat, int num_rows, int num_cols )
{
	int i, j;
	double partsum;
	if(!mat)
		return E_PARAMETER;
	for( i = 0 ; i < num_rows ; ++i ){
		partsum = 0.0;
		for( j = 0 ; j < num_cols ; ++j)
			partsum += mat[i*num_cols+j];
		partsum = 1.0 / partsum;
		for( j = 0 ; j < num_cols ; ++j)
			mat[i*num_cols+j]*=partsum;
	}
	return E_SUCCESS;
}

double hmmp_normalize_arr ( dbl_array *arr, int num_items )
{
	int i;
	double partsum = 0.0;
	if(!arr)
		return E_PARAMETER;
	for ( i = 0 ; i < num_items ; ++i )
		partsum += arr[i];
	partsum = 1.0 / partsum;
	for ( i = 0 ; i < num_items ; ++i )
		arr[i]*=partsum;
	return partsum;
}

int hmmp_init_dbl_dataset ( double *p, int size, double value )
{
	int i;
	if(!p)
		return E_PARAMETER;
	for ( i = size ; i-- ; )
		p[i] = value;
	return E_SUCCESS;
}

/// deltaP = log((Pnew-Pold)/Pnew))
double hmmp_delta_logp ( double logP_old, double logP_new )
{
	double ret;
	if(logP_old > logP_new ){
		ret=logP_old;
		logP_old=logP_new;
		logP_new=ret;
	}
	if(logP_old < HMMP_PRECISION && logP_old > -HMMP_PRECISION)
		return log(HMMP_DBL_MAX);
	ret = logP_new - logP_old;
	return ret;
}
/// log(1/(d1*d2*...dn)) = -log(d1) -log(d2) ... -log(dn)
/**
* No support for zero input of a divisor.
*/
double hmmp_log_of_divisors ( double *divisors, int num_div )
{
	double ret;
	int i;
	if(! divisors ) 
		return E_PARAMETER;
	ret = 0.0;
	for ( i = 0 ; i < num_div ; ++i ){
		if ( divisors[i]<HMMP_PRECISION && divisors[i]>-HMMP_PRECISION)
			return HMMP_DBL_MAX;
		ret -= log(divisors[i]);
	}
	return ret;
}
int hmmp_data_log_scale ( double *data, int count )
{
	int i;
	if(!data)
		return E_PARAMETER;
	for ( i = 0 ; i < count ; ++i ){
		if ( data[i] < HMMP_PRECISION )
			data[i] = -HMMP_DBL_MAX;
		else
			data[i] = log ( data[i] );
	}
	return E_SUCCESS;
}
int hmmp_model_log_param ( hmmp_Model *m )
{
	if(!m)
		return E_PARAMETER;
	hmmp_data_log_scale(m->initial,m->num_states);
	hmmp_data_log_scale(m->transition,m->num_states*m->num_states);
	hmmp_data_log_scale(m->emission,m->num_states*m->num_symbols);
	return E_SUCCESS;
}

int hmmp_multiplication_overflow(unsigned int *arr, int count){
	double dbl_test = 1.0;
	int	int_test = 1;
	int i;
	if ( !arr || count < 2 ){
		return E_ARGUMENT;
	}
	for ( i = 0 ; i < count ; ++i ){
		dbl_test *= (double) arr[i];
		int_test *= arr[i];
	}
	for( i = 1 ; i < count; ++i ){
		int_test /= arr[i];
	}
	if ( dbl_test > (double) (UINT_MAX-1) )
		return 1;
	return 0;
}
int hmmp_model_copy ( hmmp_Model *dest, hmmp_Model *source )
{
	if(!dest || !source )
		return E_PARAMETER;
	dest->model_id = source->model_id;
	dest->num_states = source->num_states;
	dest->num_symbols = source->num_symbols;
	dest->prior = source->prior;
	memcpy ( dest->initial, source->initial, source->num_states * sizeof ( dbl_array ) );
	memcpy ( dest->transition, source->transition,
			source->num_states*source->num_states * sizeof ( dbl_matrix ));
	memcpy ( dest->emission, source->emission, 
		source->num_states*source->num_symbols * sizeof ( dbl_matrix));
	return E_SUCCESS;
}
double hmmp_model_logprobability ( dbl_array *forward_scale, int length ){
	if ( !forward_scale || !length )
		return 0.0;
	return hmmp_log_of_divisors(forward_scale, length);
};