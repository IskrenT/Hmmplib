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
#ifndef HMMP_DATAPROC_H
#define HMMP_DATAPROC_H
/** @file
*	This hmmp_dataproc.h contains small data computations and transformations used
*	throughout the library.
*/
#include "hmmp_datatypes.h"

/// Transpose a matrix.
/** 
*	@param[in,out] mat Address of the matrix to be transposed
*	@param[in]	num_rows Number of rows
*	@param[in]	num_cols Number of columns
*	@return @ref hmmp_Error Error code.
*/
int hmmp_transpose_matrix ( dbl_matrix *mat, int num_rows, int num_cols );

/// Normalize the values in each column of a matrix to sum up to 1.0.
/** 
*	@param[in,out] mat Address of the matrix to be modified
*	@param[in]	num_rows Number of rows
*	@param[in]	num_cols Number of columns
*	@return @ref hmmp_Error Error code.
*/
int hmmp_normalize_columns ( dbl_matrix *mat, int num_rows, int num_cols );

/// Normalize the values in each row of a matrix to sum up to 1.0.
/** 
*	@param[in,out] mat Address of the matrix to be modified
*	@param[in]	num_rows Number of rows
*	@param[in]	num_cols Number of columns
*	@return @ref hmmp_Error Error code.
*/
int hmmp_normalize_rows ( dbl_matrix *mat, int num_rows, int num_cols );

/// Initialize a double dataset with the same value for all elements.
/** 
*	@param[in,out] p Address of the array to be initialized
*	@param[in]	size Number of values in the array
*	@param[in]	value The value to be assigned to each member of the array
*	@return @ref hmmp_Error Error code..
*/
int hmmp_init_dbl_dataset ( double *p, int size, double value );

/// Normalize an array of real numbers to 1.0.
/** 
*	This function is used to apply the scaling at the model parameters during
*	the [forward](@ref hmmp_forward_alg()) and [backward](@ref hmmp_backward_alg()) 
*	algorithms.
*	It is also used after reestimation of the initial parameters for multiple sequences
*	in the learning hmmp_baum_welch() algorithm.
*	
*	@param[in,out]	arr Address of the array of real numbered values to be manipulated
*	@param[in]	num_items Number of values in the array
*	@return the normalizing constant applied to the array
*/
double hmmp_normalize_arr ( dbl_array *arr, int num_items );

/// Calculate the differance between two logarithmic probabilities.
/** 
*	This function is used to calculate the distance between two
*	logarithmic probabilities. Since there isn't an unique quantative
*	measure of this distance, it is suggested to choose a definition
*	depending on the needs of application. For some models and
*	systems of models this particular implementation of the distance function
*	might be inappropriate. In this case a replacement of this function is
*	recommended.
*
*	This function is implemented according to the following definition:
*
*	log(delta P) = |log(new P) - log(old P)| = |log(new P / old P)|
*
*	This function is only used in the hmmp_baum_welch() algorithm, to determine
*	a breakpoint at which to stop repeating the algorithm.
*	
*	@param[in]	logP_old Prior logarithmic probability for the model
*	@param[in]	logP_new New obtained probability for the model
*	@return the module of the logarithm of the ratio between the two probabilities
*	or module of the numerical distance between the two logarithmic probabilities
*/
double hmmp_delta_logp ( double logP_old, double logP_new );

/// Calculate the logarithm  of 1/(d1*d2*...dn)
/** 
*	This function is used to calculate the logarithm of multiple divisors
*	using the following formula:
*
*	log(1/(d1*d2*...dn)) = -log(d1) -log(d2) ... -log(dn)
*
*	It can be used to calculate the model logarithmic probability
*	from the scaling factors of either the hmmp_forward_alg() or hmmp_backward_alg().
*	A renamed wrapper function that does this is - hmmp_model_logprobability()
*
*	@param[in]	divisors Address of a double array containing the divisors
*	@param[in]	num_div Number of divisors
*	@return The desired logarithmc value
*/
double hmmp_log_of_divisors ( double* divisors, int num_div );

/// Change data values to logarithmic scale
/** 
*	This function uses the 'log()' function from <math.h>.
*	It tansforms the data into logarithmic scale at base 'e'
*	natural logarithm.
*
*	When a value is in the inteval [ 0 ; HMMP_PRECISION )
*	it is transformed to: -HMMP_DBL_MAX So in logarithmic
*	scale it approaches 0 from above.
*
*	@param[in]	data The data to be manipulated
*	@param[in]	count The number of data items
*	@return @ref hmmp_Error Error code.
*/
int hmmp_data_log_scale ( double *data, int count );

/// Change the model parameters to logarithmic scale
/** 
*	This function uses hmmp_data_log_scale().
*	It tansforms the model parameters to logarithmic scale at base 'e'
*	natural logarithm.
*
*	When a value is in the inteval [ 0 ; HMMP_PRECISION )
*	it is transformed to: -HMMP_DBL_MAX So in logarithmic
*	scale it approaches 0 from above.
*
*	@param[in]	m The model to be manipulated.
*	@return @ref hmmp_Error Error code.
*/
int hmmp_model_log_param ( hmmp_Model *m );

/// Check if multiplication between unsigned integers will overflow.
/** 
*	This function is used extensively troughout the library to sequre the input to memory
*	allocation modules.
*	@param[in]	arr Array of the integers to be multiplied.
*	@param[in]	count Number of integers in the array.
*	@return If multiplication overflows: 1 'True', if not: 0 'False'
*/
int hmmp_multiplication_overflow(unsigned *arr, int count);

/// Create a copy of a model into an empty model
/** 
*	**Note:** Before using this function an empty model with the same number of states
*	and symbols should be created using: hmmp_create_model()
*	alternatively instead of this function, hmmp_create_model_copy() can be used.
*
*	@param[in]	dest address of destination model
*	@param[in]	source address of source model
*	@return @ref hmmp_Error Error code.
*/
int hmmp_model_copy ( hmmp_Model *dest, hmmp_Model *source );

/// Obtain the model probability given a sequence with the result from the forward algorithm.
/**
*	The scaling factors are result of the hmmp_forward_alg() or hmmp_forward() and are
*	valid for one sequence and one model.
*	The scaling factors should be stored in a double array of T size, where T is the length
*	of the observed sequence. Using the scaling factors the probability of the model given
*	the sequence can be obtained.
*
*	@param[in] forward_scale Address of the array containing the scaling factors.
*	@param[in] length			The length of the observed sequence and the array
*	@return The logarithm of the probability for the model to produce the sequence.
*			0.0 on invalid input.
*/
double hmmp_model_logprobability ( dbl_array *forward_scale, int length );



#endif