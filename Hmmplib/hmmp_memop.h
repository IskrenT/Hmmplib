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
#ifndef HMMP_MEMOP_H
#define HMMP_MEMOP_H
/** @file
*	This header includes declarations for the memory operation functions.
*	These functions can be used to safely create, initialize and destroy instances of:
*	a model, a sequence, multiple models and sequences, matrices, arrays.
*	Check @ref datafilelayout "Data Layout" to for reference regarding user defined types
*	and their usage.
*/
#include "hmmp_datatypes.h"
/// Create an empty model.
/**
*	In order to use the model it has to be initialized. Function declarations for 
*	initialization can be found in hmmp_generate.h and hmmp_file.h.
*	Always delete model after use with the model delete function listed below.
*	@param[in] num_states	Number of possible states in the model.
*	@param[in] num_symbols	Number of possible symbols emitted by the model.
*
*	@return	Address of the new model in heap. Zero 0 ( NULL ) on failure.
*	
*	@see hmmp_delete_model()
*/
hmmp_Model *hmmp_create_model ( int num_states, int num_symbols );

/// Safely delete a model.
/**
*	@param[in] model Adress of the previously created model to be deleted.
*	@return @ref hmmp_Error Error code..
*	@see hmmp_create_model()
*/
int hmmp_delete_model ( hmmp_Model *model );

/// Create an copy of an existing model.
/**
*	Always delete model after use with the model delete function listed below.
*
*	@param[in] source	The model to be coppied.
*	@return				Address of the new model in heap. Zero 0 ( NULL ) on failure.
*
*	@see hmmp_delete_model()
*/
hmmp_Model *hmmp_create_model_copy ( hmmp_Model *source );

/// Create an array of emtpy models with the same number of states and symbols.
/**
*	Always delete models after use with the models delete function listed below.
*	
*	@param[in] count		Number of models to be created.
*	@param[in] num_states	Number of possible states each model can have.
*	@param[in] num_symbols	Number of possible symbols emitted by the models.
*	@return		Address of the new model in heap. Zero 0 ( NULL ) on failure.
*
*	@see hmmp_delete_arr_models()
*/
hmmp_Model *hmmp_create_arr_models ( int count, int num_states, int num_symbols );

/// Delete an array of previously created models.
/** 
*	@param[in]	arr		Address of the array to be deleted.
*	@param[in]	count	Number of models in the array to be deleted.
*	@return	0 Success.
*/
int hmmp_delete_arr_models(hmmp_Model *arr, int count);

/// Create an uninitialized array ( vector ) of real double precision floating-point numbers.
/** 
*	@param[in]	num_elements Number of elements in the array.
*	@return Address of the new array in heap. Zero 0 ( NULL ) on failure.
*	@see hmmp_delete_dbl_array()
*/
dbl_array *hmmp_create_dbl_array(unsigned num_elements);

/// Delete an array ( vector ) of real double precision floating-point numbers.
/** 
*	@param[in]	arr	Address of the array to be deleted.
*	@return @ref hmmp_Error Error code..
*/
int hmmp_delete_dbl_array(dbl_array* arr);

/// Create an uninitialized matrix of real double precision floating-point numbers.
/** 
*	CAUTION: User responsibility to input proper matrix size and interpret the recieved
*	array as a matrix in the required dimentions. For information on how to interpret
*	the transition and emission matrices check: @ref datafilelayout "Data Layout"
*
*	@param[in]	num_elements Number of elements in the matrix.
*	@return Address of the new matrix in heap. Zero 0 ( NULL ) on failure.
*	@see hmmp_delete_dbl_matrix()
*/
dbl_matrix *hmmp_create_dbl_matrix(unsigned num_elements);

/// Delete a previously created matrix of real double precision floating-point numbers.
/** 
*	@param[in]	mat Address of the matrix to be deleted.
*	@return @ref hmmp_Error Error code..
*	@see hmmp_create_dbl_matrix()
*/
int hmmp_delete_dbl_matrix(dbl_matrix* mat);

/// Create an uninitialized array (vector) of integers.
/** 
*	@param[in]	num_elements Number of elements in the array.
*	@return Address of the new array in heap. Zero 0 ( NULL ) on failure.
*	@see hmmp_delete_int_array()
*/
int_array *hmmp_create_int_array(unsigned num_elements);

/// Delete a previously created array (vector) of integers.
/** 
*	@param[in]	arr Address of the array to be deleted.
*	@return @ref hmmp_Error Error code..
*	@see hmmp_create_int_array()
*/
int hmmp_delete_int_array(int_array* arr);

/// Create an uninitialized matrix of integers.
/** 
*	CAUTION: User responsibility to input proper matrix size and interpret the recieved
*	array as a matrix in the required dimentions. For information on how to interpret
*	the transition and emission matrices check: @ref datafilelayout "Data Layout"
*
*	@param[in]	num_elements Number of elements in the matrix.
*	@return Address of the new matrix in heap. Zero 0 ( NULL ) on failure.
*	@see hmmp_delete_int_matrix()
*/
int_matrix *hmmp_create_int_matrix(unsigned num_elements);

/// Delete a previously created matrix of integers.
/** 
*	@param[in]	mat Address of the matrix to be deleted.
*	@return @ref hmmp_Error Error code..
*	@see hmmp_create_int_matrix()
*/
int hmmp_delete_int_matrix(int_matrix* mat);

/// Create an empty sequence.
/**
*	In order to use the sequence it has to be initialized. Function declarations for 
*	initialization can be found in hmmp_generate.h and hmmp_file.h.
*	Always delete sequences after use with the sequence delete function listed below.
*
*	@param[in] length	Desired lenth of the new sequence
*	
*	@return	Address of the new sequence in heap. Zero 0 ( NULL ) on failure.
*	
*	@see hmmp_delete_sequence()
*/
hmmp_Sequence *hmmp_create_sequence( int length );

/// Safely delete a sequence.
/**
*	@param[in]	seq	Adress of the previously created sequence to be deleted.
*	@return @ref hmmp_Error Error code..
*	@see hmmp_create_model()
*/
int hmmp_delete_sequence(hmmp_Sequence *seq);

/// Create an array of empty sequences with the same length and symbol space
/**
*	Always delete sequences after use with the hmmp_memop.h function suggested below.
*	
*	@param[in] count		Number of sequences to be created
*	@param[in] length		The maximum length of each sequence
*	@param[in] num_symbols	Number of different possible symbols each sequence
*	@return		Address of the resulting array of sequences. Zero 0 ( NULL ) on failure
*
*	@see hmmp_delete_arr_seq()
*/
hmmp_Sequence *hmmp_create_arr_seq(int count, int length, int num_symbols);

/// Delete an array of previously created sequences.
/** 
*	@param[in]	arr		Address of the array to be deleted.
*	@param[in]	count	Number of sequences in the array to be deleted.
*	@return	0 Success.
*/
int hmmp_delete_arr_seq(hmmp_Sequence *arr, int count);

/// Check if overflow will occur in the multiplication of argument of memory allocation functions
/** 
*	When working with models with large state spaces and/or testing against very long sequences
*	overflow might occur when calling a memory allocation(creation) function from hmmp_memop.h
*	Generally the user is responsible to send a proper value as an argument to these functions.
*	Still the easiest way to use the allocation functions is simply calling them with the argument
*	being calculated in the call like so:
*	
*	    hmmp_create_...(states * sequence_length * sizeof(...) * ...);
*
*	the expression inside the argument of these calls can overflow. The allocation functions
*	does not provide means to check for such overflow. In order to check for a possible overflow
*	beforehand the following statement can be used before any of the allocation function calls:
*
*		if ( ret = hmmp_memop_overflow( number_of_states, sequence_length ) )
*			\\ ... HANDLE THE ERROR \\ ;
*			
*	'ret' will be 0 or 'False' when overflow can not occur with **proper usage** of the 
*	allocation functions.
*	
*
*	@param[in]	num_states	Number of states in the used model
*	@param[in]	seq_length	Length of the sequence you test the model against
*	@return	
*	Returns 0 or 'False' when overflow can not occur with **proper usage** of the 
*	allocation functions for the model and the sequence.
*	Returns E_MEM_OVERFLOW_L1 when the overflow will occur creating the variables:
*	'alfa', 'beta', 'gamma', 'xi'.
*	Returns E_MEM_OVERFLOW_L2 when the overflow will occur creating the variable 'xi'.
*/
int hmmp_memop_overflow(int num_states, int seq_length);

#endif