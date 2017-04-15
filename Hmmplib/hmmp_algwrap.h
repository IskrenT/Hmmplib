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
#ifndef HMMP_ALGWRAP_H
#define HMMP_ALGWRAP_H

/** @file
*	This file contains "wrappers" to the algorithmic solutions. These wrappers
*	allocate memory for the necessery internal variables for the algorithms.
*	Only use the plain algorithmic solutions if reusing the memory containers for
*	different steps is necessery. For the plain algorithmic solutions check:
*	 hmmp_alg.h . For the more general abstract parralel solutions to HMM check:
*	 hmmp_general.h 
*
*	The functions in this files are suitable for executing any of the algorithms
*	once. If multiple executions are needed consider using the solutions in
*	 hmmp_general.h
*
*	**Some of the algorithms needs pre-allocated containers for their output data,
*	please read the function and parameter descriptions carefully!**
*/

#include "hmmp_datatypes.h"

///Execute the forward algorithm on a model and a sequence including memory allocation.
/**
*	The forward algorithm is used to obtain the forward variables. From the forward
*	variables the probability of the model to have produced the sequence can be obtained.
*	For more information check @ref general .
*
*	To use this function pointers to a matrix and an array must be declared.
*	At the function call provide the addresses of these pointers.
*	After a return of the function these pointers will contain the addresses to the
*	resuting data structures. If the function failed to allocate memory for these data
*	structures, the values of the pointers will be 0 ( NULL ) and the return value
*	of the function will be an hmmp_Error code.
*
*	After using this function the user is responsible to destroy and free the memory
*	allocated to the resulting data structures. This can be done using the appropriate
*	functions declared in hmmp_memop.h .
*
*	@param[in] model	The model to execute algorithm on
*	@param[in] seq		The sequence to execute algorithm on
*	@param[out] o_addr_alfa	Address of a pointer that would receive the address of 
*							the resulting matrix of forward variables.
*	@param[out] o_scale_arr Address of a pointer that would receive the address of the 
*							 resulting array containing the scaling factors 
*	@return @ref hmmp_Error type
*	@see Use hmmp_forward_alg() instead if want to allocate the memory yourself.
*			Use hmmp_evaluate_sequences() if you wan't to use the forward algorithm
*			for evaluation.
*/
int hmmp_forward (	hmmp_Model *model, hmmp_Sequence seq, dbl_matrix **o_addr_alfa,
					dbl_array **o_scale_arr );

///Execute the backward algorithm on a model and a sequence including memory allocation.
/**
*	The backward algorithm is used to obtain the backward variables. It is mostly used
*	in the forward-backward algorithm and the Baum-Welch algorithm. It's result can also
*	be used to obtain the probability of the model to have produced the sequence.
*	For more information check @ref general .
*
*	To use this function pointers to a matrix and an array must be declared.
*	At the function call provide the addresses of these pointers.
*	After a return of the function these pointers will contain the addresses to the
*	resuting data structures. If the function failed to allocate memory for these data
*	structures, the values of the pointers will be 0 ( NULL ) and the return value
*	of the function will be an error code.
*
*	After using this function the user is responsible to destroy and free the memory
*	allocated to the resulting data structures. This can be done using the appropriate
*	functions declared in hmmp_memop.h .
*
*	@param[in] model	The model to execute algorithm on
*	@param[in] seq		The sequence to execute algorithm on
*	@param[out] o_addr_beta	Address of a pointer that would receive the address of 
*							the resulting matrix of backward variables.
*	@param[out] o_scale_arr Address of a pointer that would receive the address of the 
*							 resulting array containing the scaling factors 
*	@return @ref hmmp_Error type
*	@see	Use hmmp_backward_alg() instead if want to allocate the memory yourself.
*			Use hmmp_baum_welch() if you wan't to use the backward algorithm
*			for learining in the Baum-Welch algorithm.
*			For evaluation use hmmp_evaluate_sequences() that utilizes the 
*			forward algorithm.
*/
int hmmp_backward(	hmmp_Model *model, hmmp_Sequence seq, dbl_matrix **o_addr_beta,
					dbl_array **o_scale_arr);

///Execute the Viterbi algorithm on a model and a sequence including memory allocation.
/**
*	Viterbi's algorithm is used for decoding, or finding the best matching state sequence
*	to a sequence of observable symbols. For a configuration of a single model and 
*	multiple sequence use hmmp_decode() instead.
*
*	To use this function the input model should have logarithmic probabilities in the
*	parameters. Use hmmp_model_log_param() to logarithmize the model parameters.
*	Use a copy of the model if you don't want to manipulate the original model.
*	To create a copy hmmp_create_model_copy() can be used.
*
*	**Note: This function does not requre allocation of the internal data structures used
*			By the Viterbi algorithm, BUT it DOES REQUIRE allocation of the containers
*			that will hold the results.**
*
*	@param[in] log_model	The model with logarithmized parameters
*	@param[in] obs_seq		The sequence to execute algorithm on
*	@param[out] o_state_seq	Address of pre-allocated uninitialized sequence to hold result
*	@param[out] o_logP		Address of a double variable to hold the resulting probability
*	@return @ref hmmp_Error type
*	@see	Use hmmp_viterbi_alg() if you want to allocate internal variables yourself.
			Use hmmp_decode() if you want to use viterbi algorithm for decoding with
			configuration of multiple-sequences ( or one ) and one model.
*/
int hmmp_viterbi(	hmmp_Model log_model,  hmmp_Sequence obs_seq, hmmp_Sequence *o_state_seq,
				double *o_logP );

#endif