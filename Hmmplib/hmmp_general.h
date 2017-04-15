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
#ifndef HMMP_GENERAL_H
#define HMMP_GENERAL_H
/** @file
*	This file contains declaration for the general high-abstraction solutions of HMM.
*	Check the examples in @ref general for detailed explaination how to use the functions
*	declared in this module.
*	Here you can find solutions for:
*		- Evaluating multiple models against a single sequence.
*		- Evaluating a single model against multiple sequences.
*		- Decoding multiple sequences with a single model.
*		- Learning with multiple sequences.
*
*	This module is implemented with different levels of parallelism using OpenMP
*	If you wish to utilize the parallel implementations make sure OpenMP is available and
*	and turned on in your environment configuration settings.
*	The OpenMP header must be available and allowed to be included in your project:
*	**omp.h** must be available for linking with **hmmp_general.c** and **hmmp_alg.c**!
*/
#include "hmmp_datatypes.h"
/// Use Hmmplib for decoding with a single model and multiple sequences.
/**
*	This function is a general high-abstraction solution of the decoding problem in 
*	Hidden Markov Model ( HMM ) implemented using the Viterbi algoirthm.
*	For more information on what is decoding and examples of how to use the function 
*	check: @ref decoding
*
*	This function operates on multiple sequences and a single model.
*	Multi-core domain decomposition parallelism is implemented using OpenMP.
*	To specify the desired number of threads for the algorithm change
*	the global variable HMMP_NUM_THREADS.
*	
*	**Note:** After using the results of the function, use the appropriate functions from
*	hmmp_memop.h to delete ( deallocate ) the data structures holding the reuslts.
*
*	@param[in] model The model to operate with
*	@param[in] observ_array Adress of an array of observable sequences to operate with
*	@param[in] num_obs Number of observable sequences in the array
*	@param[out] o_state_array	Address of a pointer to receive the resulting array of 
								state sequences
*	@param[out] o_logPS	Address of a pointer to receive the resulting array of probabilities.
*						Each element is the logarithmic probability of the state sequence to 
						emit the corresponding observable sequence.
*	@return @ref hmmp_Error Error code.
*/
int hmmp_decode(hmmp_Model model, hmmp_Sequence *observ_array, int num_obs, 
				hmmp_Sequence **o_state_array, dbl_array **o_logPS );

/// Use Hmmplib for evaluation with multiple models and a single sequence.
/**
*	This function is a general high-abstraction solution of the evaluation problem in 
*	Hidden Markov Model ( HMM ) implemented using the Forward algoirthm.
*	For more information on what is evaluation and examples how to use the function 
*	check: @ref evaluation
*
*	This function operates on multiple models and a single sequence.
*	Multi-core domain decomposition parallelism is implemented using OpenMP.
*	To specify the desired number of threads for the algorithm change
*	the global variable HMMP_NUM_THREADS.
*	
*	**Note:** After using the results of the function, use the appropriate functions from
*	hmmp_memop.h to delete ( deallocate ) the data structures holding the reuslts.
*
*	@param[in] arr_models Address of an array of models to operate with
*	@param[in] num_models Number of models in the array
*	@param[in] observ_seq The observable sequence to operate with
*	@param[out] o_logP	Address of a pointer to receive the resulting array of logarithmic
*						probabilities. Each probability corresponds to a model and is
*						equal to the probability of this model to produce the sequence.
*	@return @ref hmmp_Error Error code.	
*/
int hmmp_evaluate_models(hmmp_Model *arr_models, int num_models, 
						 hmmp_Sequence observ_seq, dbl_array **o_logP );

/// Use Hmmplib for evaluation with a single models and multiple sequence.
/**
*	This function is a general high-abstraction solution of the evaluation problem in 
*	Hidden Markov Model ( HMM ) implemented using the Forward algoirthm.
*	For more information on what is evaluation and examples how to use the function 
*	check: @ref evaluation
*
*	This function operates on a single models and multiple sequences.
*	Multi-core domain decomposition parallelism is implemented using OpenMP.
*	To specify the desired number of threads for the algorithm change
*	the global variable HMMP_NUM_THREADS.
*	
*	**Note:** After using the results of the function, use the appropriate functions from
*	hmmp_memop.h to delete ( deallocate ) the data structures holding the reuslts.
*
*	@param[in] model		The model to operate with
*	@param[in] observ_arr Adress of an array of observable sequences to operate with
*	@param[in] num_obs		Number of observable sequences in the array
*	@param[out] o_logP	Address of a pointer to receive the resulting array of logarithmic
*						probabilities. Each probability corresponds to a sequence and is
*						equal to the probability of the model to produce this sequence.
*	@return @ref hmmp_Error Error code.	
*/
int hmmp_evaluate_sequences(hmmp_Model model, hmmp_Sequence *observ_arr, int num_obs,
							dbl_array **o_logP );


/// Use Hmmplib for learning with a single model and multiple sequences.
/**
*	This function is a general high-abstraction solution of the learning problem in 
*	Hidden Markov Model ( HMM ) implemented using the Baum-Welch algorithm.
*	For more information on what is learning and examples how to use the function 
*	check: @ref learning
*
*	This function operates on a single model and multiple sequences.
*	For this setup, the reestimation is done using the modified advanced reestimation
*	formulas, including multiple observed sequences. Traditionally, in this case the
*	initial parameters of the model are not reestimated, but in the current implementation
*	a normalized summation over all reestimated initial parameters is done.
*
*	Multi-core parallelism is implemented using OpenMP. The parallelism is of type
*	task decomposition in sections of the Baum-Welch algorithm. The recommended number
*	of threads is 2 so that the two highest complexity concurrent tasks in the algorithm
*	can work independently, namely the [Forward](@ref hmmp_forward_alg()) and the 
*	[Backward](@ref hmmp_backward_alg()) algorithms.
*	Optionally more threads can be used since the later stages of the Baum-Welch algorithm
*	use domain decomposition ( hmmp_bwa_gamma_alg() , hmmp_bwa_xi_alg() , hmmp_bwa_reest_alg() )
*
*	To specify the desired number of threads for the algorithm change
*	the global variable HMMP_NUM_THREADS ( 2 to (number of cores) recommended ).
*	
*	**Note:** The model parameters will change after executing this function.
*
*	@param[in,out] model The address of the model to operate on ( model paramaters will change )
*	@param[in] seq_arr	Adress of an array of observable sequences to operate with
*	@param[in] num_seq	Number of observable sequences in the array
*	@param[in] max_steps Maximum number of times to repeat the Baum-Welch algorithm
*	@param[in] delta_p	Minimum differance between the old and the new probability of the model,
						above which, the algorithm will keep executing. See hmmp_delta_logp().
*	@return @ref hmmp_Error Error code.	
*/
int hmmp_baum_welch ( hmmp_Model *model,hmmp_Sequence *seq_arr, int num_seq , int max_steps, double delta_p );

#endif