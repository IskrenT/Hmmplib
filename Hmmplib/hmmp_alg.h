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
#ifndef HMMP_ALG_H
#define HMMP_ALG_H
/** @file
*	This file contains only the algorithmic solutions, in order to use any of them
*	it is requred to have an existing initialized model and sequence, and to have
*	created ( allocated memory ) for the output variables using hmmp_memop.h.
*
*	**Failure to provide valid model, sequence, allocated containers for the result
*	is not handled by the algorithms and can result in runtime memory errors.**
*
*	These algorithms are unsafe and not suggested for general and computational use. 
*	Only use these functions if you are developing a new system of proccessing HMM or
*	if you wish to reuse the same containers for multiple calculations.
*
*	For general and computational use the functions declared in hmmp_general.h or 
*	hmmp_algwrap.h.
*/
#include "hmmp_datatypes.h"
///Execute the forward algorithm on a model and a sequence.
/**
*	**This function does not include memory allocation!!!**
*
*	The forward algorithm is used to obtain the forward variables. From the forward
*	variables the probability of the model to have produced the sequence can be obtained.
*	The dimentions of the matrix containing the forward variables should be:\n
*	N x T ( number of states x sequence length ).
*	For more information check @ref general .
*	
*	@param[in] model	The model to execute algorithm on
*	@param[in] seq		The sequence to execute algorithm on
*	@param[out] o_alfa	Address of pre-allocated matrix to store the resulting variables
*	@param[out] o_alfa_scale Adress of pre-allocated array to store the scaling factors 
*							 used at each step of the computation
*	@return @ref hmmp_Error Error code.
*	@see Use hmmp_forward() instead if you do not want to allocate memory yourself
*/
int hmmp_forward_alg (	hmmp_Model model, hmmp_Sequence seq, dbl_matrix *o_alfa,
						dbl_array *o_alfa_scale );
///Execute the backward algorithm on a model and a sequence.
/**
*	**This function does not include memory allocation!!!**
*
*	The backward algorithm is used to obtain the backward variables. It is mostly used
*	in the forward-backward algorithm and the Baum-Welch algorithm. It's result can also
*	be used to obtain the probability of the model to have produced the sequence.
*	The dimentions of the matrix containing the backward variables should be:\n
*	N x T ( number of states x sequence length ).
*	For more information check @ref general .
*
*	@param[in] model	The model to execute algorithm on
*	@param[in] seq		The sequence to execute algorithm on
*	@param[out] o_beta	Address of pre-allocated matrix to store the resulting variables
*	@param[out] o_beta_scale Adress of pre-allocated array to store the scaling factors 
*							 used at each step of the computation
*	@return @ref hmmp_Error Error code.
*	@see Use hmmp_backward() instead if you do not want to allocate memory yourself
*/
int hmmp_backward_alg ( hmmp_Model model, hmmp_Sequence seq, dbl_matrix *o_beta,
					    dbl_array *o_beta_scale );

///Execute the backward rescaling algorithm.
/**
*	This algorithm is meant to be used after the results from the backward and the forward
*	algorithms are available. It changes the backward variables scaling to match the
*	scaling of the forward algorithm. This is done to reduce computation of the next steps
*	in the Baum-Welch algorithm. This happens, because the terms cancel out with the
*	denominators in the formulas, namely the probability of the model to produce the sequence,
*	but they cancel out only if the same scaling factors for the forward and the backward
*	variables are present.
*
*	This function is part of the Baum-Welch algorithm used for learining: hmmp_baum_welch()
*
*	@param[in,out] beta		Address of the backward variables matrix
*	@param[in] num_states	Number of states in the model
*	@param[in] seq_len		Length of the observed sequence
*	@param[in] scale_alfa	Address of the array containing the scaling factors of the 
							forward variables
*	@param[in] scale_beta	Address of the array containing the scaling factors of the 
							backward variables
*	@return @ref hmmp_Error Error code.
*/
int hmmp_backward_rescale( dbl_matrix *beta, int num_states, int seq_len,
						   dbl_array *scale_alfa, dbl_array *scale_beta );

///Vitebri's algorithm for finding the best matching state sequence to a sequence of symbols.
/**
*	**This function does not include memory allocation!!!**
*
*	Vitebri's algorithm is used for decoding, or finding the best matching state sequence
*	to a sequence of observable symbols. For computation and working configuration it is
*	recommended to use: hmmp_viterbi() or hmmp_decode().
*
*	@param[in] model	The model used for decoding
*	@param[in] seq		The sequence to be decoded
*	@param[in] backtrack Integer backtracking ( by index ) N * T matrix
*	@param[in] mu		Inner probability variables in a 2 * N matrix
*	@param[out] o_state_seq The resulting state sequence with highest probability
*	@param[out] o_logP	The logarithmic probability of the resulting state sequence
						given the observed sequence of symbols
*	@return @ref hmmp_Error Error code.
*/
int hmmp_viterbi_alg( hmmp_Model model,	hmmp_Sequence seq, int_array *backtrack, 
					 dbl_matrix *mu,	hmmp_Sequence *o_state_seq, double *o_logP );

///Part of the Baum-Welch algorithm: Finding the forward-backward variable ( gamma )
/**
*	**This function does not include memory allocation!!!**
*
*	This function is used in the Baum-Welch algorithm. It can also be used in different
*	configurations for learning. Only use this function if you understand how to interpret
*	the result, and be sure to provide valid input as this function is not safe!!!
*	The dimentions of the matrix containing the forward-backward variables should be:\n
*	N x T ( number of states x sequence length ).
*
*	@param[out] o_gamma	The resulting forward-backward variables are saved in this matrix
*	@param[in] alfa		The forward variables.
*	@param[in] beta		The backward variables ( with scaling factors from the forward alg. )
*	@param[in] alfa_scale An array containing the scaling factor for the forward algorithm
*	@param[in] num_states Number of states in the current model
*	@param[in] seq_length The length of the observed sequence
*	@return @ref hmmp_Error Error code.
*/
int hmmp_bwa_gamma_alg ( dbl_matrix *o_gamma, dbl_matrix *alfa, dbl_matrix *beta, 
						dbl_array *alfa_scale, int num_states, int seq_length );

///Part of the Baum-Welch algorithm: Finding the xi variables.
/**
*	**This function does not include memory allocation!!!**
*
*	This function is used in the Baum-Welch algorithm. Only use this function if you 
*	understand how to interpret	the result, and be sure to provide valid input as this 
*	function is not safe!!!
*
*	Finding the xi variables is only possible after obtaining the forward ( alfa )
*	and the backward ( beta ) variables.
*
*	The dimentions of the matrix containing the xi variables should be:\n
*	N x N x (T-1) ( number of states x number of states x (sequence length-1) ).
*
*	@param[out] o_xi	The resulting xi variables are saved in this N*N*T 3d matrix
*	@param[in] alfa		The forward variables.
*	@param[in] beta		The backward variables ( with scaling factors from the forward alg. )
*	@param[in] model	The current working model
*	@param[in] seq		The observed sequence
*	@return @ref hmmp_Error Error code.
*/
int hmmp_bwa_xi_alg ( dbl_matrix *o_xi, dbl_matrix *alfa, dbl_matrix *beta, 
					 hmmp_Model model, hmmp_Sequence seq );

///Part of the Baum-Welch algorithm: Re-estimating the model parameters
/**
*	**This function does not include memory allocation!!!**
*	All output variables must be allocated before using this function.
*
*	This function is used in the Baum-Welch algorithm. Only use this function if you 
*	understand how to interpret	the result, and be sure to provide valid input as this 
*	function is not safe!!!
*
*	Separating the results to four different variables - (nominators and denominators)
*	of the resulting rescaled parameters - (transition and emission) is done to enable
*	accumulation of these results over multiple sequences, following the modified advanced
*	reestimation formulas. Traditionally these formulas does not reestimate the initial
*	parameters as the model in this case provides an initial state with probability 1.
*	The current implementation though provides a reestimate over the initial parameters
*	implemented within the learning hmmp_baum_welch() algorithm.
*
*	**Note** This function does not provide 
*	@param[in] model	The current working model
*	@param[in] seq		The observed sequence
*	@param[in] xi		The previously obtained xi variables
*	@param[in] gamma	The previously obtained forward-backward ( gamma ) variables
*	@param[out] o_a_num	Outputs the nominators of the new transition matrix
*	@param[out] o_b_num Outputs the nominators for the new emission matrix
*	@param[out] o_a_denom Outputs the denominators for the new transition matrix
*	@param[out] o_b_denom Outputs the denominators for the new emission matrix
*	@return @ref hmmp_Error Error code.
*/
int hmmp_bwa_reest_alg ( hmmp_Model model,		hmmp_Sequence seq,
						 dbl_matrix *xi,		dbl_matrix *gamma,
						 dbl_matrix *o_a_num,	dbl_matrix *o_b_num,
						 dbl_array *o_a_denom,	dbl_array *o_b_denom	);
#endif