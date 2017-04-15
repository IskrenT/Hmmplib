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
#ifndef HMMP_GENERATE
#define HMMP_GENERATE

#include "hmmp_datatypes.h"

/// Generate a model perfectly matching to a sequence.
/**	
*	The resulting model has probability 1.0 of emitting the sequence
*
*	**Note:** Each time step is matched with a state in the model. Therefore it is not
*	practical to use this function for sequences with a significant size.
*	This function is designed for testing purposes only.
*
*	@param[in] seq The sequence that the resulting model will match.
*	@return The address of the newly created model, matching the sequence.
*			Returns 0 ( NULL ) if the model could not be created.
*/
hmmp_Model *hmmp_gen_from_seq ( hmmp_Sequence *seq );

///	Generate an array of models with random parameters and the same state/symbol space.
/**
*	@param[in] count		Number of models to create
*	@param[in] num_states	Number of states for every model
*	@param[in] num_symbols	Number of symbols for every model
*	@param[in] seed			Integer seed for the random function
*	@return Returns the address of the newly created array of models.\n
*			Returns 0 ( NULL ) if the array could not be created.
*/
hmmp_Model *hmmp_gen_random_models ( int count, int num_states, int num_symbols, int seed );

///	Generate an array of random sequences with the same length and symbol space.
/**
*	@param[in] count		Number of sequences to create
*	@param[in] cardinality	Number of different possible symbols ( symbol space )
*	@param[in] length		The length of the desired sequence
*	@param[in] seed			Integer seed for the random function
*	@return Returns the address of the newly created array of sequences.\n
*			Returns 0 ( NULL ) if the array could not be created.
*/
hmmp_Sequence *hmmp_gen_random_sequences ( int count, int cardinality, int length, int seed );

#endif