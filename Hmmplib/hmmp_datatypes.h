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
#ifndef HMMP_DATATYPES_H
#define HMMP_DATATYPES_H
/** @file hmmp_datatypes.h 
*	This file contains all the custom data type definitions of Hmmplib. In this file you
*	can find declaration of external constants defined in dataproc.c and error codes enum
*	definition.
*/

/// One dimentional array is used to represent all complex types.
typedef double dbl_matrix;

/// One dimentional array is used to represent all complex types.
typedef double dbl_array;

/// One dimentional array is used to represent all complex types.
typedef int int_array;

/// One dimentional array is used to represent all complex types.
typedef int int_matrix;

/// An instance of this structure represents one hidden markov model.
/** @see hmmp_Model */
struct s_hmmp_Model {
	dbl_matrix *transition;	///< Pointer to model's transition probabilities matrix.
	dbl_matrix *emission;	///< Pointer to model's emission probabilities matrix.
	dbl_array *initial;		///< Pointer to model's initial probabilities vector.
	double prior;	///< The last calculated probability of the model against a sequence.
	int num_states;			///< Number of states in the model.
	int num_symbols;		///< Number of observable symbols of the model.
	int model_id;			///< Model's specific ID used for distinction between models.
};
/// Definition of a model type ommiting the 'struct' keyword.
/** @see s_hmmp_Model */
typedef struct s_hmmp_Model hmmp_Model;
/** @var s_hmmp_Model::transition 
*	Transition matrix mapping as follows:
*	i = transition from state
*	j = transition to state
*	Matrix {Aij} is row-major represented in memory:
*	transition[i*N + j] = probability for transition from state 'i' to state 'j'
 */
/** @var s_hmmp_Model::emission
*	Emission matrix mapping as follows:    
*	j = index of state
*	k = index of symbol to be emmited
*	M = number of symbols.
*	Matrix {Bjk} is row-major represented in memory: 
*	emission[j*M + k] = probability of emitting symbol 'k' from state 'j'
 */

/// An instance of this structure represents one squence
/** This structure can be used for both observable symbol sequences and hidden state
*	transition sequences.*/
/** @see hmmp_Sequence */
struct s_hmmp_Sequence {
	int seq_id;			///< Sequence's specific ID
	int length;			///< Sequence's length ( number of time steps )
	int cardinality;	///< Number of different symbols/states possible.
	int_array *sequence;///< Pointer to an array containing the sequence.
};
/// Definition of a sequence type ommiting the 'struct' keyword.
/** @see s_hmmp_Sequence */
typedef struct s_hmmp_Sequence hmmp_Sequence;

#undef HMMP_DBL_MAX
/** Maximum representable double value used in logarithmic calculations.\n
*	Default value: HMMP_DBL_MAX = 1.7976931348623157e+308 defined in **hmmp_dataproc.c**
*/
extern const double HMMP_DBL_MAX;

#undef HMMP_PRECISION
/**	Minimum distance between two double precision floating point numbers.\n
*	Default value HMMP_PRECISION = 2.2204460492503131e-016 defined in **hmmp_dataproc.c**
*/
extern const double HMMP_PRECISION;	

#undef HMMP_NUM_THREADS
/**	Specify the number of threads for OpenMP multi-threaded solutions in hmmp_general.h\n
*	Use **2 or more** for hmmp_baum_welch()\n
*	Use (Number of available processor cores) otherwise.
*/
extern int HMMP_NUM_THREADS;	
/// Hmmplib error codes.
/**
*	Library error codes. Function ___ to translate to string.
*	@see hmmp_Error
*/
enum _hmmp_lib_error{
	E_SEQUENCE  =  2,	///< Success, but sequence is inconsistent!
	E_FILE_EOF  =  1,	///< Success, and EOF reached!
	E_SUCCESS	=  0,	///< Function successfully completed!
	E_VARIABLES	= -1,	///< HMMvariables is NULL!
	E_INCOMPLETE= -2,	///< Number of States or Symbols is Zero!
	E_PARAMETER	= -3,	///< One or more of the pointers of the parameters is NULL!
	E_FILE_OPEN	= -4,	///< Error opening file!
	E_FILE_READ	= -5,	///< Error reading file!
	E_FILE_FORMAT=-6,	///< Error file format is corrupt!
	E_ALLOCATION= -7,	///< Memory allocation failure!
	E_ARGUMENT	= -8,	///< One or more of the call argument values are not expected!
	E_MEM_OVERFLOW_L1 = -9, ///< Cannot initialize: alfa, beta, gamma, xi variables.
							///< Number of states and/or sequence lenght will cause overflow.
	E_MEM_OVERFLOW_L2 = -10 ///< Cannot initialize: xi variables.
							///< Number of states and/or sequence lenght will cause overflow.
};
/// Definition of error type ommiting the 'enum' keyword.
/** @see _hmmp_lib_error */
typedef enum _hmmp_lib_error hmmp_Error;
#endif
