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
#ifndef HMMP_FILE_H
#define HMMP_FILE_H
/** @file
*	This file contains all available functions for file I/O operations of the library.
*	The layout format in the files is explained in [File Format](@ref fileformat).
*	There are available functions for: saving/loading - multiple models, multiple sequences,
*	real number data, integer data.
*/
#include "hmmp_datatypes.h"

/// Save multiple models in one file
/** 
*	If file does not exist it will be created. If the file exist it's contents will be
*	discarded!
*
*	@param[in]	filename String containing the file name and relative location
*	@param[in]	arr_models Address of an existing array of models
*	@param[in]	num_models Number of models in the array
*	@return @ref hmmp_Error Error code.
*	@see [File format and layout](@ref fileformatmodels)
*/
int hmmp_save_models ( char *filename, hmmp_Model *arr_models, int num_models );

/// Load multiple models from one file
/** 
*	@param[out]	arr_models Address of uninitialized pointer designated to hold 
*							the address of the resulting array of models
*	@param[in]	filename String containing the name and relative location of the file
*	@param[in]	max_num Maximum number of models to load from the file
*	@return Number of models loaded or @ref hmmp_Error Error code.
*	@see [File format and layout](@ref fileformatmodels)
*/
int hmmp_load_models ( hmmp_Model **arr_models, char *filename, int max_num  );

/// Save multiple sequences in one file
/** 
*	If file does not exist it will be created. If the file exist it's contents will be
*	discarded!
*	@param[in]	filename String containing the file name and relative location
*	@param[in]	arr_seq Address of an existing array of sequences
*	@param[in]	num_seq Number of sequences in the array
*	@return @ref hmmp_Error Error code.
*	@see [File format and layout](@ref fileformatsequences]
*/
int hmmp_save_sequences ( char *filename, hmmp_Sequence *arr_seq, int num_seq );

/// Load multiple sequences from one file
/** 
*	@param[out]	addr_seqp Address of uninitialized pointer designated to hold the address of 
*						the resulting array of sequences
*	@param[in]	filename String containing the name and relative location of the file
*	@param[in]	max_num Maximum number of sequences to load from the file
*	@return Number of sequences loaded or @ref hmmp_Error Error code.
*	@see [File format and layout](@ref fileformatsequences)
*/
int hmmp_load_sequences ( hmmp_Sequence** addr_seqp, char *filename, int max_num );

/// Save real numbered data in 10-base scientific notation to a file.
/** 
*	If file does not exist it will be created. If the file exist it's contents will be
*	discarded!
*
*	@param[in]	filename String containing the file name and relative location
*	@param[in]	data Address of an existing array of real numbers
*	@param[in]	data_count Number of elements/values in the array
*	@param[in] val_per_line Number of elements per line to be printed in the file
*	@return @ref hmmp_Error Error code.
*	@see [File format and layout](@ref fileformatdata)
*/
int hmmp_save_real ( char *filename, double *data, int data_count, int val_per_line );

/// Load real numbered data from a file
/** 
*	@param[in,out]	load_to Address of an existing array that will receive the data
*	@param[in]	filename String containing the file name and relative location
*	@param[in]	max_count Maximum number of values to load ( use the size of the array )
*	@return The number of values loaded from the file or @ref hmmp_Error Error code.
*	@see [File format and layout](@ref fileformatdata)
*/
int hmmp_load_real ( double *load_to, char *filename, int max_count );

/// Save integer data to a file.
/** 
*	If file does not exist it will be created. If the file exist it's contents will be
*	discarded!
*
*	@param[in]	filename String containing the file name and relative location
*	@param[in]	data Address of an existing array of integers
*	@param[in]	data_count Number of elements/values in the array
*	@param[in] val_per_line Number of elements per line to be printed in the file
*	@return @ref hmmp_Error Error code.
*	@see [File format and layout](@ref fileformatdata)
*/
int hmmp_save_int ( char *filename, int *data, int data_count, int val_per_line );

/// Load integer data from a file.
/** 
*	@param[in,out]	load_to Address of an existing array that will receive the data
*	@param[in]	filename String containing the file name and relative location
*	@param[in]	max_count Maximum number of values to load ( use the size of the array )
*	@return The number of values loaded from the file or @ref hmmp_Error Error code.
*
*	@see [File format and layout](@ref fileformatdata)
*/
int hmmp_load_int ( int *load_to, char *filename, int max_count );

#endif