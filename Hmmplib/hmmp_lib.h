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
#ifndef HMMP_LIB_H
#define HMMP_LIB_H
/** @file hmmp_lib.h
*	@brief Main header include. This header includes all modules headers.
*
*	This header is created for convinience when incorporating Hmmplib in your project.\n
*	**CAUTION:** Please do not rely on the completeness of this header. Include all 
*	nessecery header files in your project.
*
*/
/** @file hmmp_datatypes.h @brief Contains all datatype definitions: structures, enums.*/
#include "hmmp_datatypes.h" // DONE COMMENTING FOR DOXYGEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/** @file hmmp_alg.h @brief Contains the algorithmic implementations of HMM. */
#include "hmmp_alg.h" // DONE COMMENTING FOR DOXYGEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/** @file hmmp_algwrap.h 
*	@brief Memory allocation and calls to the algorithmic implementations of HMM.*/
#include "hmmp_algwrap.h"// DONE COMMENTING FOR DOXYGEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/** @file hmmp_general.h @brief Contains all the library high-level entry points.*/
#include "hmmp_general.h"

/** @file hmmp_generate.h @brief Generating random or 1 sequence models. Initialization */
#include "hmmp_generate.h" // DONE COMMENTING FOR DOXYGEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/** @file hmmp_file.h @brief Contains all file I/O functions and stucture explanations. */
#include "hmmp_file.h" // DONE COMMENTING FOR DOXYGEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/** @file hmmp_dataproc.h
*	@brief Contains mathematical manipulations of model properties. Copy function. */
#include "hmmp_dataproc.h" // DONE COMMENTING FOR DOXYGEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/** @file hmmp_memop.h
*	@brief Contains memory allocation and deallocation functions for models and sequences.*/
#include "hmmp_memop.h" // DONE COMMENTING FOR DOXYGEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#endif