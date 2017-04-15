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
#include "hmmp_file.h"
#include "hmmp_datatypes.h"
#include <stdio.h>
#include <malloc.h>
int hmmp_save_models ( char *filename, hmmp_Model *arr_models, int num_models )
{
	FILE *file;
	int i, j, n_times_n, n_times_m;
	char flag_broken = 0;
	if ( !arr_models || !filename )
		return E_PARAMETER;
	if(!(file = fopen ( filename , "w" )))
		return E_FILE_OPEN;
	fprintf( file, "num_models: %d\n\n", num_models);
	for ( i = 0 ; i < num_models ; ++i ){
		fprintf( file, "model_id: %d\nnum_states: %d\nnum_symbols: %d\nprior: %.*E\n",
				arr_models[i].model_id,		arr_models[i].num_states, 
				arr_models[i].num_symbols, 17,arr_models[i].prior );
		fprintf( file, "initial: ");
		for ( j = 0 ; j < arr_models[i].num_states ; ++j )
			fprintf(file, "%.*E ", 17, arr_models[i].initial[j] );

		fprintf( file, "\ntransition: ");
		n_times_n = arr_models[i].num_states * arr_models[i].num_states; 
		for ( j = 0 ; j < n_times_n ; ++j )
			fprintf(file, "%.17E ", arr_models[i].transition[j] );

		fprintf( file, "\nemission: ");
		n_times_m = arr_models[i].num_states * arr_models[i].num_symbols;
		for ( j = 0 ; j < n_times_m ; ++j )
			fprintf(file, "%.*E ", 17, arr_models[i].emission[j] );
		fprintf(file, "\n\n");
	}
	fclose(file);
	return E_SUCCESS;
}

int hmmp_load_models ( hmmp_Model **arr_models, char *filename, int max_num  )
{
	FILE *file;
	hmmp_Model *p_models;
	int i, j, num_models, model_id, num_states, num_symbols, n_times_n, n_times_m;
	double prior;
	char flag_broken = 0;

	if ( !arr_models || !filename )
		return E_PARAMETER;
	if(!(file = fopen ( filename , "r" )))
		return E_FILE_OPEN;

	if ( fscanf(file,"num_models: %d\n", &num_models ) != 1 )
		return E_FILE_READ;

	if ( max_num && max_num < num_models )
		num_models = max_num;
	p_models = ( hmmp_Model *)malloc(num_models * sizeof ( hmmp_Model ));
	if ( !p_models )
		return E_ALLOCATION;

	for ( i = 0 ; i < num_models ; ++i ){
		if ( fscanf(file,"model_id:%d\nnum_states:%d\nnum_symbols:%d\nprior:%lf\n", 
			&model_id, &num_states, &num_symbols, &prior ) != 4 ){
			flag_broken = 1;
			break;
		}
		p_models[i].initial = (dbl_matrix *)malloc( sizeof(double)*num_states );
		if ( !p_models[i].initial ){
			flag_broken = 2;
			break;
		}
		p_models[i].transition = (dbl_matrix *)malloc( sizeof(double)*num_states*num_states );
		if ( !p_models[i].transition ){
			free(p_models[i].initial);
			flag_broken = 2;
			break;
		}
		p_models[i].emission = (dbl_matrix *)malloc( sizeof(double)*num_states*num_symbols );
		if ( !p_models[i].emission ){
			free(p_models[i].initial);
			free(p_models[i].transition);
			flag_broken = 2;
			break;
		}		
		p_models[i].model_id = model_id;
		p_models[i].num_states = num_states;
		p_models[i].num_symbols = num_symbols;
		p_models[i].prior = prior;

		if ( fscanf(file,"initial:") == EOF ){
			flag_broken = -1;
			break;
		}
		if ( flag_broken == -1 )
			break;
		for ( j = 0 ; j < num_states ; ++j ){
			if (fscanf(file, "%lf", &p_models[i].initial[j] ) != 1 ){
				flag_broken = -1;
				break;
			}
		}
		if ( flag_broken == -1 )
			break;
		if ( fscanf(file,"\ntransition:") == EOF ){
			flag_broken = -1;
			break;
		}
		n_times_n = num_states * num_states;
		for ( j = 0 ; j < n_times_n ; ++j ){
			if (fscanf(file, "%lf", &p_models[i].transition[j] ) != 1 ){
				flag_broken = -1;
				break;
			}
		}
		if ( flag_broken == -1 )
			break;
		if ( fscanf(file,"\nemission:") == EOF ){
			flag_broken = -1;
			break;
		}
		n_times_m = num_symbols * num_states;
		for ( j = 0 ; j < n_times_m ; ++j ){
			if (fscanf(file, "%lf", &p_models[i].emission[j] ) != 1 ){
				flag_broken = -1;
				break;
			}
		}
		if ( flag_broken == -1 )
			break;
		if ( fscanf(file,"\n\n") == EOF ){
			flag_broken = -1;
			break;
		}
	}
	fclose(file);
	if ( flag_broken == 1 || flag_broken == 2 || flag_broken == -1){
		num_models = i;
		if ( flag_broken == -1 )
			++num_models;
		for ( i = 0 ; i < num_models ; ++i ){
			free ( p_models[i].initial );
			free ( p_models[i].transition );
			free ( p_models[i].emission );
			p_models[i].model_id = -1;
		}
		free ( p_models );
		if ( flag_broken == 2 )
			return E_ALLOCATION;
		return E_FILE_FORMAT;
	}
	*arr_models = p_models;
	return num_models;
}


int hmmp_save_sequences ( char *filename, hmmp_Sequence *arr_seq, int num_seq )
{
	FILE *file;
	int i, j;
	char flag_cardinality = 0;
	if(!arr_seq || !filename)
		return E_PARAMETER;
	if(!(file = fopen ( filename , "w" )))
		return E_FILE_OPEN;
	fprintf( file, "num_sequences: %d\n\n", num_seq);
	for ( i = 0 ; i < num_seq ; ++i ){
		fprintf( file, "seq_id: %d\nlength: %d\ncardinality: %d\n",
				arr_seq[i].seq_id, arr_seq[i].length, arr_seq[i].cardinality );
		fprintf( file, "sequence: ");
		for ( j = 0 ; j < arr_seq[i].length ; ++j ){
			if(arr_seq[i].sequence[j] >= arr_seq[i].cardinality)
				flag_cardinality = 1;
			fprintf(file, "%d ", arr_seq[i].sequence[j] );
		}
		fprintf(file, "\n\n");
	}
	fclose(file);
	if(flag_cardinality)
		return E_SEQUENCE;
	return E_SUCCESS;
}

int hmmp_load_sequences ( hmmp_Sequence **addr_seqp, char *filename, int max_num )
{
	FILE *file;
	hmmp_Sequence *sequences;
	int i, j, num_seq, seq_id, length, cardinalty;
	char flag_broken = 0;
	
	if(!addr_seqp || !filename )
		return E_PARAMETER;

	if(!(file = fopen ( filename , "r" )))
		return E_FILE_OPEN;

	if ( fscanf(file,"num_sequences: %d\n", &num_seq ) != 1 )
		return E_FILE_READ;

	if ( max_num && max_num < num_seq )
		num_seq = max_num;
	sequences = ( hmmp_Sequence *)malloc(num_seq * sizeof ( hmmp_Sequence ));
	if ( !sequences )
		return E_ALLOCATION;

	for ( i = 0 ; i < num_seq ; ++i ){
		if ( fscanf(file,"seq_id:%d\nlength:%d\ncardinality:%d\n", 
			&seq_id, &length, &cardinalty ) != 3 ){
			flag_broken = 1;
			break;
		}
		sequences[i].sequence = (int_array *)malloc( sizeof(int_array)*length );
		if ( !sequences[i].sequence){
			flag_broken = 2;
			break;
		}
		sequences[i].seq_id = seq_id;
		sequences[i].length = length;
		sequences[i].cardinality = cardinalty;
		if ( fscanf(file,"sequence:") == EOF ){
			flag_broken = -1;
			break;
		}
		for ( j = 0 ; j < length ; ++j ){
			if (fscanf(file, "%d", &sequences[i].sequence[j] ) != 1 ){
				flag_broken = -1;
				break;
			}
			if ( sequences[i].sequence[j] >= cardinalty ){
				sequences[i].cardinality = cardinalty = sequences[i].sequence[j];
			}
		}
		if ( flag_broken == -1 )
			break;
		if ( fscanf(file,"\n\n") == EOF ){
			flag_broken = -1;
			break;
		}
	}
	fclose(file);
	if ( flag_broken == 1 || flag_broken == 2 || flag_broken == -1){
		num_seq = i;
		if ( flag_broken == -1 )
			++num_seq;
		for ( i = 0 ; i < num_seq ; ++i ){
			free ( sequences[i].sequence );
			sequences[i].seq_id = -1;
		}
		free ( sequences );
		if ( flag_broken == 2 )
			return E_ALLOCATION;
		return E_FILE_FORMAT;
	}
	*addr_seqp = sequences;
	return num_seq;
}


int hmmp_save_real ( char *filename, double *data, int data_count, int val_per_line )
{
	FILE *file1;
	int i;
	if(!data || !filename)
		return E_PARAMETER;
	if(!(file1 = fopen ( filename , "w" )))
		return E_FILE_OPEN;

	for ( i = 0 ; i < data_count ; ++i ){
		fprintf(file1, "%.*E ", 17, data[i] );
		if(!((i+1) % val_per_line ))
			fprintf(file1, "\n");
	}
	fclose(file1);
	return E_SUCCESS;
}

int hmmp_load_real ( double *load_to, char *filename, int max_count )
{
	FILE *file1;
	int i;
	if(!load_to || !filename)
		return E_PARAMETER;
	if(!(file1 = fopen ( filename , "r" )))
		return E_FILE_OPEN;
	for ( i = 0 ; i < max_count ; ++i ){
		if (fscanf(file1, "%lf", &load_to[i] ) == EOF )
			break;
	}
	fclose(file1);
	return i;
}

int hmmp_save_int ( char *filename, int *data, int data_count, int val_per_line )
{
	FILE *file1;
	int i;
	if(!data || !filename )
		return E_PARAMETER;
	if(!(file1 = fopen ( filename , "w" )))
		return E_FILE_OPEN;

	for ( i = 0 ; i < data_count ; ++i ){
		fprintf(file1, "%d ", data[i] );
		if(!((i+1) % val_per_line ))
			fprintf(file1, "\n");
	}
	fclose(file1);
	return E_SUCCESS;
}

int hmmp_load_int ( int *load_to, char *filename, int max_count )
{
	FILE *file1;
	int i;
	if(!load_to || !filename)
		return E_PARAMETER;
	if(!(file1 = fopen ( filename , "r" )))
		return E_FILE_OPEN;
	for ( i = 0 ; i < max_count ; ++i ){
		if (fscanf(file1, "%d", &load_to[i] ) == EOF )
			break;
	}
	fclose(file1);
	return i;
}