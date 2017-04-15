\page general General Usage
  \tableofcontents
  
\section started Getting Started
	\subsection creation Creating empty instances and manual initialization from file.
		The simplest way of using Hmmplib is to create instances of the models
		and the sequences. That can be done using the functions declared in hmmp_memop.h.
		 \code{.c}
			#include "hmmp_file.h"
			#include "hmmp_memop.h"
			// ...
			hmmp_Model *m;
			m = hmmp_create_model(10,5); // 10 states 5 symbols
			hmmp_load_real(m->initial,"initial.txt",10);
			hmmp_load_real(m->transition,"transition.txt",10*10);
			hmmp_load_real(m->emission,"emission.txt",10*5);
			m->prior = 0.0;
			m->model_id = 1;
			
			// DO ACTUAL WORK
			
			hmmp_delete_model(m);
		\endcode
		The same can be done for creating a sequence:
		 \code{.c}
			hmmp_Sequence *s;
			s = hmmp_create_sequence(50); // This sequence is 50 symbols long.
			s->seq_id = 1;
			s->length = 50;
			s->cardinality = 10;		// The sequence has 10 distinct possible symbols.
			hmmp_load_int(s->sequence,"sequence1.txt",50);
			
			// DO ACTUAL WORK
			
			hmmp_delete_sequence(s);
		\endcode
		Similarly hmmp_create_arr_models() and hmmp_create_arr_seq() can be used to create arrays of uninitialized models and sequences.
		This however is more advanced and can lead to mistakes. The recommended approach is described in the \ref initialization "Next Subsection".\n
		The file format for simple integer or real data is just numbers separated by ' ' interval or a new line symbol '\\n'.\n
		For example transition.txt:
		\code{.txt}
			0.3 0.6 0.1
			0.2 0.2 0.6
			0.1 0.1 0.8
		\endcode
		For a detailed explaination of the data layout in memory and files check \ref datafilelayout "Data Layout"
	\subsection initialization Generalized creation and initialization from file.
	This is the more general approach from creating and initializing models and sequences from files.\n
	This is a safer and easier way. It also allows for storing multiple sequences or models in a single file.\n
	 \code{.c}
	 	#include "hmmp_file.h"
		#include "hmmp_memop.h"
		// ...
		hmmp_Model *models;
		int num_models;
		num_models = hmmp_load_models(&models, "all_models.txt",50); // 50 - maximum number of models to load
		
		// DO ACTUAL WORK
		
		hmmp_save_models("all_models_modified.txt", models, num_models);
		hmmp_delete_arr_models(models, num_models);
	\endcode
	Here is some example code for creating and initializing a single sequence using the generalized file method:
	\code{.c}
		hmmp_Sequence *single;
		hmmp_load_sequences(&single, "one_sequence.txt",1);
		
		// DO ACTUAL WORK
		
		hmmp_save_sequences("one_sequence_modified.txt", single, 1);
		hmmp_delete_sequence(single);
	\endcode
	The file format for multiple sequences or models can be seen \ref fileformat "Here".\n\n
	**Note:**\n
	`hmmp_load_models()` includes creation of the models.\n
	`hmmp_load_sequences()` includes creation of the sequences.\n
	**After using the models and the sequences to process them with the algorithms it is recommended**\n
	**to save them back in a file and delete the array of sequences/models using the following:**\n
	`hmmp_delete_sequence()`\n
	`hmmp_delete_arr_seq()`\n
	`hmmp_delete_model()`\n
	`hmmp_delete_arr_models()`\n
	
<BR><BR>

------------------------
\section evaluation Evaluation

Definition:
For a given model: λ=( A,B,π ) and a sequence of observable symbols	O<sup>T</sup>=( o<sub>1</sub>,o<sub>2</sub>,o<sub>3</sub>,…,o<sub>T</sub>  ),\n
obtain the probability of the model *λ* to have produced the sequence *O<sup>T</sup>*.\n
In mathematical terms: *P(O<sup>T</sup>│λ)= ?*

The evaluation is implemented using the Forward algorithm: hmmp_forward_alg().

In an active configuration for recognition there are multiple sequences and multiple models.\n
This can be broken down to two distinct cases:\n
- having multiple models and a single sequence.
- having a single model and multiple sequences.

The choice of which case is appropriate for the current problem is made by the user.\n
The examples below exibit how to use the implemented parallelism.  Domain decomposition is implemented for both cases.
	\subsection models Evaluating multiple models against a single sequence.
	Example Usage:\n
	**[The example code does not check for error code function returns. Although it might work with bad input, it is recommended to check function returns.]**
	\code{.c}
	#include "hmmp_lib.h" // Include all library headers.
	#include <time.h>

	int HMMP_NUM_THREADS = 4; //Specifying the number of threads for the task to run on.
	//Optimal number of threads is equal to the number of processor cores.
	// Make sure your system has OpenMP support turned ON!!!

	int evaluate_multiple_models_one_sequence()
	{
		const int NUM_OF_MODELS = 250;
		const int NUM_STATES = 20;
		const int NUM_SYMBOLS = 25;
		const int SEQUENCE_LEN = 1000000;
		
		hmmp_Model *array_models;
		hmmp_Sequence *one_sequence;
		dbl_array *log_probabilities;
		
		//Generating random models and a random sequence.[Can be loaded from file instead.]
		array_models = hmmp_gen_random_models(NUM_OF_MODELS,NUM_STATES,NUM_SYMBOLS, clock());
		one_sequence = hmmp_gen_random_sequences( 1,NUM_SYMBOLS, SEQUENCE_LEN, clock() );

		//Evaluating multiple models against a single sequence
		hmmp_evaluate_models(array_models,NUM_OF_MODELS, *one_sequence, &log_probabilities);
		
		//Saving the result in a file [optional]
		hmmp_save_real("result.txt", log_probabilities, NUM_OF_MODELS, 1);
		
		hmmp_delete_arr_models(array_models, NUM_OF_MODELS);
		hmmp_delete_arr_seq(one_sequence,1);
		hmmp_delete_dbl_array(log_probabilities);
		return 0;
	}
	\endcode
	**NOTE:** The result is stored in a double array and represents the natural logarithm of the corresponding desired probabilies.\n
	Logarithmic scale is used for all probability representations to prevent underflow.
	
	\subsection sequences Evaluating a single model against multiple sequences.
	Example Usage:\n
	**[The example code does not check for error code function returns. Although it might work with bad input, it is recommended to check function returns.]**
	\code{.c}
		#include "hmmp_lib.h" // Include all library headers.
		
		int HMMP_NUM_THREADS = 4; //Specifying the number of threads for the task to run on.
		//Optimal number of threads is equal to the number of processor cores.
		// Make sure your system has OpenMP support turned ON!!!
		
		int evaluate_one_model_multiple_sequences()
		{
			const int NUM_SEQUENCES = 250;
			const int NUM_STATES = 20;
			const int NUM_SYMBOLS = 25;
			//const int SEQUENCE_LEN = The lengths of each sequence are described in the file;
			int num_seq;
			
			hmmp_Model *one_model;
			hmmp_Sequence *arr_sequences;
			dbl_array *log_probabilities;

			//Loading the model and the sequences from files. [Can be randomly generated instead.]
			hmmp_load_models(&one_model, "onemodel.txt", 1);
			num_seq = hmmp_load_sequences(&arr_sequences, "sequences.txt", NUM_SEQUENCES );
			
			//Evaluating one model against multiple sequences
			hmmp_evaluate_sequences(*one_model, arr_sequences, num_seq, &log_probabilities);
				
			//Saving the result in a file [optional]
			hmmp_save_real("result.txt",log_probabilities,NUM_SEQUENCES,1);
			
			hmmp_delete_arr_models(one_model, 1);
			hmmp_delete_arr_seq(arr_sequences,NUM_SEQUENCES);
			hmmp_delete_dbl_array(log_probabilities);
			return 0;
		}
	\endcode
	**NOTE:** The result is stored in a double array and represents the natural logarithm of the corresponding desired probabilies.\n
	Logarithmic scale is used for all probability representations to prevent underflow.
	
<BR><BR>

------------------------
\section decoding Decoding
Definition:
For a given model: λ=( A,B,π ) and a sequence of observable symbols	O<sup>T</sup>=( o<sub>1</sub>,o<sub>2</sub>,o<sub>3</sub>,…,o<sub>T</sub>  ),\n
obtain the sequence of hidden states Q<sup>T</sup> that has the highest probability to have emitted O<sup>T</sup>\n
In mathematical terms: Q<sup>T</sup><sub>max</sub> = argmax<sub>k</sub>P(Q<sup>T</sup><sub>k</sub>|O<sup>T</sup>)

The decoding is implemented using the Viterbi algorithm: hmmp_viterbi_alg().
The current implementation supports multiple sequences and a single model, with domain decomposition parallelism.

Example usage:\n
	**[The example code does not check for error code function returns. Although it might work with bad input, it is recommended to check function returns.]**
	\code{.c}
		#include "hmmp_lib.h" // Include all library headers.
		#include <time.h> // For random generation seed
		
		int HMMP_NUM_THREADS = 4; //Specifying the number of threads for the task to run on.
		//Optimal number of threads is equal to the number of processor cores.
		// Make sure your system has OpenMP support turned ON!!!
		
		int decode_multiple_sequences()
		{
			const int SEQUENCE_length = 1000000;
			const int NUM_OF_MODELS = 1;
			const int NUM_OF_STATES = 25;
			const int NUM_OF_SYMBOLS = 20;
			const int NUM_OF_SEQUENCES = 24;
			hmmp_Model *m;
			hmmp_Sequence *obs_sequences, *state_sequences;
			dbl_array *logP_sequences;

			//Generating random input data. [Can be loaded from file instead.]
			m = hmmp_gen_random_models ( NUM_OF_MODELS,NUM_OF_STATES, NUM_OF_SYMBOLS, clock() );
			obs_sequences = hmmp_gen_random_sequences ( NUM_OF_SEQUENCES , NUM_OF_SYMBOLS, SEQUENCE_length, clock() );
			
			//Decoding the sequences with the model, results in 'state_sequences' and 'logP_sequences'
			hmmp_decode(*m,obs_sequences,NUM_OF_SEQUENCES,&state_sequences, &logP_sequences );

			//Saving the result in a file [optional]
			hmmp_save_sequences("result_hidden_seq.txt", state_sequences, NUM_OF_SEQUENCES);
			hmmp_save_real("correspoding_log_prob.txt", logP_sequences, NUM_OF_SEQUENCES, 1);
			
			hmmp_delete_arr_models ( m , NUM_OF_MODELS );
			hmmp_delete_arr_seq ( obs_sequences, NUM_OF_SEQUENCES );
			hmmp_delete_arr_seq ( state_sequences, NUM_OF_SEQUENCES );
			hmmp_delete_dbl_array( logP_sequences);
		}
	\endcode

For each observable sequence there is an output of a hidden state sequence and the correspoding logarithmic probability.
<BR><BR>

------------------------
\section learning Learning
Definition:
For a given model: λ=( A,B,π ) and a sequence of observable symbols	O<sup>T</sup>=( o<sub>1</sub>,o<sub>2</sub>,o<sub>3</sub>,…,o<sub>T</sub>  ),\n
change the model parameters ( A,B,π ) to maximize the probability of the model to produce this sequence.

In mathematical terms: λ' = *f*( λ ,O<sup>T</sup> )= ? ,when P( λ'│O<sup>T</sup>) = max<sub>k</sub>P(λ<sub>k</sub>│O<sup>T</sup>)

Unfortunately there is no analytical solution to this problem. The most common method for solving the problem is the **Baum-Welch algorithm.**\n
The Baum-Welch algorithm uses EM ( Expectation Maximization ) mathematical technique.\n
After each completed step of learning, the Baum-Welch algorithm promises: P( λ'│O<sup>T</sup>) > P( λ│O<sup>T</sup>) \n
This eventually converges to a local maximum. When training a model it is convinient to set maximum number of steps and/or minimum distance ΔP.\n
In the example below both are used and the learning process stops when one of the conditions is satisfied.

The Baum-Welch algorthm can be broken down to 6 distinct sub-algorithms:
	1. `hmmp_forward_alg()` - Finding forward variables. ( Scaling is applied at each "time" step to prevent underflow. )
	2. `hmmp_backward_alg()` - Finding backward variables.
	3. `hmmp_backward_rescale()` - Rescaling backward variables to use forward scaling factors. This results in easy reduction in the next steps.
	4. `hmmp_bwa_gamma_alg()` - Finding gamma variables.
	5. `hmmp_bwa_xi_alg()` - Finding xi variables.
	6. `hmmp_bwa_reest_alg()` - Using gamma and xi to find the new model parameters.
The current implementation supports multiple sequences and a single model and multiple steps. The parallelism is as follows:\n
- Computation for each sequence is completed linearly ( one by one ). <i>Planning to implement domain decomposition in the future.</i>
- Computation for each step is completed linearly after the previous step. <i>The algorithm doesn't allow parallelism here.</i>

- Steps 1. and 2. Are completed in parallel.
- Step 3. requres step 1. and step.2 and is completed linearly ( lower complexity, no impact on performance )
- Step 4. and 5. require the previous steps and are completed in parallel, containing domain decomposition within them.
- Step 6. requres the previous steps and is completed in parallel with domain decomposition.

**NOTE:** The highest complexity is at step 1. and 2. as well as step 5. Since there is no domain decomposition implemented in step 1. and 2. 
in regard to multiple sequences, the optimal configuration for parralelism is 2 threads for step 1 and step 2,and 'number of processor cores' threads for step 4, step 5 and step 6.

Example usage:\n
	**[The example code does not check for error code function returns. Although it might work with bad input, it is recommended to check function returns.]**
	\code{.c}
	#include "hmmp_lib.h" // Include all library headers.
	#include <time.h> // For random generation seed
	
	int HMMP_NUM_THREADS = 2; //Specifying the number of threads for the task to run on.
	// Make sure your system has OpenMP support turned ON!!!
		
	int learn_example()
	{
		const double DELTA_P = 0.001; // minimum probability distance = log(Pnew/Pold)
		const int MAX_STEPS = 20;
		const int NUM_STATES = 25;
		const int NUM_SYMBOLS = 25;
		const int NUM_SEQUENCES = 5;
		const int SEQ_LENGTH = 100000;
				
		hmmp_Model *model;
		hmmp_Sequence *arr_s;
		int steps;

		//Generating random data. [Can be loaded from file instead.]
		model = hmmp_gen_random_models ( 1,NUM_STATES,NUM_SYMBOLS, clock());
		arr_s = hmmp_gen_random_sequences (NUM_SEQUENCES,NUM_SYMBOLS,SEQ_LENGTH, clock());

		//Execute the Baum-Welch algorithm for learning with multiple repetition steps.		
		steps = hmmp_baum_welch(model,arr_s,NUM_SEQUENCES,MAX_STEPS,DELTA_P);
				
		// Save the model in a file after learning. [Optional]
		hmmp_save_models("manipulated_model.txt", model, 1);
		
		hmmp_delete_arr_models ( model,1 );
		hmmp_delete_arr_seq ( arr_s,NUM_SEQUENCES );
		return 0;
	}
	\endcode
The variable `steps` contains the number of steps completed while learning. \n
It will be either equal to MAX_STEPS or lower if the probability distance function's result is lower than ΔP.

------------------------
\section datafilelayout Data Layout in Memory and Files.
	In Hmmplib all multi-dimentional variables data is presented in a single linear array. For a proper use of the variables,
	parameters and files, an introduction to the chosen memory layout is required.\n\n
	This section explains the data layout of the model matrices ( transition, emission ) and the model variables.\n
	Here you can also find argumentation of why the certain layout is chosen.
	\subsection modelparam Layout of data in the model parameters
		\subsubsection initial Initial parameters
			The model's initial parameters represent the probability of the model to begin emmision in any of the states.\n
			If the model has *N* number of possible states, then there are *N* number of initial probabilities and they are simply stored in an array.
		\subsubsection transition Transition parameters
			The model's transition parameters represent the probabilities of the model to transition from the current state to a new state 
			( can be the same as the previous one ) in exactly one step ( and with exactly one emission after the transition ).\n
			If the model has *N* number of possible states, then <b>there are <i>N*N</i> number of transition probabilities</b>
			( from each state to each other state and to itself ).
			
			Traditionally the transition parameters are presented in an <i>N*N</i> matrix: <b>A = {a<sub>ij</sub>}</b>. \n
			In the program implementation they are presented linearly in an array.
			Each row of the matrix is layed down consequetivly in memory after the previous one. Also known as <b>Row-major layout</b>.\n
			The index of the row **i** represents: transition from a state.\n
			The index of the column **j** within a row represents: transition to a state.
			
			Example:
			\code{.txt}
			Number of states: 3
			Transition probabilities:
			0.05 0.1 0.85 0.15 0.2 0.65 0.4 0.4 0.2
			\endcode
			**Note:** The transition probabilities within a row ( from a state ) follow the standard stochastic constraints. - add up to 1.0 \n
			In this example the first 3 variables represent transition from the first state. The next 3 - from the second state and so on...\n
			0.05 : from state 1 -> state 1\n
			0.85 : from state 1 -> state 3\n
			0.15 : from state 2 -> state 1\n
			0.4  : from state 3 -> state 1 and state 2\n
			0.2  : from state 2 -> state 2 and from state 3 -> state 3\n
			For memory access to a particular transition parameter the following expression can be used:\n
			<b>transition[i*N + j] = probability for transition from state 'i' to state 'j'</b>\n
			The layout for the transition parameters is chosen for convinience, it localizes the access during the backward algorithm and not during the forward algorithm. Since the backward and the forward algorithms have similar usage and both have to be executed in the Baum-Welch algorithm readability and understanding the layout is more important than the layout itself. \n
		\subsubsection emmision Emission parameters
			The model's emission parameters represent the probabilities of the model to emmit any of the possible symbols after making a transition.\n
			**Note:** The emission probabilities also follow the standard stochastic constraints for each state. This means that an emission will happen after a state transition with probability 1.0.\n
			In a discrete implementation of HMM such as Hmmplib a symbol is solely defined by its index.\n
			If the model has *N* possible states and *M* possible symbols, then <b> there are <i>N*M</i> number of emission probabilities.</b>
			( emission from each state of each symbol ).
			
			Traditionally the emission parameters are presented in an <i>N*M</i> matrix: <b>B = {b<sub>jk</sub>}</b>. \n
			Each row of the matrix is layed down consequetivly in memory after the previous one. Also known as <b>Row-major layout</b>.\n
			The index of the row **j** represents: current state.\n
			The index of the column **k** within a row represents: index of one of the possible symbols.
			\code{.txt}
				Number of states: 3
				Number of symbols: 2
				Emission probabilities:
				0.5 0.5 0.3 0.7 0.2 0.8
			\endcode
			0.5 : from state 1 emitting symbol 1 or 2\n
			0.3 : from state 2 emitting symbol 1\n
			0.7 : from state 2 emitting symbol 2\n
			0.2 : from state 3 emitting symbol 1\n
			0.8 : from state 3 emitting symbol 2

			For memory access to a particular emission parameter the following expression can be used:\n
			<b>emission[j*M + k] = emission probability of symbol with index 'k' from state with index 'j'.</b>\n
			**Note:** Emission always occurs after a transition, or after choosing initial state.

	\subsection variablelayout Layout of data in multi-dimentional variables
		\subsubsection forwardvariables Forward variables
			The forward variables are defined as:\n
			<b>α<sub>t</sub>(i) = P( о<sub>1:t</sub>, q<sub>t</sub> = i │ λ )</b>\n
			Where:
			- t - current step ( in a sequence )
			- о<sub>1:t</sub> = о<sub>1</sub>, о<sub>2</sub>, ..., о<sub>t</sub> - observed symbols in the sequence O<sup>T</sup> until the current step t.
			- q<sub>t</sub> = i - index of current state.
			- λ - the model processed against the sequence
			
			In words: <b>α<sub>t</sub>(i)</b> is the probability of the model to be in state <b>i</b> in step <b>t</b> and to have produced the symbols of the observable sequence before that.
			
			Since the forward variables <b>α<sub>t</sub>(i)</b> have two indexes, they are usually placed in a matrix.
			
			Every step 't' the forward variables <b>α<sub>t</sub>(i)</b> are calculated using the forward variables
			from the previous step <b>α<sub>t-1</sub>(i)</b>. One way to optimize the access to these variables is to provide high locality in memory.\n
			This means that the variables in each step 't' have to be adjasent:
			
			α<sub>1</sub>(1), α<sub>1</sub>(2), ..., α<sub>1</sub>(i), - the first step\n
			α<sub>2</sub>(1), α<sub>2</sub>(2), ..., α<sub>2</sub>(i), - the second step\n
			..\n
			α<sub>T</sub>(1), α<sub>T</sub>(2), ..., α<sub>T</sub>(i) - the last step
			
			Traditionally the matrix is placed the opposite way in the mathematical model.\n
			<b>Hmmplib uses a column major layout for the the forward variables in regard to states, both in memory and data files.</b>
			
			The index of the row **t** represents: time steps ( consequetive steps of observation ).\n
			The index of the column **i** the state of the model for the current forward variable.
			
			For memory access to a particular forward variable <b>α<sub>t</sub>(i)</b> the following expression can be used: 
			<b>alfa[t*N + i]</b> ,when N is the number of states.

			**Note:** All forward variables use scaling to 1 at each step. This scaling accumulates since the algorithm uses the variables from the last step as an input.
			In order to obtain the actual values the following formula can be used:
			
			<b>α<sub>t</sub>(i) = ^α<sub>t</sub>(i) / ( c<sub>1</sub> * c<sub>2</sub> * ... * c<sub>t</sub> )</b>\n
			Where:
			α<sub>t</sub>(i) are the unscaled forward variables\n
			^α<sub>t</sub>(i) are the scaled forward variables\n
			c<sub>x</sub> are the scaling values at step x
			
		\subsubsection backwardvariables Backward variables
			The Backward variables are defined as:\n
			<b>β<sub>t</sub>(i) = P( о<sub>(t+1):T</sub> │ q<sub>t</sub> = i, λ )</b>\n

			Where:
			- t - current step ( in a sequence )
			- о<sub>(t+1):T</sub> = о<sub>t+1</sub>, о<sub>t+2</sub> ..., о<sub>T</sub> - observed symbols in the sequence O<sup>T</sup> after the current step t.
			- q<sub>t</sub> = i - index of current state.
			- λ - the model processed against the sequence
			
			In words: <b>β<sub>t</sub>(i)</b> is the probability of the model to emit the rest of the sequence **о<sub>(t+1):T</sub>**, given that the state of the model in time step **t** is **i**.

			Similarly to the @ref forwardvariables, the layout of the Backward variables optimizes access with memory localization.
			The only difference is that the (Backward algorithm)[@ref hmmp_backward_alg] operates backwards through the "time" steps:
			
			β<sub>1</sub>(1), β<sub>1</sub>(2), ..., β<sub>1</sub>(i), - the first step\n
			β<sub>2</sub>(1), β<sub>2</sub>(2), ..., β<sub>2</sub>(i), - the second step\n
			...\n
			β<sub>T</sub>(1), β<sub>T</sub>(2), ..., β<sub>T</sub>(i) - the last step	<- Algorithm starts here\n
			
			Traditionally the matrix is placed the opposite way in the mathematical model.
			Hmmplib uses a column major layout for the the backward variables in regard to states, both in memory and data files.

			The index of the row **t** represents: time steps ( consequetive steps of observation ).\n
			The index of the column **i** the state of the model for the current backward variable.

			For memory access to a particular backward variable <b>β<sub>t</sub>(i)</b> the following expression can be used:
			<b>beta[t*N + i]</b> ,when N is the number of states.
			
			**Note:** All backward variables use scaling to 1 at each step. This scaling accumulates since the algorithm uses the variables from the last step as an input.
			In order to obtain the actual values the following formula can be used:
			
			<b>β<sub>t</sub>(i) = ^β<sub>t</sub>(i) / ( d<sub>t</sub> * d<sub>t+1</sub> * ... * d<sub>T</sub> )</b>\n
			Where:
			β<sub>t</sub>(i) are the unscaled backward variables\n
			^β<sub>t</sub>(i) are the scaled backward variables\n
			d<sub>x</sub> are the scaling values at step x
			
		\subsubsection gammavariable Baum-Welch Gamma variable
			The forward-backward variables or gamma variables are defined as:
			
			**γ<sub>t</sub>(i) = P( q<sub>t</sub> = i | O<sup>T</sup>  , λ )**
			
			Given a sequence and a model **γ<sub>t</sub>(i)** is the probability of the model to be in state **i** at time step **t**.
			With conditional probabilities the equation is transformed to:
			
			**γ<sub>t</sub>(i) = ( α<sub>t</sub>(i) * β<sub>t</sub>(i) ) / P( O<sup>T</sup> │ λ )**
			
			Since the scaled forward and backward variables are used, the denominator is eliminated as it is equal to the product of the scaling factors.
			Only one term remains as it overlaps in both variables, that is **c<sub>t</sub>** - the scaling factor at step **t**. The above formula transforms to:
			
			**γ<sub>t</sub>(i) = ^α<sub>t</sub>(i) * ^β<sub>t</sub>(i) / c<sub>t</sub>**
			
			Where: ^α, ^β are the scaled variables with the same scaling factors.
			
			The gamma variables are a result of computing the forward and the backward algorithms.
			Similarly to the @ref forwardvariables and the @ref backwardvariables, the forward-backward variables or gamma variables
			are layed in memory such that for one time step **t** all the variables are adjasent in regard to states **i**.
			
			γ<sub>1</sub>(1), γ<sub>1</sub>(2), ..., γ<sub>1</sub>(i), - the first step\n
			γ<sub>2</sub>(1), γ<sub>2</sub>(2), ..., γ<sub>2</sub>(i), - the second step\n
			...\n
			γ<sub>T</sub>(1), γ<sub>T</sub>(2), ..., γ<sub>T</sub>(i) - the last step
		
			The index of the row t represents: time steps ( consequetive steps of observation ).
			The index of the column i the state of the model for the current backward variable.

			For memory access to a particular forward-backward variable <b>γ<sub>t</sub>(i)</b> the following expression can be used:
			<b>gamma[t*N + i]</b> ,when N is the number of states.
			
			**Note:** All forward-backward variables are calculated using the same scaling factors for both the forward and backward input variables.
			This means that before using the hmmp_bwa_gamma_alg() rescaling of the backward variables is neccessery using hmmp_backward_rescale().\n

		\subsubsection xivariable Baum-Welch Xi variable
			Xi variable is defined as:
			
			**ξ<sub>t</sub>(i,j) = P( q<sub>t</sub> = i , q<sub>t+1</sub> = j | O<sup>T</sup> , λ )**

			or the probability for the model to be in state **i** at time step **t** and in state **j** in time step **t+1**, while emitting the sequence **O<sup>T</sup>**.
			Using conditional probabilities this equation is transformed to:
			
			**ξ<sub>t</sub>(i,j) = ( α<sub>t</sub>(i) * a<sub>ij</sub> * b<sub>jo(t+1)</sub> * β<sub>t+1</sub>(j) ) / P( O^T │ λ )**
	
			Since the scaled forward and backward variables are used the denominator is eliminated as it is equal to the product of the scaling factors.
			The equation simlifies to:
			
			**ξ<sub>t</sub>(i,j) = α<sub>t</sub>(i) * a<sub>ij</sub> * b<sub>jo(t+1)</sub> * β<sub>t+1</sub>(j)**
			
			The xi variable has three indices and therefore it is a three dimentional variable. In the xi variables are stored linearly in an array.
			The layout of the variables in memory is chosen to acommodate for access and initialization.
			- The index with the smallest memory increment ( adjasent in memroy ) is: **t** - time steps.
			- The index with the medium memory increment ( T number of variables apart ) is: **j** - the state in **t+1**.
			- The index with the largest memory increment ( T*N number of variables apart ) is: **i** - the state in **t**.
			
			Visually it looks like:
			![Caption text](baum_welch_ksi_layout.png)
			
			For memory access to a specific xi variable **ξ<sub>t</sub>(i,j)** the following expression can be used:\n
			**xi[ i*N*T + j*N + t]** ,when **N** is the number of states and **T** is the length of the observed sequence.
			
		\subsubsection viterbivariable Viterbi algorithm variables
			The Vitebri algorithm is used for decoding, or finding a best matching sequence of states to a sequence of observable symbols.\n
			There are three defined variables and data structures used by the algorithm. Because these variables serve a temporary role within the algorithm
			it is not neccessery to provide a detailed description as the algorithm can be used indirectly in hmmp_decode() and hmmp_vitebri() both provide
			abstraction of the implementations and memory allocation for these variables. Here is a brief explaination over each variable:
			
			**μ<sub>t</sub>(q<sub>t</sub>) = max [ P( q<sub>1:t</sub> , o<sub>1:t</sub>) ] over [ q<sub>1:(t-1)</sub> ]** where:
			
			- **μ<sub>t</sub>(q<sub>t</sub>)** - the joined probability of the partial state sequence and partial 
			observed sequence at the most probable state in every step from step 1 to t, given the model λ and the observed sequence O<sup>T</sup>.
			- **q<sub>t</sub>** - the state at time step **t**.
			- **q<sub>1:t</sub>** = q<sub>1</sub>, q<sub>2</sub>, ..., q<sub>t</sub>
			- **o<sub>1:t</sub>** = o<sub>1</sub>, o<sub>2</sub>, ..., o<sub>t</sub>
			
			With conditional probabilities the equation is transformed to:
			
			**μ<sub>t</sub>(q<sub>t</sub>) = ( max [ ⁡μ<sub>t-1</sub>(q<sub>t-1</sub>) * а<sub>q(t-1)i</sub> ] over [ q<sub>t-1</sub> ] )* b<sub>io(t)</sub>**
			
			The algorithm operates over the obvious recursion in the equation. Since the **μ** variables do not contribute to the final result only the
			variables at the current and the previous steps are stored in a 2xN matrix, linearly in memory.
			
			At each step of the algorithm the for each state a backtracking index is stored. This index indicates what was the preceding state which lead
			to a maximum probability at the current state. These indices are stored in a N*T integer matrix layed consequetivly in memory.
			
			At the end of the algorithm a backtracking procedure begins at the state with highest probability in the last step, recursively tracking it's predecessors.
			The result of the backtracking is the desired sequence of states with highest probability of emitting the observed sequence, the probability of which is stored in \n
			**max<sub>qt</sub>[ μ<sub>T</sub>(q<sub>t</sub>) ]**
	\subsection fileformat File format for using the file I/O interface
		Hmmplib uses files in text format with 3 distinct ways to structure the files explained below:
		\subsubsection fileformatdata Raw data format
			Raw data is integers or real numbers stored in either scientific or regular notation. Raw data is printed in high precision 10 base scientific notation.
			Between each number either a new-line '\\n' or a space ' ' is mandatory, this allows for visual structure, but it is ignored in the file I/O functions hmmp_file.h.
			
			Example 1 ( Structured real numbers data in scientific notation ) - Can be used for model parameters or variables.
			\code
				1.00000000000000010E-001 2.00000000000000010E-001 2.99999999999999990E-001 
				4.00000000000000020E-001 5.00000000000000000E-001 5.99999999999999980E-001 
				6.99999999999999960E-001 8.00000000000000040E-001 9.00000000000000020E-001 
			\endcode
			
			Example 2 ( Data in regular notation ) - Can be used for model parameters or variables.
			\code
				0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0
			\endcode
			
			Example 3 ( Integer data ) - Can be used to present observational or state sequences.
			\code
				1 2 3 3 5 2 1 1 4 4 3 11 3 11 4 2 2 1 4 2
			\endcode
		\subsubsection fileformatmodels Muliple models in one file
			Storing multiple models in a file follows a strict file structure.
			
			Example:
			\code {.txt}
				num_models: 3

				model_id: 12
				num_states: 3
				num_symbols: 4
				prior: 0.5
				initial: 1 0 0 
				transition: 0.1 0.1 0.8 0.2 0.2 0.6 0.3 0.3 0.4
				emission: 0.1 0.4 0.4 0.1 0.2 0.2 0.2 0.2 0.4 0.4 0.1 0.1

				model_id: 15
				num_states: 4
				num_symbols: 4
				prior: 0.3
				initial: 0.5 0.4 0 0.1
				transition: 0.1 0.1 0.6 0.2 0.2 0.4 0.3 0.1 0.1 0.1 0.1 0.7 0.2 0.1 0.2 0.5
				emission:  0.1 0.4 0.4 0.1 0.2 0.2 0.2 0.2 0.4 0.4 0.1 0.1 0.2 0.2 0.1 0.5

				model_id: 16
				num_states: 4
				num_symbols: 4
				prior: 0.3
				initial: 4.00000000000000020E-001 5.00000000000000000E-001 0.00000000000000000E+000 1.00000000000000010E-001 
				transition:
				1.00000000000000010E-001 2.00000000000000010E-001 2.99999999999999990E-001 4.00000000000000020E-001
				2.00000000000000010E-001 5.99999999999999980E-001 1.00000000000000010E-001 1.00000000000000010E-001
				2.99999999999999990E-001 2.99999999999999990E-001 2.99999999999999990E-001 1.00000000000000010E-001
				2.50000000000000000E-001 2.50000000000000000E-001 2.50000000000000000E-001 2.50000000000000000E-001 
				emission: 
				2.50000000000000000E-001 2.50000000000000000E-001 1.00000000000000010E-001 4.00000000000000020E-001
				2.00000000000000010E-001 2.00000000000000010E-001 2.00000000000000010E-001 2.00000000000000010E-001
				4.00000000000000020E-001 4.00000000000000020E-001 1.00000000000000010E-001 1.00000000000000010E-001
				1.00000000000000010E-001 2.99999999999999990E-001 2.99999999999999990E-001 1.00000000000000010E-001 
			\endcode
			**Note:**
			- The first parameter is the number of models in the file - 'num_models'
			- There are two "\n\n" newline symbols after and before a model.
			- 'prior' represents the last probability evaluated for the model.
			- The only way to distinguish models is using unique 'model_id'.
			
		\subsubsection fileformatsequences Multiple sequences in one file
			Storing multiple sequences follows a similar structure to storing multiple models.
			
			Example:
			\code
				num_sequences: 3

				seq_id: 1
				length: 11
				cardinality: 31
				sequence: 2 5 1 10 30 30 10 2 5 10 2

				seq_id: 2
				length: 9
				cardinality: 11
				sequence: 2 5 4 4 1 4 1 10 10

				seq_id: 3
				length: 13
				cardinality: 16
				sequence: 2 5 4 4 1 4 1 10 10 15 15 15 15
			\endcode
		
			**Note:**
			- The first parameter is the number of sequences in the file - 'num_sequences'
			- There are two "\n\n" newline symbols after and before a sequences.
			- 'cardinality' represents the number of different symbols available to the sequence.
			- A symbol is defined by it's index.
			- The only way to distinguish models is using unique 'sequence_id'.
			
<BR><BR>

------------------------