/**
		F.H.P.C. Assingment 2
		@file cleaned_t_b_all_option_1.cc
		@brief Parallel fill of an array and serial prefix sum
		
		@author Pietro Morichetti
		@date 17/12/2019
		@version 1.1 
*/

#include <iostream>
#include <omp.h>

int main(int argc, char**argv){
	
	constexpr std::size_t  N {10000}; // set the size
	long int S {0}; // the counter
	double *array = (double*)malloc(N*sizeof(double));
	
	int register n_thread {omp_get_thread_num()}; // number of threads taken from the environment
	std::size_t register dim_sub_array {N/n_thread}; // size of each portion of the array
	long long int R {(N%n_thread != 0) ? N%n_thread : 0}; // verify if exist some remainder
	
	#pragma omp parallel num_threads(n_thread) reduction(+:S) // start parallel region
	{
		const int register my_thread_id = omp_get_thread_num();
		
		// ------- PARALLEL FILLS OF THE ARRAY -------
		
		if(my_thread_id < R || R == 0){ // threads with TID less than R have to manage it (their portion of array is increased by one), otherwise if no remains
			if(R != 0){++dim_sub_array;} //  increase the portion size by one
			for(std::size_t j =  0; j < dim_sub_array; ++j){
				*(array + my_thread_id*dim_sub_array + j) = j;
			}
		}else{
			for(std::size_t j = 0; j < dim_sub_array; ++j){
				*(array + my_thread_id*dim_sub_array + j + R) = j;
			}
		}
		
		// ------- PERFORMS SERIAL PARALLEL PREFIX SUM ON EACH PORTION OF THE ARRAY -------
		
		if(my_thread_id < R || R == 0){ // threads with TID less than R have to manage it (their portion of array is increased by one), otherwise if no remains
			for(std::size_t j =  1; j <= dim_sub_array; ++j){
				S += *(array + my_thread_id*dim_sub_array + j - 1);
			}
		}else{
			for(std::size_t j = 1; j <= dim_sub_array; ++j){
				S += *(array + my_thread_id*dim_sub_array + j + R - 1);
			}
		}
		
	} // end parallel region
	
	free(array);
	return 0;
}
