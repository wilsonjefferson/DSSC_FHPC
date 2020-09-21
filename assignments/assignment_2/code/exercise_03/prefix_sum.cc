/**
		F.H.P.C. Assingment 2
		@file prefix_sum.cc
		@brief Implementation of the prefix sum algorithm
		
		@author Pietro Morichetti
		@date 17/12/2019
		@version 1.1 
*/

#include <iostream>
#include <omp.h>

int main(int argc, char**argv){

	const std::size_t  N {1000000}; // sets the size
	double *array = (double*)malloc(N*sizeof(double));
	
	#ifndef _OPENMP
		for(std::size_t j = 0; j < N; ++j){ // fills arrays
			*(array + j) = j + 1;
		//	std::cout << *(array + j) << "  ";
		}
		std::cout << "\n\n";
	#endif
	
	#ifdef _OPENMP
	const int n_thread {(argc > 1) ? atoi(*(argv + 1)) : 4};
	double* partial_prefix_sum = new double[n_thread + 1]{0};
	
	#pragma omp places(cores)
	#pragma omp parallel num_threads(n_thread) // START PARALLEL REGION
	{
		const int register my_thread_id = omp_get_thread_num();
		int register local_partial_prefix {0};
		int register local_offset {0};
		
		#pragma omp for schedule(static)// everyone fills the array
		for(std::size_t j = 0; j < N; ++j){
			*(array + j) = j + 1;
		}
		
		//--------- SERIAL PARALLEL PREFIX SUM ON THE MAIN ARRAY-----------//
		
		#pragma omp for schedule(static)
		for(std::size_t j =  0; j < N; ++j){
			local_partial_prefix += *(array + j); // local_partial_prefix is local variable
			*(array + j) = local_partial_prefix; 
		}
		*(partial_prefix_sum + my_thread_id + 1) = local_partial_prefix; // last element into the prefix array
		
		#pragma omp barrier
		
		//--------- SERIAL PREFIX SUM ON THE PREFIX ARRAY -----------//
		
		for(std::size_t j =  0; j < my_thread_id + 1; ++j){
			local_offset += *(partial_prefix_sum + j); // everyone determines own partial prefix sum value
		}
		
		//--------- SERIAL UPDATE OF THE MAIN ARRAY -----------//
		
		#pragma omp for schedule(static)
		for (int j = 0; j < N; ++j){
			*(array + j) += local_offset;
		}	
	} // END PARALLEL REGION
	
	//--------- PREFIX SUM COMPLETED -----------//
	
/*	for(std::size_t  j = 0; j < N; ++j){ // print prefix sum
		std::cout  << j << "	" << *(array + j) << "\n";
	}
	std::cout << "\n";*/
	
	free(partial_prefix_sum);
	#endif
	
	// ------- SERIAL MODE -------
	
	#ifndef _OPENMP
		long int prefix_sum {*(array)}; // step 0: set first value of prefix sum
		for(std::size_t  i = 1; i < N; i++){ // calc the prefix sum
			prefix_sum += *(array + i);
			*(array + i) = prefix_sum;
		}
	#endif
	
	free(array);
	return 0;
}