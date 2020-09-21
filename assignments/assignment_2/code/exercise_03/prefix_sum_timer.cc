/**
		F.H.P.C. Assingment 2
		@file prefix_sum_timer.cc
		@brief Implementation of the prefix sum algorithm (timer version)
		
		@author Pietro Morichetti
		@date 17/12/2019
		@version 1.1 
*/

#include <iostream>
#include <omp.h>
#include <ctime>
#include <time.h>

int main(int argc, char**argv){

	const std::size_t  N {1000000000}; // sets the size
	double *array = (double*)malloc(N*sizeof(double));
	
	// ------- clock parameters ---------//
	std::clock_t c_start {0};
	std::clock_t c_end {0};
	double duration {0};
   //-----------------------------------//


	#ifndef _OPENMP	
//                 c_start = std::clock();
		for(std::size_t j = 0; j < N; ++j){ // fills arrays
			*(array + j) = j + 1;
		//	std::cout << *(array + j) << "  ";
		}
//                c_end = std::clock();
//		std::cout << "\n\n";
	#endif
	
	#ifdef _OPENMP
	const int n_thread {(argc > 1) ? atoi(*(argv + 1)) : 4};
    double* partial_prefix_sum = (double*)calloc(n_thread + 1, sizeof(double));
         	
	int counter {n_thread};
	
	#pragma omp places(cores)
	#pragma omp parallel num_threads(n_thread) firstprivate(c_start)// START PARALLEL REGION
	{
		const int register my_thread_id = omp_get_thread_num();
		int register local_partial_prefix {0};
		int register local_offset {0};
		
		
//		c_start = std::clock();
		#pragma omp for schedule(static)// everyone fills the array
		for(std::size_t j = 0; j < N; ++j){
			*(array + j) = j + 1;
		}
		*(partial_prefix_sum + my_thread_id + 1) = local_partial_prefix; // last element into the prefix array
		
		/*#pragma omp critical // clock setting
		{
			if(counter == 1){
				c_end = std::clock();
			}else{
				counter--;
			}
		}*/
		
		#pragma omp barrier
		
		
		c_start = std::clock();
		for(std::size_t j =  0; j < my_thread_id + 1; ++j){
			local_offset += *(partial_prefix_sum + j); // everyone determines its own partial prefix sum value
		}
		
		#pragma omp for schedule(static)
		for (int j = 0; j < N; ++j){
			*(array + j) += local_offset;
		}
		
		#pragma omp critical // clock setting
		{
			if(counter == 1){
				c_end = std::clock();
			}else{
				counter--;
			}
		}
	} // END PARALLEL REGION
	
	free(partial_prefix_sum);
	#endif
	
	// ------- SERIAL MODE -------
	
	#ifndef _OPENMP
		c_start = std::clock();
		long int prefix_sum {*(array)}; // step 0: set first value of prefix sum
		for(std::size_t  i = 1; i < N; i++){ // calc the prefix sum
			prefix_sum += *(array + i);
			*(array + i) = prefix_sum;
		}
		c_end = std::clock();
	#endif

        //-------- clock operation -----------//        
        duration = (c_end - c_start ) / (double) CLOCKS_PER_SEC;
        std::cout << (double)N*2*sizeof(double) / (1024*1024*1024) / duration << "\n";
        //-------------------------------------//
        	
	free(array);
	return 0;
}
