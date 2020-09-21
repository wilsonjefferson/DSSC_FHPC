/**
		F.H.P.C. Assingment 2
		@file cleaned_t_b_all_option_2.cc
		@brief All threads fill the array and perform the sum.
		
		@author Pietro Morichetti
		@date 17/12/2019
		@version 1.1 
*/

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define _GNU_SOURCE
#define N 10000000 // size of the problem

int main(int argc, char **argv){
	
   long int S = 0;
   int* array = (int*)malloc(N * sizeof(int));	
	
   #if defined(_OPENMP)
   if(argc > 1){
      omp_set_num_threads(atoi(*(argv + 1))); // set the number of threads
   }
   #endif

   #pragma omp parallel for // everybody fills the array
   for(int ii = 0; ii < N; ++ii){
      array[ii] = ii;
   }

   #pragma omp parallel for reduction(+:S) // everyone add up and reduce on S
   for(int ii = 0; ii < N; ++ii){
      S += array[ii];
   }

   free(array);
   return 0;
}
