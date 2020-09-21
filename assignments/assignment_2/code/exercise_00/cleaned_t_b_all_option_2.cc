/**
		F.H.P.C. Assingment 2
		@file cleaned_t_b_all_option_2.cc
		@brief Fill and sum by means of mpi reduction
		
		@author Pietro Morichetti
		@date 17/12/2019
		@version 1.1 
*/

#include <stdlib.h>
#include <stdio.h>
#include <sched.h>
#include <omp.h>
#include <iostream>

#define N 100

int main(int argc, char **argv){
	
   #if defined(_OPENMP)
   if(argc > 1){
      omp_set_num_threads(atoi(*(argv + 1))); 
   }
   #pragma omp places(cores)
   #endif

   long int S = 0;
   int* array = (int*)malloc(N * sizeof(int));
   
   #pragma omp parallel for reduction(+:S)
   for(int ii = 1; ii <= N; ++ii){
      array[ii] = ii;
      S += array[ii];
   }

   free(array);
   return 0;
}
