/**
		F.H.P.C. Assingment 2
		@file cleaned_t_b_one.cc
		@brief thread master fills the array, others makes sum.
		
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
    int *array = (int*)malloc(N * sizeof(int));

  #if defined(_OPENMP)
  if(argc > 1){
     omp_set_num_threads(atoi(*(argv + 1)));  // set the number of threads
  }
  #endif

  for(int ii = 0; ii < N; ++ii){ // thread 0 fills the array
     array[ii] = ii;
  }

  #pragma omp parallel for reduction(+:S) // everyone adding up and reduct on S
  for(int ii = 0; ii < N; ++ii){
     S += array[ii];
  }

   free(array);
   return 0;
}
