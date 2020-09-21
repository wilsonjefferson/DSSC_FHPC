/**
		F.H.P.C. Assingment 2
		@file cleaned_t_b_one_timer.cc
		@brief Fill and sum by means of mpi reduction (timer version)
		
		@author Pietro Morichetti
		@date 17/12/2019
		@version 1.1 
*/

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <ctime>
#define N 1000000000 // size of the problem

int main(int argc, char **argv){

  long int S = 0;
  int *array = (int*)malloc(N * sizeof(int));

  // ------- clock parameters ---------//
  std::clock_t c_start {0};
  std::clock_t c_end {0};
  double duration {0};
  //----------------------------------//
  
  #if defined(_OPENMP)
  if(argc > 1){
     omp_set_num_threads(atoi(*(argv + 1)));  // set the number of threads
  }
  #endif

//  c_start = std::clock();
  for(int ii = 1; ii <= N; ++ii){ // thread 0 fills the array
     array[ii] = ii;
  }
//  c_end = std::clock();
  
  c_start = std::clock();
  #pragma omp parallel for reduction(+:S) // everyone adding up and reduct on S
  for(int ii = 1; ii <= N; ++ii){
     S += array[ii];
  }
  c_end = std::clock();
  
  // --------- clock operation ----------//
  duration = (c_end - c_start ) / (double) CLOCKS_PER_SEC;
  printf("%6.3g\n", (double)N*2*sizeof(double) / (1024*1024*1024) /duration);
  //----------------------------------//

   free(array);
   return 0;
}
