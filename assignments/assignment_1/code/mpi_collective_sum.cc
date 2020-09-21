/**
		F.H.P.C Assignment 1
		@file mpi_collective_sum.cc
		@brief The sum by means of mpi collective technique
		
		@author Pietro Morichetti
		@date 08/11/2019
		@version 1.1 
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char *argv[]){
	
	 int myid, numprocs = 0;
	 unsigned long long int i, N, dim, partial_sum, final_sum = 0;
	 unsigned long long int* array;
	 double R, start_comm, end_comm, start_comp, end_comp, t_comm_1 = 0.0;
	 std::fstream fd;
	
     MPI_Init(&argc,&argv);
     MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	 MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	 
	 if(myid == 0){
			
			// read the dimension of the problem from a file
			start_comm = MPI_Wtime();
			fd.open("./number_operation.txt", std::ios::in);
			fd >> dim;
			fd.close();
			end_comm = MPI_Wtime();
			//printf("dim: %ld\n", dim);
			printf("T_read: %10.8f \n", end_comm - start_comm);
			
			if ( dim <=1) {
				printf("error: dim minor or equal to 1\n");
				exit(-1) ;
			} 
	 }
	 start_comm = MPI_Wtime(); // start communication time
	 //root sends to everyone the size of the problem
	 MPI_Bcast(&dim, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	 end_comm = MPI_Wtime();// end communication time
	 t_comm_1 = end_comm - start_comm; // first communication time
	 start_comp = MPI_Wtime(); // computation time
	 N=dim/numprocs; R=dim%numprocs; // determine the size of subproblem (and the rest)
	 //printf("%d: N = %ld, R = %2.2f\n", myid, N, R);
	 if(myid<R){N++;} // process increase the size in order to include the rest
	 array = new unsigned long long int[N]; // everyone allocate the array with its own size
	 for(i=1; i<= N; i++){
				if(myid < R){
					*(array-1+i) =i+N*myid;
				}else{  // process how have to include the rest
					*(array-1+i)=i+N*myid+R; // correct starts ...
				}
				//printf("%d: array[%ld] = %ld\n", myid, i-1, array[i-1]);
	   }
	   partial_sum = 0; // safe operation
	   for(i=0; i<N; i++){
			partial_sum += *(array+i);
	   }
	   
	   //printf("%d: partial_sum = %lld\n", myid, partial_sum);
	   delete [] array;
	   end_comp = MPI_Wtime(); // end computation
	   // root reduces the sums
	   start_comm = MPI_Wtime(); // start communication
	   MPI_Reduce(&partial_sum, &final_sum, 1, MPI_LONG_LONG, 
							MPI_SUM, 0, MPI_COMM_WORLD);
		// NOTE: the MPI_SUM is considered as a part of the communication section,
		// not a computational instruction
							
	   end_comm = MPI_Wtime(); // end communication
		t_comm_1 = t_comm_1 + (end_comm - start_comm); // second communication time
		
	   if(myid==0){printf("final_sum: %lld\n", final_sum);}
	   printf("T_comm %d: %10.8f\n", myid, t_comm_1);
	   printf("T_comp %d: %10.8f\n", myid, end_comp - start_comp);
		
	   MPI_Finalize();
	   return 0;
}