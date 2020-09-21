/**
		F.H.P.C. Assingment 1
		@file mpi_sum.cc
		@brief The sum by means of mpi technique
		
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
	 unsigned long long int i, dim, N, partial_sum = 0;
	 unsigned long long int* array;
	 double R, start_comm, end_comm, start_comp, end_comp, t_comm_1 = 0.0;
	 MPI_Status status;
	 std::fstream fd;

     MPI_Init(&argc,&argv);
     MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	 MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	 
	 if(myid == 0){
			
			// read the dimension of the problem from a file
			start_comm = MPI_Wtime(); // start read time (master)
			fd.open("./number_operation.txt", std::ios::in);
			fd >> dim;
			fd.close();
			end_comm = MPI_Wtime(); // end read master (master)
			//printf("dim: %ld\n", dim);
			printf("T_read: %10.8f \n", end_comm - start_comm);
				
			if ( dim <=1 ) {
				printf("error: dim minor or equal to 1\n");
				exit(-1) ;
			} 
			
			start_comm = MPI_Wtime(); // start first communication time (master)
			for(i=1; i<numprocs; i++){ // root sends to everyone the size of the problem
				MPI_Ssend(&dim, 1, MPI_LONG_LONG, i, 10, MPI_COMM_WORLD);
				//printf("%d: sended to %d\n", myid, i);
			}
	 }else{
			// processes receive from the root...
			start_comm = MPI_Wtime(); // start first communication time (slaves)
			MPI_Recv(&dim, 1, MPI_LONG_LONG, 0, 10, MPI_COMM_WORLD, &status);
			//printf("%d: recieved from 0: %ld\n", myid, dim);
	 }
	  end_comm = MPI_Wtime();
	  t_comm_1 = end_comm - start_comm;
	  start_comp = MPI_Wtime(); // start computation time
	  N=dim/numprocs; R=dim%numprocs; // determine the size of subproblem (and the rest)
	  //printf("%d: N = %ld, R = %2.2f\n", myid, N, R);
	  
	  if(myid < R){N++;} // the process increases the size to include the rest
	  array = new unsigned long long int[N]; // everyone allocate the array with its own size
	  
	  for(i=1; i <= N; i++){
	  for(i=1; i <= N; i++){
			if(myid < R){
				*(array-1+i) =i+N*myid;
			}else{  // process how have to include the rest
				*(array-1+i)=i+N*myid+R; // correct starts ...
			}
			//printf("%d: array[%ld] = %ld\n", myid, i-1, array[i-1]);
	   }
	  
	  for(i=0; i<N; i++){
		partial_sum += *(array+i);
	  }
	   
	   //printf("%d: partial_sum = %lld\n", myid, partial_sum);
	   delete [] array;
	   
		if(myid==0){
			array = new unsigned long long int[numprocs-1]; // buffer of the sums
			start_comm = MPI_Wtime(); // start second communication time
			for(i=1; i<numprocs; i++){
					MPI_Recv((array+i-1), 1, MPI_LONG_LONG, i, 10, MPI_COMM_WORLD,&status) ;
					//printf("received from %ld: %ld\n", i, array[i-1]);
			}
			end_comm = MPI_Wtime(); // end second communication time
			t_comm_1 = t_comm_1 + (end_comm - start_comm); // all the communication time
			for(i=0; i<numprocs-1; i++){
					partial_sum += *(array+i);
			}
			delete [] array;
			end_comp = MPI_Wtime(); // end computation
			printf("final sum: %lld\n", partial_sum);
		}else{
			end_comp = MPI_Wtime(); // end computation time
			start_comm = MPI_Wtime(); // start second communication time
			MPI_Ssend(&partial_sum, 1 , MPI_LONG_LONG, 0, 10 , MPI_COMM_WORLD);
			//printf("sended to 0: %ld\n", partial_sum);
			end_comm = MPI_Wtime(); // end second communication time
			t_comm_1 = t_comm_1 + (end_comm - start_comm); // all the communication time
		}
		
		printf("T_comm %d: %10.8f\n", myid, t_comm_1);
		
		if(myid==0){printf("T_comp %d: %10.8f\n", myid, end_comp - start_comp - t_comm_1);}
		printf("T_comp %d: %10.8f\n", myid, end_comp - start_comp);
		
		MPI_Finalize();
		return 0;
}