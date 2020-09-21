/**
		F.H.P.C Assignment 1
		@file serial_sum.cc
		@brief naive sum method
		
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

int main(){
	 unsigned long long int i, dim,sum = 0;
	 unsigned long long int* array;
	 std::fstream fd;
	 
	 fd.open("./number_operation.txt", std::ios::in);
	fd >> dim;
	fd.close();
	
	if ( dim <=1 ) {
			printf("error: dim minor or equal to 1\n");
			exit(-1) ;
	} 
	
	array = new unsigned long long int[dim];
	for(i=1; i<= dim; i++){
		*(array-1+i) = i;
	}
	for(i=0; i<dim; i++){
		sum += *(array+i);
	 }
	delete [] array;
	// printf("final sum = %lld\n", sum); // if you want to see the final sum
	return 0;
}