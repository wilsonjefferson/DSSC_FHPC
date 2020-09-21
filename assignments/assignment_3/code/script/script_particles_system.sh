#		F.H.P.C. Assingment 3
#		@file script_particles_system.sh
#		@brief script to manage the executables.
		
#		@author Pietro Morichetti
#		@date 20/01/2020
#		@version 1.1 


#!/bin/bash
typeset -i num_mpi_process=${1}
typeset -i num_omp_process=${2}
typeset -i code_version=${3}

# setting the threads number from the environment
export OMP_NUM_THREADS=${num_omp_process} 

echo "------ GENERATE THE PARTICLES SYSTEM ------"
echo ""
mpic++ ./code/exercise_02/generate_mpi.cc -std=c++0x -o ./code/exercise_02/a.out
mpirun -np ${num_mpi_process} ./code/exercise_02/a.out
# the generation is in MPI version
echo ""
echo "------ START THE EVOLUTION OF THE SYSTEM ------"
case ${code_version} in
	0) # starts the MPI version
		echo "                 (MPI version)"
		echo ""
		mpic++ ./code/exercise_02/evolution_mpi.cc -std=c++0x -o ./code/exercise_02/b.out # add -std=c++11 on Ulisse
		mpirun -np ${num_mpi_process} ./codes/exercise_02/b.out ;;
	1) # starts the OpenMP version
		echo "                (OpenMP version)"
		echo ""
		c++ ./code/exercise_02/evolution_omp.cc -std=c++0x -fopenmp ./code/exercise_02/b.out # add -std=c++11 on Ulisse
		./codes/exercise_02/b.out ;;
	2) # starts the hybrid MPI and OpenMP version
		echo "               (Hybrid version)"
		echo ""
		mpic++ ./code/exercise_02/evolution_mpi_omp.cc -std=c++0x -fopenmp ./code/exercise_02/b.out # add -std=c++11 on Ulisse
		mpirun -np ${num_mpi_process} ./code/exercise_02/b.out ;;
	*) # catch not valid input
		echo "${code_version} is not permitted, you have to choose among these comands:"
		echo "- 0: to compile and execute the MPI code version"
		echo "- 1: to compile and execute the OpenMP code version"
		echo "- 2: to compile and execute the Hybrid code version"
		echo "" ;;
esac
 
