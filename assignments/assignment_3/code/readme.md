# MACHINE SPECIFICATIONS

## test machine (local environment)
Acer Aspire E 15, E5-571G-597D, Intel Core i5-4210U, 1.70 GHz
https://www.acer.com/ac/en/LA/content/model/NX.MLZST.001

## effective machine
Ulisse

# COMMENTS ON THE REPORT/CODES

- generation phase (an interpretation of text): each task considers only a portion of the cube particles space,
  in which builds "some" particles.
- evolution phase (an interpretation of the text): each task considers only "some" particles that lives in a 
  portion of the cube particles system assignmented.
- SOCKET_NUM_THREAD_MAX is the maximum number of threads per socket and it is based on the Ulysse number, if      you set a number that is greater than this constant means that you are involved in hyper-threading; so, the     programs will adjust the number of threads equals to SOCKET_NUM_THREAD_MAX, in order to avoid the hyper-    threading.

# MODULE AND SEVERAL VERSIONS

- compilers: c++ (GCC) 4.4.7, mpic++ (from openmpi module)
- module: openmpi/1.8.3/gnu/4.9.2(default)

# EXERCISE 2

## compilation
 - generation of particles system(MPI version)		mpic++ ./codes/exercise_02/generate_mpi.cc -std=c++0x -o ./codes/exercise_02/a.out
 - evolution of the particles system(MPI version) 	mpic++ ./codes/exercise_02/evolution_mpi.cc -std=c++0x -o ./codes/exercise_02/b.out
 - evolution of the particles system(OpenMP version)	c++ ./codes/exercise_02/evolution_omp.cc -std=c++0x -fopenmp ./codes/exercise_02/b.out
 - evolution of the particles system(Hybrid version)  	mpic++ ./codes/exercise_02/evolution_mpi_omp.cc -std=c++0x -fopenmp ./codes/exercise_02/b.out

## execution
 - generation of particles system(MPI version)		mpirun -np ${num_mpi_process} ./codes/exercise_02/a.out
 - evolution of the particles system(MPI version) 	mpirun ./codes/exercise_02/b.out
 - evolution of the particles system(OpenMP version)	./codes/exercise_02/b.out
 - evolution of the particles system(Hybrid version)  	mpirun -np ${num_mpi_process} ./codes/exercise_02/b.out
