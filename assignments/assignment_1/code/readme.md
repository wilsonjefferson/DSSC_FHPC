# MACHINE SPECIFICATIONS
## test machine (local environment)
Acer Aspire E 15, E5-571G-597D, Intel Core i5-4210U, 1.70 GHz
https://www.acer.com/ac/en/LA/content/model/NX.MLZST.001

## effective machine
Ulisse

# MODULE AND SEVERAL VERSIONS
- compilers: gcc (GCC) 4.9.2, g++ (GCC) 4.9.2, mpic++ (from openmpi module)
- module: openmpi/1.8.3/gnu/4.9.2(default)
- ulisse version: default

# SECTION 2

## question 2.1
compile the serial version: gcc pi.c -o serial_pi.x
run the serial version: time ./pi.x 1000000

how to load the module: module load openmpi
conpile the parallel version:  mpicc ./mpi_pi.c -o mpi_pi.x
run the parallel version: /usr/bin/time mpirun ./mpi_pi.x 1000000

# SECTION 3

## naive procedure
how to load the module: module load openmpi
compile the mpi_sum: mpic++ ./mpi_sum.cc -o mpi_sum.x
run mpi_sum: /usr/bin/time ./mpi_sum.x

## collective procedure
how to load the module: module load openmpi
compile the mpi_sum: mpic++ ./mpi_collective_sum.cc -o mpi_collective_sum.x
run mpi_sum: /usr/bin/time ./mpi_collective_sum.x


# SECTION 4

I used the ./codes/script_sec_4.sh (code folder) in order to collect all the data from the mpi_sum program.
I used the ./codes/serial_sum.cc (code folder) in order to answer the least question of the section.
