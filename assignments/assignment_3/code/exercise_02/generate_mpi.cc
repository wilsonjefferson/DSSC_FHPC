/**
		F.H.P.C. Assingment 3
		@file generate_mpi.cc
		@brief mpi processes generate the particles system
		
		@author Pietro Morichetti
		@date 20/01/2020
		@version 1.1 
*/

#include <iostream>
#include <mpi.h>
#include <random>
#include <stdio.h>
#include <math.h>

#define REF_PROCESS 0 // shows the operations performed by one specific mpi process
#define SUMMARY 1 // shows the results information of the particles system ( general information, print of array, ....)
#define CHECK 0 // shows all support informations (malloc,  ....)

typedef struct {double coord[3];} vect;
typedef struct {vect pos; vect vel; double E;} particle;
typedef struct {int info[3]{8, 0, 0}; double evolution_time[1]{0};} metanode;

bool operator==(const vect& p1, const vect& p2){
	return p1.coord[0] == p2.coord[0] && p1.coord[1] == p2.coord[1] && p1.coord[2] == p2.coord[2];
}

void build_mpi_data_type(MPI_Datatype& MPI_particle_type, MPI_Datatype& MPI_metanode_type){
	
	// build MPI_particle_type
	const std::size_t  num_members_1 {3};
	int blocklengths_1[num_members_1] = {3, 3, 1};
	
	MPI_Aint displacements_1[num_members_1] = {offsetof(particle, pos), offsetof(particle, vel), offsetof(particle, E)};
	MPI_Datatype types_1[num_members_1] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	
	MPI_Type_create_struct(num_members_1, blocklengths_1, displacements_1, types_1, &MPI_particle_type);
	MPI_Type_commit(&MPI_particle_type);
	
	//----------
	
	// build MPI_inode_type
	const std::size_t  num_members_2 {2};
	int blocklengths_2[num_members_2] = {3, 1};
	
	MPI_Aint displacements_2[num_members_2] = {offsetof(metanode, info), offsetof(metanode, evolution_time)};
	MPI_Datatype types_2[num_members_2] = {MPI_INTEGER, MPI_DOUBLE};
	
	MPI_Type_create_struct(num_members_2, blocklengths_2, displacements_2, types_2, &MPI_metanode_type);
	MPI_Type_commit(&MPI_metanode_type);
}

double norm(const vect& p){
	return pow(p.coord[0]*p.coord[0] + p.coord[1]*p.coord[1] + p.coord[2]*p.coord[2], 0.5);
}

int main(int argc, char *argv[]){
	
	const int N_p {20}; // number of particles
	int myid, N_t {0}; // N_t = number of tasks i.e. processors
	
	//------- DEFINE MPI SETTING ------//
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&N_t);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Datatype MPI_particle_type, MPI_metanode_type;
	MPI_File fp;
	//----------------------------------//
	
	if(N_t > N_p){ // check the number of "workers"
		if(SUMMARY) if(myid == REF_PROCESS) 
			std::cout << "------------------------------------------------------------\n"
								 << "------------------------------------------------------------\n\n"
								 << "ERROR: the dimension of the problem is less than the number of processes, \n"
								 << "it means some of them are unsed.\n" 
								 << "Set an appropriate number of \"workers\" for this problem.\n\n"
								 << "------------------------------------------------------------\n"
								 << "------------------------------------------------------------\n";
		MPI_Finalize();
		return 1;		
	}
	
	// mass of the system
	const double M_max {100};
	double m {M_max/N_p};
	
	// number of particles managed by one single processor and the potential reminder
	int R {N_p%N_t};
	int num_part_per_task {N_p/N_t}; 
	if(myid < R){++num_part_per_task;} 
	
	particle * array_particle = (particle*)malloc(num_part_per_task*sizeof(particle));
	if(array_particle == NULL){
		if(SUMMARY) if(myid == REF_PROCESS) std::cout << "ERROR: the process " << myid << "is not able to use the \"malloc\" function.\n";
		MPI_Finalize();
		return 1;
	}
	
	// setting the start read/write in the file, for each process
	MPI_Offset offset {(myid < R) ? 20 + myid*num_part_per_task*((int)sizeof(particle)) : 
														 20 + (myid*num_part_per_task + R)*(int)sizeof(particle)};

	if(SUMMARY) if(myid == REF_PROCESS) std::cout  << "The speaker is " << myid << "\n"
																					<< "N_p = " << N_p << ", N_t = " << N_t << "\n"
																					<< "(my) num_part_per_task = " << num_part_per_task 
																					<< ", R = " << R << ", (my) offeset = " << offset << "\n";
	
	// information about the system
	metanode meta_data {}; // first block of bytes that contain general information about the data_ic file
	meta_data.info[1] = N_p; // set this specific block of bytes as the number of particles
	
	//------- DEFINE NEW DATATYPE ------//
	build_mpi_data_type(MPI_particle_type, MPI_metanode_type);
	
	//------- UNIFORM AND NORMAL DISTRIBUTION PARAMETERS ------//
	double start_range {(double)myid/N_t}; // define the start point of the sheet space system
	double offset_range {(double)1/N_t}; // define the end point of the sheet space system
	
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> uniform_setted_interval(start_range, start_range + offset_range); // defines the coordinates range
    std::uniform_real_distribution<double> uniform(0, 1); // defines the coordinates range
	std::normal_distribution<double> normal(0.0, 0.01); // defines the velocity distribution
	//-------------------------------------------------------------//
	
	//------- STARTS TO BUILD PARTICLES ------//
	for(int i = 0; i < num_part_per_task; ++i){ // each process builds our particles
		array_particle[i].pos.coord[0] = uniform_setted_interval(gen);
		array_particle[i].pos.coord[1] = uniform(gen);
		array_particle[i].pos.coord[2] = uniform(gen);
		
		for(int k = 0; k < i; ++k){ // check if just exist
			if(array_particle[i].pos == array_particle[k].pos){
					if(SUMMARY)if(myid == REF_PROCESS) std::cout << "invalid particle\n";
					--i; 
					break;
			}
		}
		
		array_particle[i].vel.coord[0] = normal(gen);
		array_particle[i].vel.coord[1] = normal(gen);
		array_particle[i].vel.coord[2] = normal(gen);
		
		array_particle[i].E = 0.5*m*norm(array_particle[i].vel)*norm(array_particle[i].vel);						 
	}
	//------------------------------------------//
	
	// open file
	if(MPI_File_open(MPI_COMM_WORLD, "data_ic", MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &fp) == 1){
		if(SUMMARY) if(myid == REF_PROCESS) std::cout << "ERROR: the process " << myid << "is not able to open the \"data_ic\" file.\n";
		MPI_Finalize();
		return 1;
	}
	
	// process zero writes the meta_data of the file
	if(myid == 0) MPI_File_write(fp, &meta_data, 1, MPI_metanode_type, &status);
	
	// each process writes our particles starting from specific position in the file
	MPI_File_write_at_all(fp, offset, array_particle, num_part_per_task, MPI_particle_type, &status);
	
	if(CHECK){ // check if the particles are corretcly written in the file
		//------- CHECKING WRITTEN INFORMATION ------//
		delete[] array_particle;
		array_particle = (particle*)calloc(N_p, sizeof(particle));
		
		MPI_File_read_at_all (fp, 20, array_particle, N_p, MPI_particle_type, &status);
		
		if(myid == REF_PROCESS){
			std::cout << "------- print array of positions -------\n";
			for(int i = 0; i < N_p; ++i){
				std::cout << i << ": ";
				for(int j = 0; j < 3 ; ++j){
					std::cout << array_particle[i].pos.coord[j] << " ";
				}
				std::cout << "\n";
			}
			std::cout << "----------------------------------------\n";
		}
		//----------------------------------------------//
	}
	
	delete[] array_particle;
	MPI_Type_free(&MPI_particle_type);
	MPI_Type_free(&MPI_metanode_type);
	MPI_File_close(&fp);
	MPI_Finalize();
	return 0;
}
