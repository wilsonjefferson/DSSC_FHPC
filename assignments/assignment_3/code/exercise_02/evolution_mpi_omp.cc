/**
		F.H.P.C. Assingment 3
		@file evolution_mpi_omp.cc
		@brief mpi and omp processes collaborate in the evolution of the particles system
		
		@author Pietro Morichetti
		@date 20/01/2020
		@version 1.1 
*/

#include <iostream>
#include <omp.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>

#define SOCKET_NUM_THREAD_MAX 9 // number of cores allowed for one CPU
#define REF_PROCESS 0 // shows the operations performed by one specific mpi process
#define REF_THREAD 0 // shows the operations performed by one specific thread
#define SUMMARY 1 // shows the results information of the particles system (final positions, final velocity, ....)
#define UPDATE 0 // shows all operations performed on the particles (positions, velocity, ....)
#define CHECK 0 // shows all support informations (malloc, group, comm, ....)

typedef struct {double coord[3]{0, 0, 0};} vect; // spatial coordinates
typedef struct {vect pos; vect vel; double E;} particle; // particle features
typedef struct {int info[3]{8, 0, 0}; double evolution_time[1]{0};} metanode;
typedef struct { // universal constants and differents measures of system times
	const double M_max {100};
	double m {1};
	const double G {pow(10, -6)};
	const int iter_max {3};
	const double epsilon {0.05};
	
	double delta_t_min[1]{0};
	double delta_t[1]{0};
	double tmp_delta_t_min[1]{0};
	const int time_to_write_checkpoint {3}; // each 3 evolution of the system save a checkpoint file
} evolution_paramiter;
typedef struct {int size_array_group[4] {-1, -1, -1, -1}; int *array_group[4]{};} mpi_group_parameter; // groups variables

vect& operator+=(vect& p1, const vect& p2){
	p1.coord[0] += p2.coord[0];
	p1.coord[1] += p2.coord[1];
	p1.coord[2] += p2.coord[2];
	return p1;
}

vect operator-(const vect& p1, const vect& p2){
	vect tmp;
	tmp.coord[0] = p1.coord[0] - p2.coord[0];
	tmp.coord[1] = p1.coord[1] - p2.coord[1];
	tmp.coord[2] = p1.coord[2] - p2.coord[2];
	return tmp;
}

vect& operator/=(vect& p1, const double& p2){
	p1.coord[0] /= p2;
	p1.coord[1] /= p2;
	p1.coord[2] /= p2;
	return p1;
}

vect& operator*=(vect& p1, const double& p2){
	p1.coord[0] *= p2;
	p1.coord[1] *= p2;
	p1.coord[2] *= p2;
	return p1;
}

vect operator/(const vect& p1, const double& p2){
	vect tmp;
	tmp.coord[0] = p1.coord[0]/p2;
	tmp.coord[1] = p1.coord[1]/p2;
	tmp.coord[2] = p1.coord[1]/p2;
	return tmp;
}

vect operator+(const vect& p1, const vect& p2){
	vect tmp;
	tmp.coord[0] = p1.coord[0] + p2.coord[0];
	tmp.coord[1] = p1.coord[1] + p2.coord[1];
	tmp.coord[2] = p1.coord[2] + p2.coord[2];
	return tmp;
}

vect operator*(const vect& p1, const double& p2){
	vect tmp;
	tmp.coord[0] = p1.coord[0]*p2;
	tmp.coord[1] = p1.coord[1]*p2;
	tmp.coord[2] = p1.coord[2]*p2;
	return tmp;
}

std::ostream& operator<<(std::ostream& os, const vect& p){
	os << p.coord[0] << " " << p.coord[1] << " " << p.coord[2] << "\n";
	return os;
}

double norm(const vect& p){
	return pow(p.coord[0]*p.coord[0] + p.coord[1]*p.coord[1] + p.coord[2]*p.coord[2], 0.5);
}

int check_workers(const int& N_t_process, int&  N_t_thread, const int&  N_p, const int&  myid_process){
	
	int N_t_thread_optimal {(N_p/N_t_process < SOCKET_NUM_THREAD_MAX) ? N_p/N_t_process : SOCKET_NUM_THREAD_MAX};

	if(N_t_process*N_t_thread > N_p && N_t_thread == 1){ // catch exceptions
		if(myid_process == REF_PROCESS){
			std::cout << "------------------------------------------------------------\n"
								 << "------------------------------------------------------------\n\n"
								 << "ERROR: the dimension of the problem is less than the number of cores, it means some cores are unused.\n" 
								 << "Set an appropriate number of \"workers\" for this problem.\n\n"
								 << "------------------------------------------------------------\n"
								<< "------------------------------------------------------------\n";}
			return 1;
	} else if(N_t_process*N_t_thread > N_p && N_t_thread > SOCKET_NUM_THREAD_MAX && N_t_process*SOCKET_NUM_THREAD_MAX > N_p){
		N_t_thread = N_t_thread_optimal;
		if(myid_process == REF_PROCESS) 
			std::cout << "------------------------------------------------------------\n"
								 << "------------------------------------------------------------\n\n"
								 << "WARNING: the execution will involve in hyper-threading mode and it could be inefficient.\n" 
								 << "The program sets the number of threads equal to " << N_t_thread << " (OPTIMAL)\n\n"
								 << "------------------------------------------------------------\n"
								 << "------------------------------------------------------------\n";
	}else if(N_t_thread != N_t_thread_optimal){
		if(myid_process == REF_PROCESS) 
			std::cout << "------------------------------------------------------------\n"
								 << "------------------------------------------------------------\n\n"
								 << "WARNING: number of threads is not optimal for this problem.\n" 
								 << "It could be better to set it to "<< N_t_thread_optimal << " for this specific problem.\n\n"
								 << "------------------------------------------------------------\n"
								 << "------------------------------------------------------------\n";
	}
	return 0;
}

void build_mpi_data_type(MPI_Datatype& MPI_particle_type, MPI_Datatype& MPI_metanode_type){
	
	// build MPI_particle_type
	const std::size_t  num_members_1 {3};
	int blocklengths_1[num_members_1] = {3, 3, 1};
	
	// array of addresses
	MPI_Aint displacements_1[num_members_1] = {offsetof(particle, pos), offsetof(particle, vel), offsetof(particle, E)};
	// array of Data types
	MPI_Datatype types_1[num_members_1] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	
	MPI_Type_create_struct(num_members_1, blocklengths_1, displacements_1, types_1, &MPI_particle_type);
	MPI_Type_commit(&MPI_particle_type);
	
	//----------
	
	// build MPI_metanode_type
	const std::size_t  num_members_2 {2};
	int blocklengths_2[num_members_2] = {3, 1};
	
	MPI_Aint displacements_2[num_members_2] = {offsetof(metanode, info), offsetof(metanode, evolution_time)};
	MPI_Datatype types_2[num_members_2] = {MPI_INTEGER, MPI_DOUBLE};
	
	MPI_Type_create_struct(num_members_2, blocklengths_2, displacements_2, types_2, &MPI_metanode_type);
	MPI_Type_commit(&MPI_metanode_type);
}

int specific_group(const int& myid_process, const int& R_process, const int R_artificial, const int size_array_group, const int value_array_group, 
										  const int i, mpi_group_parameter& group_magnitude){
	
	if (myid_process < R_process && R_artificial > 0){
		group_magnitude.size_array_group[i] = size_array_group;
		group_magnitude.array_group[i] = (int*)malloc(group_magnitude.size_array_group[i]*sizeof(int));
		if(group_magnitude.array_group[i] == NULL){
			if(CHECK) if (myid_process == REF_PROCESS) std::cout << "ERROR: *group_magnitude.array_group[" << i << "] was not allocated.\n";
		}
		if(CHECK) if (myid_process == REF_PROCESS) std::cout << "group_magnitude.array_group[" << i << "]: ";
		for(int j = 0; j < group_magnitude.size_array_group[i]; ++j){
			group_magnitude.array_group[i][j] = j + value_array_group;
			if(CHECK) if (myid_process == REF_PROCESS) std::cout << group_magnitude.array_group[i][j] << " ";
		} 
		if(CHECK) if (myid_process == REF_PROCESS) std::cout << "\n";
	}else{return 0;}
	return 1;
}

int build_group(const int& N_t_process, const int& myid_process, const int& R_process, 
								   MPI_Group* group, mpi_group_parameter& group_magnitude,MPI_Comm* comm){
	
	// setting candidates for each subgroups
//	                                 		<myid_process>,   	  		<R_process>, 			<R_artificial>,				<size_array_group>,     <value_array_group>, 	<index_group>, <group_magnitude>
	if(specific_group(        myid_process,   	  	R_process, 							1,						            R_process,                      0, 				    0, 		group_magnitude) == 1)
		if(CHECK) if (myid_process == REF_PROCESS) std::cout << "WARNING: MPI process " << myid_process << " is able to enter in the " << 0 << " group.\n";

	if(specific_group(-myid_process - 1, 	   	  -R_process, 							1,	      N_t_process - R_process,		  R_process,    		    1, 		group_magnitude) == 1)
		if(CHECK) if (myid_process == REF_PROCESS) std::cout << "WARNING: MPI process " << myid_process << " is able to enter in the " << 1 << " group.\n";
		
	if(specific_group(   myid_process- 1,    	    R_process,			R_process,			      		  	  R_process + 1, 						  0, 	 			2, 		group_magnitude) == 1)
		if(CHECK) if (myid_process == REF_PROCESS) std::cout << "WARNING: MPI process " << myid_process << " is able to enter in the " << 2 << " group.\n";
		
	if(specific_group(      -myid_process, 	-R_process + 2, 			R_process,	N_t_process - R_process + 1,	R_process - 1,	 		    3, 		group_magnitude) == 1)
		if(CHECK) if (myid_process == REF_PROCESS) std::cout << "WARNING: MPI process " << myid_process << " is able to enter in the " << 3 << " group.\n";
	
	// define subgroups and communicators
	for(int i = 0; i < 4; ++i){
		if(group_magnitude.size_array_group[i] != -1){
			MPI_Group_incl(group[0], group_magnitude.size_array_group[i], group_magnitude.array_group[i], &group[i+1]);
			MPI_Comm_create_group(MPI_COMM_WORLD, group[i+1], 0, &comm[i]);
		}
	}
	
	return 0;
}

void deallocate_memories(MPI_Group* group, MPI_Datatype& MPI_particle_type, MPI_Datatype& MPI_metanode_type, MPI_File& fp_1, MPI_File& fp_2,
												 const particle* array_particle, const particle* array_all_particle,
												 mpi_group_parameter& group_magnitude,MPI_Comm* comm){
	
	for(int i = 0; i < 4; ++i){
		if(group_magnitude.size_array_group[i] != -1){
				MPI_Group_free(&group[i+1]); MPI_Comm_free(&comm[i]); delete [] group_magnitude.array_group[i];
		}
	}
	
	MPI_Group_free(&group[0]);
	MPI_Type_free(&MPI_particle_type);
	MPI_Type_free(&MPI_metanode_type);
	
	MPI_File_close(&fp_1);
	MPI_File_close(&fp_2);
	
	delete[] array_particle;
	delete[] array_all_particle;
}

void force_computing(const int& particle_index, const int& start_point, const int& end_point,  
									 const particle* array_all_particle, particle* array_all_particle_process, 
									 const int& myid_process, const int& myid_thread, vect& force_i){
	
	vect tmp_coord;
	double tmp_force_coord_norm;
	
	for(int j = start_point; j < end_point; ++j){ // updates system considering "all" particles
		// r2
		if(UPDATE)
				if(myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << end_point << "-" << j << ": r2 = " << array_all_particle[j].pos;
		// r_1 - r_2
		tmp_coord = array_all_particle[particle_index].pos - array_all_particle[j].pos;
		if(UPDATE)
			if(myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << end_point << "-" << j << ": r1 - r2 = " << tmp_coord;
		
		// || r_1 - r_2 ||^3
		tmp_force_coord_norm = norm(tmp_coord)*norm(tmp_coord)*norm(tmp_coord);
		if(UPDATE)
			if(myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << end_point << "-" << j << ": ||r1 - r2||^3 = " <<   tmp_force_coord_norm << "\n";
		
		// (r_1 - r_2)/(|| r_1 - r_2 ||^3)
		tmp_coord /= tmp_force_coord_norm;
		if(UPDATE)
			if(myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << end_point << "-" << j <<": (r1 - r2)/||r1 - r2||^3 = " << tmp_coord;					
		
		// partial force on the single particle 
		force_i +=  tmp_coord;
		if(UPDATE)
			if(myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << end_point << "-" << j << ": partial force = " << force_i;		
		
		// for the energy: sum of 1/(|| r_1 - r_2 ||)
		array_all_particle_process[particle_index].E += 1/norm(tmp_coord);
		if(UPDATE)
			if(myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << end_point << "-" << j << ": partial energy = " << array_all_particle_process[particle_index].E << "\n";		
	}
}

void evolution_computation (vect& force_i, vect& acceleration_i, const evolution_paramiter& system, 
												 const particle* array_all_particle, particle* array_all_particle_process,
												 const int& N_p, const int& myid_process, const int& myid_thread, const int& i, const int& start_particle){
	
	// force of the single particle
	force_i *= (N_p - 1)*system.m*system.m*system.G;
	if(SUMMARY)
		if(myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << i << ": force = " << force_i;
	
	// acceleration of the single particle
	acceleration_i = force_i/system.m;
	if(SUMMARY)
		if( myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << i << ": acc = " << acceleration_i;
	
	// delta velocity of the single particle
	array_all_particle_process[i - start_particle].vel =  array_all_particle[i].vel + acceleration_i*system.delta_t[0];
	if(SUMMARY)
		if(myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << i << ": velocity = "  << array_all_particle_process[i - start_particle].vel;
	
	// delta position of the single particle
	array_all_particle_process[i - start_particle].pos = array_all_particle[i].pos + (array_all_particle_process[i - start_particle].vel + array_all_particle[i].vel)*system.delta_t[0];
	if(SUMMARY)
		if(myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << i << ": position = " << array_all_particle_process[i - start_particle].pos;
	
	// energy of the single particle
	array_all_particle_process[i - start_particle].E = 0.5*norm(array_all_particle_process[i - start_particle].vel)*norm(array_all_particle_process[i - start_particle].vel) + system.G*(N_p-1)*system.m*array_all_particle_process[i - start_particle].E;
	if(SUMMARY)
		if(myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << i << ": energy = " << array_all_particle_process[i - start_particle].E << "\n";
}

void print_array(const particle* p, const int size){
	std::cout << "------- print array of positions -------\n";
	for(int i = 0; i < size; ++i){
		std::cout << i << ": ";
		for(int j = 0; j<3; ++j){
			std::cout << p[i].pos.coord[j] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "----------------------------------------\n";
}

int main(int argc, char *argv[]){ // --------------------- MAIN ---------------------
	
	int myid_process, N_t_process {0};
	
	//------- DEFINE MPI SETTING ------//
	int mpi_provided_thread_level;
	MPI_Group group[5];
	MPI_Comm comm[4];
	MPI_Datatype MPI_particle_type, MPI_metanode_type;
	MPI_File fp_1, fp_2;
	MPI_Status status;		
	
	MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_thread_level);
	if ( mpi_provided_thread_level < MPI_THREAD_FUNNELED ) {
		printf("a problem arise when asking for MPI_THREAD_FUNNELED level\n");
		MPI_Finalize();
		exit( 1 );
	}
	MPI_Comm_size(MPI_COMM_WORLD,&N_t_process);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid_process);
	MPI_Comm_group(MPI_COMM_WORLD, &group[0]);
	//---------------------------------//
	
	// open file
	if(MPI_File_open(MPI_COMM_WORLD, "data_ic", MPI_MODE_RDONLY, MPI_INFO_NULL, &fp_1) == 1){
		if(CHECK) if(myid_process == REF_PROCESS) std::cout << "ERROR: fp_1 is not pointing to the \"data_ic\" file.\n";
		MPI_Finalize(); return 1;
	}
	if(MPI_File_open(MPI_COMM_WORLD, "check_point", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_2) == 1){
		if(CHECK) if(myid_process == REF_PROCESS) std::cout << "ERROR: fp_2 is not pointing to the \"check_point\" file.\n";
		MPI_Finalize(); return 1;		
	}
	
	//------- DEFINE NEW DATATYPE ------//
	build_mpi_data_type(MPI_particle_type, MPI_metanode_type);	
	
	// information about the system
	metanode meta_data {}; // first block of bytes that contain general information about the data_ic file
	
	// setting number of particles
	MPI_File_read_all(fp_1, &meta_data, 1, MPI_metanode_type, &status);
	int N_p {meta_data.info[1]};
	
	// setting number of particles per process
	int num_part_per_process {N_p/N_t_process}; // number of particles managed by one single processor
	int R_process {N_p%N_t_process};

	// variables used to perform the collective sending operations
	int num_part_per_process_smaller {num_part_per_process};
	int num_part_per_process_greater {(R_process == 0) ? num_part_per_process : num_part_per_process + 1};

	if(myid_process < R_process){++num_part_per_process;}
	
	// interval of particles for each process
	int start_particle {(myid_process < R_process) ? myid_process*num_part_per_process : myid_process*num_part_per_process + R_process};
	int end_particle {(myid_process < R_process) ? num_part_per_process*(myid_process + 1) : (myid_process + 1)*num_part_per_process + R_process};	
	
	// get the maximum number of threads (it's equivalent to OMP_NUM_THREADS setted in the environment), and the optimal one
	int N_t_thread {omp_get_max_threads()};
	
	//-------------------- CHECKING NUMBER OF "WORKERS" --------------------//
	if(check_workers(N_t_process, N_t_thread, N_p, myid_process) == 1) {MPI_Finalize(); return 1;}
	
	// setting number of particles per thread
	int R_thread {num_part_per_process%N_t_thread};
	
	//--------------------------- BUILDING GROUPS ---------------------------//
	mpi_group_parameter group_magnitude{}; // default constructor
		
	if(build_group(N_t_process, myid_process, R_process, group, group_magnitude, comm) == 1) {MPI_Finalize(); return 1;}
	MPI_Barrier(MPI_COMM_WORLD); // to be sure that all processes have builded the groups
	
	// times of the particles system
	evolution_paramiter system{}; // default constructor
	system.m = system.M_max/N_p;
	
	// setting different array
	particle *array_all_particle = (particle*)calloc(N_p, sizeof(particle));
	if(array_all_particle == NULL){
		if(CHECK) if(myid_process == REF_PROCESS) std::cout << "ERROR: *array_all_particle was not allocated.\n";
		MPI_Finalize(); return 1;}
		
	particle *array_all_particle_process = (particle*)calloc(num_part_per_process, sizeof(particle));
	if(array_all_particle_process == NULL){
		if(CHECK) if(myid_process == REF_PROCESS) std::cout << "ERROR: *array_all_particle_process was not allocated.\n";
		MPI_Finalize(); return 1;}
	
	// setting the start read/write in the file, for each process
	MPI_Offset offset_process {(myid_process < R_process) ? 20 + myid_process*num_part_per_process*(int)sizeof(particle) : 
																									          20 + (myid_process*num_part_per_process + R_process)*(int)sizeof(particle)};
																									
	if(CHECK) if(myid_process == REF_PROCESS) std::cout << "start thread region\n";
	
	#pragma omp parallel num_threads(N_t_thread) // START PARALLEL REGION
	{	
		int myid_thread {omp_get_thread_num()};
		double potenzial_delta_t_min {0};
		
		if(SUMMARY)
			if(myid_thread == REF_THREAD && myid_process == REF_PROCESS) 
				 std::cout <<"The speaker is thread " << REF_THREAD << " in process " << REF_PROCESS << "\n"
									  << "N_p = " << N_p << ", N_t_process = " << N_t_process << ", N_t_thread = " << N_t_thread << "\n"
									  << "R_process = " << R_process << ", R_thread = " << R_thread << "\n"
									  << "num_part_per_process = " << num_part_per_process	<< "\n"
									  << "num_part_per_process_smaller = " << num_part_per_process_smaller << "\n"
									  << "num_part_per_process_greater = " << num_part_per_process_greater << "\n"
									  << "num_part_per_process_greater = " << num_part_per_process_greater << "\n"
									  << "start_particle_process = " << start_particle << ", end_particle_process = " << end_particle << "\n";
		
		if(CHECK) if(myid_thread == REF_THREAD && myid_process == REF_PROCESS) std::cout << "	start iter loop\n";	
		
		for(int iter = 0; iter < system.iter_max; ++iter){
		
			#pragma omp master // each master thread reads the file
			if(iter == 0) MPI_File_read_at_all(fp_1, 20, array_all_particle, N_p, MPI_particle_type, &status);
			
			system.tmp_delta_t_min[0] = 5; // define the maximum interval of time between two sequence of evolution
			if(SUMMARY) if(myid_process == REF_PROCESS && myid_thread == REF_THREAD) print_array(array_all_particle, N_p);			
			if(CHECK) if(myid_thread == REF_THREAD && myid_process == REF_PROCESS) std::cout << "		start omp loop\n";	
			
			// customization of the omp reduction to apply it on the struct type
			#pragma omp declare reduction(min : evolution_paramiter : \
					omp_out.tmp_delta_t_min[0] = (omp_in.tmp_delta_t_min[0] < omp_out.tmp_delta_t_min[0]) ? omp_in.tmp_delta_t_min[0] : omp_out.tmp_delta_t_min[0]) \
					initializer(omp_priv(omp_orig))	
			
			#pragma omp for reduction(min : system) schedule(static)
			for(int i = start_particle; i < end_particle; ++i){ // makes for each particle that has to be considered
				vect force_i, acceleration_i;
				
				if(SUMMARY) if(myid_process == REF_PROCESS &&  myid_thread == REF_THREAD) std::cout << i << ": read = " << array_all_particle[i].pos;																							   
				
				//----------- COMPUTING FORCE CONSIDERING ALL PARTICLES ------------//			
				if(CHECK) if(myid_thread == REF_THREAD && myid_process == REF_PROCESS) std::cout << "			start force 1 loop\n";
				
				force_computing(i - start_particle, start_particle, i,  array_all_particle, array_all_particle_process, 
												myid_process, myid_thread, force_i);
				
				force_computing(i - start_particle, i+1, end_particle,  array_all_particle, array_all_particle_process, 
												myid_process, myid_thread, force_i);
				
				if(CHECK) if(myid_thread == REF_THREAD && myid_process == REF_PROCESS) std::cout << "			end force 2 loop\n";	
				
				//----------- COMPUTING THE CHARACTERISTICS OF THE PARTICLE ------------//
				evolution_computation(force_i, acceleration_i, system, 
															array_all_particle, array_all_particle_process,
															N_p,  myid_process, myid_thread, i, start_particle);
				
				// new potential delta_t for the next iteration
				potenzial_delta_t_min = system.epsilon*norm(array_all_particle_process[i - start_particle].vel)/norm(acceleration_i);
				if(SUMMARY) if(myid_process == REF_PROCESS && myid_thread == REF_THREAD) std::cout << i << ": potenzial_delta_t_min = " << potenzial_delta_t_min << "\n";
				if(SUMMARY)if(myid_process == REF_PROCESS && myid_thread == REF_THREAD) std::cout << i << ": system.tmp_delta_t_min = " <<  system.tmp_delta_t_min[0] << "\n";
				if(potenzial_delta_t_min < system.tmp_delta_t_min[0]) system.tmp_delta_t_min[0] = potenzial_delta_t_min;
				if(SUMMARY) if(myid_process == REF_PROCESS && myid_thread == REF_THREAD) std::cout << i << ": system.tmp_delta_t_min = " <<  system.tmp_delta_t_min[0] << "\n";
			}	// END THREADS COMPUTATION PART
						
			if(CHECK) if(myid_thread == REF_THREAD && myid_process == REF_PROCESS) std::cout << "		end omp loop\n";	
									
			#pragma omp master // defines next delta t and update the history time of the system
			{
				// master threads determine the lower delta t min from all the particles under the domain of the own process
				MPI_Allreduce(system.tmp_delta_t_min, system.delta_t_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
				
				system.delta_t[0] = 0.5*system.delta_t_min[0];
				meta_data.evolution_time[0] += system.delta_t[0];
			}
			
			if(SUMMARY) if(myid_process == REF_PROCESS && myid_thread == REF_THREAD) std::cout << "delta_t = " << system.delta_t[0] << "\n";
			
			// update system
			#pragma omp for schedule(static)
			for(int i = start_particle; i < end_particle; ++i){		
				array_all_particle[i] = array_all_particle_process[i - start_particle];
			}
			
			if(CHECK) if(myid_thread == REF_THREAD && myid_process == REF_PROCESS) std::cout << "		(after first barrier)\n";
						
			#pragma omp master // REGION PROTECTED FOR THE MASTER THREADS
			{
				//------------------- COLLECTIVE SUB GROUPS OPERATIONS -------------------//
				// master threads send the modified portion of the system to the other process
				
				// all gather operation in group_0 and group_1
				if (group_magnitude.size_array_group[0] != -1) MPI_Allgather(&array_all_particle[start_particle], num_part_per_process, MPI_particle_type, &array_all_particle[0], num_part_per_process, MPI_particle_type, comm[0]);
				if (group_magnitude.size_array_group[1] != -1) MPI_Allgather(&array_all_particle[start_particle], num_part_per_process, MPI_particle_type, &array_all_particle[N_p - num_part_per_process*(N_t_process - R_process)], num_part_per_process, MPI_particle_type, comm[1]);
				
				// the root node for group_2 and group_3 broadcasts the evaluated portion of array_all_particle in the group_0 and group_1 -> everyone have new array_all_particle updated
				if (group_magnitude.size_array_group[2] != -1) MPI_Bcast(&array_all_particle[N_p - num_part_per_process_smaller*(N_t_process - R_process)], num_part_per_process_smaller*(N_t_process - group_magnitude.size_array_group[2] + 1), MPI_particle_type, group_magnitude.size_array_group[2] - 1, comm[2]);
				if (group_magnitude.size_array_group[3] != -1) MPI_Bcast(&array_all_particle[0], num_part_per_process_greater*(N_t_process - group_magnitude.size_array_group[3] + 1), MPI_particle_type, 0, comm[3]);
				
				//-----------------------------------------------------------------------------//	
				
				if(iter%system.time_to_write_checkpoint == 0){ // it's time to write a checkpoint file
					MPI_File_write(fp_2, &meta_data, 1, MPI_metanode_type, &status);
					MPI_File_write_at_all(fp_2, offset_process, array_all_particle, num_part_per_process, MPI_particle_type, &status);			
				}
			} // END PROTECTED REGION FOR THE MASTER THREADS
			
			#pragma omp barrier // all threads wait masters sharing the particles system
			if(CHECK) if(myid_thread == REF_THREAD && myid_process == REF_PROCESS) std::cout << "		(after second barrier)\n";
		}
		
		if(CHECK) if(myid_thread == REF_THREAD && myid_process == REF_PROCESS) std::cout << "	end iter loop\n";	
	} // END PARALLLEL REGION
	
	if(CHECK) if(myid_process == REF_PROCESS) std::cout << "end thread region\n";	
	if(myid_process == REF_PROCESS) print_array(array_all_particle, N_p);
	
	//-------------------------- DEALLOCATE MEMORIES --------------------------//
	deallocate_memories(group, MPI_particle_type, MPI_metanode_type, fp_1, fp_2,
											array_all_particle_process, array_all_particle,
											group_magnitude, comm);
	MPI_Finalize();
	return 0;
}