/**
		F.H.P.C. Assingment 3
		@file evolution_omp.cc
		@brief omp processes evolve the particles system
		
		@author Pietro Morichetti
		@date 20/01/2020
		@version 1.1 
*/

#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <math.h>

#define REF_THREAD 0 // shows the operations performed by one specific thread
#define SUMMARY 1 // shows the results information of the particles system (final positions, final velocity, ....)
#define UPDATE 0 // shows all operations performed on the particles (positions, velocity, ....)
#define CHECK 0 // shows all support informations (malloc, group, comm, ....)

typedef struct {double coord[3]{0, 0, 0};} vect;
typedef struct {vect pos; vect vel; double E;} particle;
typedef struct {int info[3]{8, 0, 0}; double evolution_time[1]{0};} metanode; // metanode of the file
typedef struct { // universal constants and differents measures of system times
	const double M_max {100};
	double m {1};
	const double G {pow(10, -6)};
	const int iter_max {3};
	const double epsilon {0.05};
	
	double delta_t_min {0};
	double delta_t {0};
	double tmp_delta_t_min [1]{0};
	double evolution_time {0};
	const int time_to_write_checkpoint {3}; // each 3 evolution of the system save a checkpoint file
} evolution_paramiter;

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

void force_computing(const int& particle_index, const int& start_point, const int& end_point,  
										 const particle* array_all_particle, particle* array_particle, 
										 const int& myid, vect& force_i){
	
	vect tmp_coord;
	double tmp_force_coord_norm;
	
	for(int j = start_point; j < end_point; ++j){ // updates system considering "all" particles
		// r2
		if(UPDATE) if(myid == REF_THREAD) std::cout << end_point << "-" << j << ": r2 = " << array_all_particle[j].pos;
		// r_1 - r_2
		tmp_coord = array_all_particle[particle_index].pos - array_all_particle[j].pos;
		if(UPDATE) if(myid == REF_THREAD) std::cout << end_point << "-" << j << ": r1 - r2 = " << tmp_coord;
		
		// || r_1 - r_2 ||^3
		tmp_force_coord_norm = norm(tmp_coord)*norm(tmp_coord)*norm(tmp_coord);
		if(UPDATE) if(myid == REF_THREAD) std::cout << end_point << "-" << j << ": ||r1 - r2||^3 = " <<   tmp_force_coord_norm << "\n";
		
		// (r_1 - r_2)/(|| r_1 - r_2 ||^3)
		tmp_coord /= tmp_force_coord_norm;
		if(UPDATE) if(myid == REF_THREAD) std::cout << end_point << "-" << j <<": (r1 - r2)/||r1 - r2||^3 = " << tmp_coord;					
		
		// partial force on the single particle 
		force_i +=  tmp_coord;
		if(UPDATE) if(myid == REF_THREAD) std::cout << end_point << "-" << j << ": partial force = " << force_i;		
		
		// for the energy: sum of 1/(|| r_1 - r_2 ||)
		array_particle[particle_index].E += 1/norm(tmp_coord);
		if(UPDATE) if(myid == REF_THREAD) std::cout << end_point << "-" << j << ": partial energy = " << array_particle[particle_index].E << "\n";		
	}
}

void evolution_computation(vect& force_i, vect& acceleration_i, const evolution_paramiter& system, 
													 const particle* array_all_particle, particle* array_particle,
													 const int& N_p, const int& myid, const int& i, const int& start_particle){
	
	// force of the single particle
	force_i *= (N_p - 1)*system.m*system.m*system.G;
	if(SUMMARY) if(myid == REF_THREAD) std::cout << i << ": force = " << force_i;
	
	// acceleration of the single particle
	acceleration_i = force_i/system.m;
	if(SUMMARY) if(myid == REF_THREAD) std::cout << i << ": acc = " << acceleration_i;
	
	// delta velocity of the single particle
	array_particle[i - start_particle].vel =  array_all_particle[i].vel + acceleration_i*system.delta_t;
	if(SUMMARY) if(myid == REF_THREAD) std::cout << i << ": velocity = "  << array_particle[i - start_particle].vel;
	
	// delta position of the single particle
	array_particle[i - start_particle].pos = array_all_particle[i].pos + (array_particle[i - start_particle].vel + array_all_particle[i].vel)*system.delta_t;
	if(SUMMARY) if(myid == REF_THREAD) std::cout << i << ": position = " << array_particle[i - start_particle].pos;
	
	// energy of the single particle
	array_particle[i - start_particle].E = 0.5*norm(array_particle[i - start_particle].vel)*norm(array_particle[i - start_particle].vel) + system.G*(N_p-1)*system.m*array_particle[i - start_particle].E;
	if(SUMMARY) if(myid == REF_THREAD) std::cout << i << ": energy = " << array_particle[i - start_particle].E << "\n";
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

int main(int argc, char *argv[]){
	
	metanode meta_data{};
	//------ READING NUMBER OF PARTICLES FROM THE FILE -----//
	FILE *fp = fopen ( "data_ic" , "r" );
	fseek(fp, 4, SEEK_SET );
	fread(&meta_data.info[1], sizeof(int), 1, fp);
	fclose(fp);
	int N_p {meta_data.info[1]}; // number of particles 
	//------------------------------------------------------------//
	
	// get the maximum number of threads (it's equivalent to OMP_NUM_THREADS setted in the environment), and the optimal one
	int N_t {omp_get_max_threads()};
	
	// number of particles managed by one single processor
	int num_part_per_task {N_p/N_t}; 
	int R {N_p%N_t};
		
	// mass of the system
	evolution_paramiter system {};
	system.m = system.M_max/N_p;
	
	particle * array_all_particle = (particle*)calloc(N_p, sizeof(particle));
	
	#pragma omp parallel num_threads(N_t) firstprivate(num_part_per_task)// START PARALLEL REGION
	{
		int myid {omp_get_thread_num()};
		double potenzial_delta_t_min {0};
		if(myid < R){++num_part_per_task;} 
		particle * array_particle = (particle*)calloc(num_part_per_task, sizeof(particle));
		
		// interval of particles for each process
		int start_particle {(myid < R) ? myid*num_part_per_task : myid*num_part_per_task + R};
		
		// start read/write point in the file, for each thread
		int offset {(myid < R) ? 20 + myid*num_part_per_task*(int)sizeof(particle) : 
												    20 + (myid*num_part_per_task + R)*(int)sizeof(particle)};
		
		if(SUMMARY) if(myid == REF_THREAD) 
				 std::cout <<"The speaker is thread " << REF_THREAD << "\n"
								  << "N_p = " << N_p << ", N_t = " << N_t  
								  << ", num_part_per_task = " << num_part_per_task
								  << ", R = " << R  << "\n";
		
		//open file
		FILE * fp_1, *fp_2;
		fp_1 = fopen( "data_ic" , "r" );
		fp_2 = fopen("check_point", "w+");
		
		for(int iter = 0; iter < system.iter_max; ++iter){
			
			// everybody read the file first time and only this time (collective thread operation)
			if(iter == 0){			
				fseek(fp_1 , offset, SEEK_SET);
				fread(&array_all_particle[start_particle], sizeof(particle), num_part_per_task, fp_1);
			}
			
			#pragma omp barrier
			
			if(myid == REF_THREAD) print_array(array_all_particle, N_p);
			system.tmp_delta_t_min[0] = 5; // define the maximum interval of time between two sequence of evolution
			
			// customization of the omp reduction to apply it on the struct type
			#pragma omp declare reduction(min : evolution_paramiter : \
					omp_out.tmp_delta_t_min[0] = (omp_in.tmp_delta_t_min[0] < omp_out.tmp_delta_t_min[0]) ? omp_in.tmp_delta_t_min[0] : omp_out.tmp_delta_t_min[0]) \
					initializer(omp_priv(omp_orig))	
			
			#pragma omp for reduction(min : system) schedule(static)
			for(int i = 0; i < N_p; ++i){ // makes for each particle that has to be considered
				vect force_i, acceleration_i;
				
				if(SUMMARY) if(myid == REF_THREAD) std::cout << i << ": read = " << array_all_particle[i].pos;
				
				//----------- COMPUTING FORCE CONSIDERING ALL PARTICLES ------------//			
				if(CHECK) if(myid == REF_THREAD) std::cout << "			start force 1 loop\n";
				
				force_computing(i, 0, i,  array_all_particle, array_particle, myid, force_i);
				force_computing(i, i+1, N_p,  array_all_particle, array_particle, myid, force_i);
				
				if(CHECK) if(myid == REF_THREAD) std::cout << "			end force 2 loop\n";	
				
				//----------- COMPUTING THE CHARACTERISTICS OF THE PARTICLE ------------//
				evolution_computation(force_i, acceleration_i, system, array_all_particle, array_particle, N_p, myid, i, start_particle);
				
				// new potential delta_t for the next iteration
				potenzial_delta_t_min = system.epsilon*norm(array_particle[i - start_particle].vel)/norm(acceleration_i);
				if(potenzial_delta_t_min < system.tmp_delta_t_min[0]) system.tmp_delta_t_min[0] = potenzial_delta_t_min;
			}	// END THREADS COMPUTATION PART
			
			#pragma omp master // define next delta t and update the history time of the system
			{
				system.delta_t = 0.5*system.tmp_delta_t_min[0];
				system.evolution_time += system.delta_t;
			}
			
			// update system
			#pragma omp for schedule(static)
			for(int i = 0; i < N_p; ++i){	
				array_all_particle[i] = array_particle[i - start_particle];
			}
			
			if(iter%system.time_to_write_checkpoint == 0){ // it's time to write a checkpoint file
				
				#pragma omp master // master thread updates the meta_data of the file
				{
					fseek(fp_2 , 0, SEEK_SET );
					fwrite(&meta_data, sizeof(metanode), 1, fp_2);
				}
				
				// everybody move our pointer to a specific position and starts to write
				fseek(fp_2, offset, SEEK_SET);
				fwrite(array_particle, sizeof(particle), num_part_per_task, fp_2);			
			}
		}
	
		fclose(fp_1);
		fclose(fp_2);
		delete[] array_particle;
	}
	
	delete[] array_all_particle;
	return 0;
}