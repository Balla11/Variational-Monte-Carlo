#include "../../include/monte_carlo.h"
#include "../../include/config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>


typedef struct {

	double (*wf)(double, double*, double*);
	double (*wf_d)(double, double*, double*);
	double (*wf_d2)(double, double*, double*);
	double* variational_parameters;
	double* system_parameters;

}WavePackage;

typedef struct {

	double laplacian_kinetic_energy;
	double potential_energy;
	double gradient_kinetic_energy;
	double* quantum_force_array;

}SystemObservables;

SystemObservables* observables_init(double* particle_positions, int system_size, WavePackage* trial) {

	SystemObservables* obs = malloc(sizeof(SystemObservables));
	obs->quantum_force_array = (double*)calloc(system_size, sizeof(double));
	obs->laplacian_kinetic_energy = 0.0;
	obs->potential_energy = 0.0;
	obs->gradient_kinetic_energy = 0.0;

}

void free_observables(SystemObservables* obs) {
	free(obs->quantum_force_array);
	free(obs);
}

void monte_carlo_simulation(SimConfig* simulation_parameters, WavePackage* trial, double* particle_positions, int system_size) {

	SystemObservables* obs = observables_init(particle_positions, system_size, trial);

	double accept_counter = 0.0;

	for (double t = 1.0; t < simulation_parameters->max_steps; t++) {

		metropolis_step(simulation_parameters, particle_positions, system_size, &accept_counter);

	}

	free_observables(obs);
}

void metropolis_step(SimConfig* sim_parameters, double* particle_positions, int system_size, double* accept_counter) {

	double* new_particle_positions = (double*)malloc(system_size * sizeof(double));
	arr_copy(particle_positions, new_particle_positions, system_size);

	int random_index = 3 * random_int(system_size / 3);

	propose_move(new_particle_positions, random_index, sim_parameters->displacement, sim_parameters->box_length);
	bool state = evaluate_move(particle_positions, new_particle_positions, random_index, system_size, accept_counter, trial);

	free(new_particle_positions);
}

void arr_copy(double* original_arr, double* copy_arr, int size) {

	for (int i = 0; i < size; i++) {
		copy_arr[i] = original_arr[i];
	}

}

void propose_move(double* pos, int random_index, double displacement, double box_length) {

	pos[random_index] += random_double(-displacement, displacement);
	pos[random_index+1] += random_double(-displacement, displacement);
	pos[random_index+2] += random_double(-displacement, displacement);

	check_pbc(pos, random_index, box_length);

}

void check_pbc(double* pos, int i, double box_length) {

	pos[i] -= box_length * rint(pos[i] / box_length);
	pos[i+1] -= box_length * rint(pos[i+1] / box_length);
	pos[i+2] -= box_length * rint(pos[i+2] / box_length);

}

bool evaluate_move(double* old_pos, double* new_pos, int i, int system_size, double* accept_counter, WavePackage* trial) {

	double old_wf_value = compute_sp_wave_function_value(old_pos, i, system_size, trial);
	double new_wf_value = compute_sp_wave_function_value(new_pos, i, system_size, trial);

	double probability = exp(old_wf_value - new_wf_value);
	double eta = random_double(0, 1);
	
	if (probability >= 1 || probability > eta) {
		
		accept_move(old_pos, new_pos, i);
		*accept_counter = *accept_counter + 1;
		return true;
	}
	return false;

}

int random_int(int max) {
	return rand() % max;
}

double random_double(double min, double max) {
	return min + (max - min) * ((double)rand() / RAND_MAX);
}

void accept_move(double* old_pos, double* new_pos, int index) {

	old_pos[index] = new_pos[index];
	old_pos[index + 1] = new_pos[index + 1];
	old_pos[index + 2] = new_pos[index + 2];

}