#include "../include/config.h"
#include "../include/shared.h"
#include "../include/monte_carlo.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void initialize_positions(double* positions, int num, int M, int N, double density);

int main(void) {

	// Choose variational parameters and simulation parameters arrays dimensions
	// I dont like this but let's keep it this way for now
	simulation_config.trial_wave_function = wave_package_init(1, 1);

	printf("Particles per side: %d\n", simulation_config.particles_per_side);
	printf("Box length: %.4f\n", simulation_config.trial_wave_function->simulation_parameters[0]);

	int particles_per_volume = simulation_config.particles_per_side * simulation_config.particles_per_side * simulation_config.particles_per_side;

	int system_size = 3 * particles_per_volume * simulation_config.bravais_lattice_factor;

	double* particles_positions = (double*)malloc(system_size * sizeof(double));
	initialize_positions(particles_positions, simulation_config.particles_per_side, simulation_config.bravais_lattice_factor, system_size, RHO);

	monte_carlo_simulation(simulation_config, simulation_config.trial_wave_function, particles_positions, system_size);

	free(particles_positions);
	return 0;
}

void initialize_positions(double* positions, int num, int M, int N, double density) {

	double b1[] = { 0, 0, 0 };
	double b2[] = { 0, 0, 0,  0.5, 0.5, 0.5 };
	double b4[] = { 0, 0, 0,  0.5, 0.5, 0,  0.5, 0, 0.5,  0, 0.5, 0.5 };


	double a = pow((double)M / density, 1. / 3.);
	int i = 0;


	if (M == 1) {
		for (int k = 0; k < num; k++) {
			for (int l = 0; l < num; l++) {
				for (int m = 0; m < num; m++) {
					for (int j = 0; j < M; j++) {
						positions[3 * i] = a * (k + b1[3 * j]);
						positions[3 * i + 1] = a * (l + b1[3 * j + 1]);
						positions[3 * i + 2] = a * (m + b1[3 * j + 2]);
						i = i + 1;
					}
				}
			}
		}
	}
	else if (M == 2) {
		for (int k = 0; k < num; k++) {
			for (int l = 0; l < num; l++) {
				for (int m = 0; m < num; m++) {
					for (int j = 0; j < M; j++) {
						positions[3 * i] = a * (k + b2[3 * j]);
						positions[3 * i + 1] = a * (l + b2[3 * j + 1]);
						positions[3 * i + 2] = a * (m + b2[3 * j + 2]);
						i = i + 1;
					}
				}
			}
		}
	}
	else if (M == 4) {
		for (int k = 0; k < num; k++) {
			for (int l = 0; l < num; l++) {
				for (int m = 0; m < num; m++) {
					for (int j = 0; j < M; j++) {
						positions[3 * i] = a * (k + b4[3 * j]);
						positions[3 * i + 1] = a * (l + b4[3 * j + 1]);
						positions[3 * i + 2] = a * (m + b4[3 * j + 2]);
						i = i + 1;
					}
				}
			}
		}
	}
	else {
		printf("\nERRORE M NON HA UN VALORE COMPATIBILE\n");
	}



	/*
	FILE* file = fopen("positions.txt", "w");
	for (int j = 0; j < N; j = j + 3) {
		//fprintf(file, "Particella %d, x: %f, y: %f, z: %f\n", j / 3, positions[j], positions[j + 1], positions[j + 2]);
		fprintf(file, "%f %f %f\n", positions[j], positions[j + 1], positions[j + 2]);
	}
	fclose(file);
	*/
}