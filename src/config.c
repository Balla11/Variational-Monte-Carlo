#include "../include/config.h"
#include "../include/physical_functions.h"
#include "../include/shared.h"

#include <stdlib.h>
#include <math.h>

SimConfig simulation_config = {
	.particles_per_side = 2,
	.bravais_lattice_factor = 4,
	.max_steps = 200000,
	.beta = 1.14,
	.trial_wave_function = NULL // Temporary solution sadly
};

WavePackage* wave_package_init(int vp_params_size, int sim_params_size) {

	WavePackage* wave_pack = malloc(sizeof(WavePackage));
	wave_pack->wf = jastrow_u;
	wave_pack->wf_d = jastrow_u;
	wave_pack->wf_d2 = jastrow_u;
	double* vp_params = (double*)malloc(vp_params_size * sizeof(double));
	double* sim_params = (double*)malloc(sim_params_size * sizeof(double));

	wave_pack->variational_parameters = vp_params;
	wave_pack->simulation_parameters= sim_params;
	// To be fixed
	vp_params[0] = simulation_config.beta;
	sim_params[0] = (double)simulation_config.particles_per_side * pow((double)simulation_config.bravais_lattice_factor / RHO, 1.0 / 3.0);
	

	return wave_pack;

};