#include "../../include/physical_functions.h"
#include "../../include/shared.h"

#include <stdio.h>





double lennard_jones_potential(double r, double* simulation_params) {

	double r2 = (SIGMA / r) * (SIGMA / r);
	double r6 = r2 * r2 * r2;
	double r12 = r6 * r6;

	double value = 4.0 * EPS * (r12 - r6);

	return value;
}

double lj_potential_shifted(double r, double* simulation_params) {

	double box_side = simulation_params[0];
	double half_box_side = box_side / 2.0;

	double value = lennard_jones_potential(r, simulation_params) - lennard_jones_potential(half_box_side, simulation_params);
	return value;
}

double jastrow_u(double r, double* variational_params, double* simulation_params) {
	
	double b = variational_params[0];
	double beta = b * SIGMA;
	double beta5 = beta * beta * beta * beta * beta;
	double r5 = 1.0 / (r * r * r * r * r);
	//printf("distance r: %.4f\n", r);
	double value = beta5 * r5;

	return value;
}

double jastrow_u_derivative(double r, double* variational_params, double* simulation_params) {

	return 0;
}

double symmetrical_cutoff_u(double r, double* variational_params, double* simulation_params) {

	double box_side = simulation_params[0];
	double half_box_side = box_side / 2.0;
	double value = jastrow_u(r, variational_params, simulation_params) + jastrow_u(box_side - r, variational_params, simulation_params) - 2.0 * jastrow_u(half_box_side, variational_params, simulation_params);

	return value;
}

double jastrow_u_derivative2(double r, double* variational_params, double* simulation_params) {

	double b = variational_params[0];
	double beta = b * SIGMA;
	double beta5 = beta * beta * beta * beta * beta;
	double r5 = 1.0 / (r * r * r * r * r);

	double value = -30.0 * (beta5 * r5) / (r * r);

	return value;
}