#include "../../include/estimators.h"
#include "../../include/config.h"


#include <math.h>




double compute_sp_wave_function_value(double* pos, int i, int system_size, WavePackage* trial) {

	double box_length = trial->simulation_parameters[0];
	double half_box_length2 = (box_length * box_length) / 4.0;

	double wf_value = 0.0;

	for (int j = 0; j < system_size; j += 3) {

		if (i == j) {
			continue;
		}
		double dx = pos[i] - pos[j];
		double dy = pos[i+1] - pos[j+1];
		double dz = pos[i+2] - pos[j+2];

		dx -= box_length * rint(dx / box_length);
		dy -= box_length * rint(dy / box_length);
		dz -= box_length * rint(dz / box_length);

		double r_ij2 = dx * dx + dy * dy + dz * dz;

		if (r_ij2 < half_box_length2) {
			double r_ij = sqrt(r_ij2);
			wf_value += trial->wf(r_ij, trial->variational_parameters, trial->simulation_parameters);
		}
	}

	return wf_value;

}