


typedef struct {

	double (*wf)(double, double*, double*);
	double (*wf_d)(double, double*, double*);
	double (*wf_d2)(double, double*, double*);
	double* variational_parameters;
	double* simulation_parameters;

}WavePackage;

WavePackage* wave_package_init(int vp_params_size, int sim_params_size);

typedef struct {

	int particles_per_side;
	int bravais_lattice_factor;
	double max_steps;
	double beta;
	WavePackage* trial_wave_function;

}SimConfig;




extern SimConfig simulation_config;