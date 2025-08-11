
typedef struct {

	int particles_per_side;

}SimConfig;

typedef struct {

	double (*wf)(double, double*, double*);
	double (*wf_d)(double, double*, double*);
	double (*wf_d2)(double, double*, double*);
	double* variational_parameters;
	double* simulation_parameters;
}WavePackage;


extern struct SimConfig simulation_config;