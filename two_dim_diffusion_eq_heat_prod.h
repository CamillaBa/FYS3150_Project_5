#pragma once
#include <omp.h>
#include "vector.h"
#include "matrix.h"

class two_dim_diffusion_eq_heat_prod {

	matrix<double> boundary; // boundary conditions
	vector<double> x; // x-coordinate indexed by i
	vector<double> y; // y-coordinate indexed by i
	double         alpha;
	int            m;
	int            n;
	double(*Q) (double, double, double); // heat production

	// coordinate data
	double         t; // time
	double         dt; // time step

public:
	matrix<double> u;

	// constructor
	two_dim_diffusion_eq_heat_prod(matrix<double> u, double h, double dt, double(*Q) (double, double, double), double t);

	// copy constructor
	two_dim_diffusion_eq_heat_prod(const two_dim_diffusion_eq_heat_prod  & copy);

	// update u using explicit Euler
	void explicit_euler_update();
};


// constructor
two_dim_diffusion_eq_heat_prod::two_dim_diffusion_eq_heat_prod(matrix<double> u, double h, double dt, double(*Q) (double, double, double), double time = 0) {
	two_dim_diffusion_eq_heat_prod::u        = u;
	two_dim_diffusion_eq_heat_prod::boundary = u;
	two_dim_diffusion_eq_heat_prod::m        = u.m;
	two_dim_diffusion_eq_heat_prod::n        = u.n;
	two_dim_diffusion_eq_heat_prod::Q        = Q;
	two_dim_diffusion_eq_heat_prod::dt       = dt;
	two_dim_diffusion_eq_heat_prod::alpha    = dt / h / h;
	two_dim_diffusion_eq_heat_prod::t        = time;

	// fix matrix that represents boundary condition
	for (int i = 1; i < m - 1; i++) {
		for (int j = 1; j < n - 1; j++) {
			boundary[i][j] = 0;
		}
	}

	//fill x and y vectors
	x = vector<double>(n);
	y = vector<double>(n);
	for (int i = 0; i < m; i++) { x[i] = i * h; }
	for (int j = 0; j < n; j++) { x[j] = j * h; }
}

// copy constructor
two_dim_diffusion_eq_heat_prod::two_dim_diffusion_eq_heat_prod(const two_dim_diffusion_eq_heat_prod  & copy) {
	boundary = copy.boundary;
	u        = copy.u;
	alpha    = copy.alpha;
	m        = copy.m;
	n        = copy.n;
	Q        = copy.Q;
	t        = copy.t;
	dt       = copy.dt;
	x        = copy.x;
	y        = copy.y;
}

// update u using explicit Euler
void two_dim_diffusion_eq_heat_prod::explicit_euler_update() {
	matrix<double> u_temp = boundary;
	#pragma omp parallel for
	for (int i = 1; i < m - 1; i++) {
		for (int j = 1; j < n - 1; j++) {
			u_temp[i][j] = alpha * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4 * u[i][j]) + u[i][j] + dt * Q(x[i], y[j], t);
		}
	}
	u = u_temp;
	t += dt;
}