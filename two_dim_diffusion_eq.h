#pragma once
#include <omp.h>
#include "matrix.h"

//================================================================================
// Two-dimensional diffusion equation
//================================================================================

class two_dim_diffusion_eq {

	matrix<double> boundary; // boundary conditions
	double         alpha;    // dt/h^2

	// dimensions of spacial solution
	int            m;
	int            n;

public:
	matrix<double> u;

	// constructor
	two_dim_diffusion_eq(matrix<double> u, double h, double dt);

	// copy constructor
	two_dim_diffusion_eq(const two_dim_diffusion_eq & copy);

	// update u using explicit Euler
	void explicit_euler_update();
};


// constructor
two_dim_diffusion_eq::two_dim_diffusion_eq(matrix<double> u, double h, double dt) {
	two_dim_diffusion_eq::u        = u;
	two_dim_diffusion_eq::boundary = u;
	two_dim_diffusion_eq::m        = u.m;
	two_dim_diffusion_eq::n        = u.n;
	two_dim_diffusion_eq::alpha    = dt / h / h;

	// find matrix that represents boundary condition
	for (int i = 1; i < m - 1; i++) {
		for (int j = 1; j < n - 1; j++) {
			boundary[i][j] = 0;
		}
	}
}

// copy constructor
two_dim_diffusion_eq::two_dim_diffusion_eq(const two_dim_diffusion_eq & copy) {
	boundary = copy.boundary;
	u        = copy.u;
	alpha    = copy.alpha;
	m        = copy.m;
	n        = copy.n;
}

// update u using explicit Euler
void two_dim_diffusion_eq::explicit_euler_update() {
	matrix<double> u_temp = boundary;

	#pragma omp parallel for
	for (int i = 1; i < m - 1; i++) {
		for (int j = 1; j < n - 1; j++) {
			u_temp[i][j] = alpha * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4 * u[i][j]) + u[i][j];
		}
	}
	u = u_temp;
}