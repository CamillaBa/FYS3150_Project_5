#pragma once
#include "vector.h"
#include "tridiag_solver.h"

class one_dim_diffusion_eq_heat_prod {
	vector<double> x; // x-coordinate indexed by i
	double         u0, uL;    // boundary conditions
	double         alpha;
	int            n;
	double(*Q) (double, double);


	// coordinate data
	double         t; // time
	double         dt; // time step

public:
	vector<double> u;

	// constructor
	one_dim_diffusion_eq_heat_prod(vector<double> u, double dx, double dt, double(*Q) (double, double), double time);

	// copy constructor
	one_dim_diffusion_eq_heat_prod(const one_dim_diffusion_eq_heat_prod & copy);

	// update u using explicit Euler
	void explicit_euler_update();
};

// constructor
one_dim_diffusion_eq_heat_prod::one_dim_diffusion_eq_heat_prod(vector<double> u, double dx, double dt, double(*Q) (double, double), double time = 0) {
	one_dim_diffusion_eq_heat_prod::u     = u;
	one_dim_diffusion_eq_heat_prod::n     = u.n;
	one_dim_diffusion_eq_heat_prod::u0    = u[0];
	one_dim_diffusion_eq_heat_prod::uL    = u[n - 1];
	one_dim_diffusion_eq_heat_prod::dt    = dt;
	one_dim_diffusion_eq_heat_prod::alpha = dt / dx / dx;
	one_dim_diffusion_eq_heat_prod::t     = time;
	one_dim_diffusion_eq_heat_prod::Q     = Q;

	//fill x vector
	x = vector<double>(n);
	for (int i = 0; i < n; i++) { x[i] = i * dx; }
}

// copy constructor
one_dim_diffusion_eq_heat_prod::one_dim_diffusion_eq_heat_prod(const one_dim_diffusion_eq_heat_prod & copy) {
	u0    = copy.u0;
	uL    = copy.uL;
	u     = copy.u;
	alpha = copy.alpha;
	n     = copy.n;
	Q     = copy.Q;
	t     = copy.t;
	dt    = copy.dt;
	x     = copy.x;
}

// update u using explicit Euler
void one_dim_diffusion_eq_heat_prod::explicit_euler_update() {
	vector<double> u_new(n);
	for (int i = 1; i < n - 1; i++) {
		u_new[i] = alpha * u[i - 1] + (1 - 2 * alpha) * u[i] + alpha * u[i + 1] + Q(x[i], t)*dt;
	}

	u_new[0]     = u0;
	u_new[n - 1] = uL;
	u            = u_new;
	t            += dt;
}