#pragma once
#include "vector.h"
#include "tridiag_solver.h"

//================================================================================
// One-dimensional diffusion equation
//================================================================================

class one_dim_diffusion_eq {
	// boundary conditions
	double u0;
	double uL;


	double alpha; // dt/ dx^2
	int    n;     // number of points in approximation
public:
	vector<double> u;

	// constructor
	one_dim_diffusion_eq(vector<double> u, double dx, double dt);

	// copy constructor
	one_dim_diffusion_eq(const one_dim_diffusion_eq & copy);

	// update u using explicit Euler
	void explicit_euler_update();

	// update u using the implicit Euler
	void implicit_euler_update();

	// update u using Crank-Nicolson
	void crank_nicolson_update();
};

// constructor
one_dim_diffusion_eq::one_dim_diffusion_eq(vector<double> u, double dx, double dt) {
	one_dim_diffusion_eq::u     = u;
	one_dim_diffusion_eq::n     = u.n;
	one_dim_diffusion_eq::u0    = u[0];
	one_dim_diffusion_eq::uL    = u[n - 1];
	one_dim_diffusion_eq::alpha = dt / dx / dx;
}

// copy constructor
one_dim_diffusion_eq::one_dim_diffusion_eq(const one_dim_diffusion_eq & copy) {
	u0    = copy.u0;
	uL    = copy.uL;
	u     = copy.u;
	alpha = copy.alpha;
	n     = copy.n;
}

// update u using explicit Euler
void one_dim_diffusion_eq::explicit_euler_update() {
	vector<double> u_new(n);
	for (int i = 1; i < n - 1; i++) {
		u_new[i] = alpha * u[i - 1] + (1 - 2 * alpha) * u[i] + alpha * u[i + 1];
	}
	u_new[0]     = u0;
	u_new[n - 1] = uL;
	u            = u_new;
}

// update u using the implicit Euler
void one_dim_diffusion_eq::implicit_euler_update() {
	double a, b, c;
	a = c = -alpha;
	b        = 1 + 2 * alpha;
	u        = sp_tri_ma_eq_so(a, b, c, u);
	u[0]     = u0;
	u[n - 1] = uL;
}

// update u using Crank-Nicolson
void one_dim_diffusion_eq::crank_nicolson_update() {
	vector<double> u_new(n);
	for (int i = 1; i < n; i++) {
		u_new[i] = alpha * u[i - 1] + (2 - 2 * alpha) * u[i] + alpha * u[i + 1];
	}
	u        = sp_tri_ma_eq_so(-alpha, 2 + 2 * alpha, -alpha, u_new);
	u[0]     = u0;
	u[n - 1] = uL;
}