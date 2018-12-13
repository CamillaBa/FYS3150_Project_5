#pragma once
#include <iostream>

#include "vector.h"
#include "tridiag_solver.h"
#include "one_dim_diffusion_eq_heat_prod.h"
#include "physics.h"

//================================================================================
// Unit tests
//================================================================================

void unit_test_tridiagsolver() {
	/* Checking that the matrix equation
	2 -1  0 0     x        1
	-1 2 -1 0     y    =   5
	0 -1 2 -1     z        8
	0 0 -1  2     y        4

	is solved correctly by the function "sp_tri_ma_eq_so".
	The correct solution is (x,y,z,y)=(7.8,14.6,16.4,10.2),
	found using https://rrefcalculator.com/.
	*/

	vector<double> d(4);
	d[0] = 1;
	d[1] = 5;
	d[2] = 8;
	d[3] = 4;
	vector<double> solution = sp_tri_ma_eq_so(-1, 2, -1, d);

	if (solution[0] - 7.8 < 0.01 &&
		solution[1] - 14.6 < 0.01 &&
		solution[2] - 16.4 < 0.01 &&
		solution[3] - 10.2 < 0.01)
	{
		std::cout << "Tridiagonal matrix solver passes unit test." << std::endl;
	}
}

void unit_test_steady_state() {
	/* Checking that the 1-dimensional heat eq with heat
	prod converges to known analytical solution within certain epsilon.

	This serves as purpose to verify that both the steady state solution
	and the numerical solver is implemented correctly.
	We need the steady state solution to solve the post enrichment problem.
	*/


	double dx = 0.01;
	double dt = 0.01*dx*dx;
	int     n = 0.8 / dx + 1;

	// initial/boundary conditions
	vector<double> T_init(n);
	T_init[0] = 8. / 9000;
	T_init[n - 1] = 1300. / 9000;

	one_dim_diffusion_eq_heat_prod eq(T_init, dx, dt, Q_pre_radioactive_enrichment_1D);
	// update time to t corresponding to 1 billion years
	for (int j = 0; j*dt <= 1; j++) { eq.explicit_euler_update(); }
	vector<double> solution_t1 = eq.u;

	// store analytical steady state
	vector<double> steady_state(n);
	for (int i = 0; i < n; i++) {
		steady_state[i] = steady_state_solution(i*dx);
	}

	// find relative error
	double absolute_error = 0;
	double norm_steady_state = 0;
	for (int i = 0; i < n; i++) {
		absolute_error += (steady_state[i] - solution_t1[i])*(steady_state[i] - solution_t1[i]);
		norm_steady_state += steady_state[i] * steady_state[i];
	}
	absolute_error = sqrt(absolute_error);
	norm_steady_state = sqrt(norm_steady_state);

	// print relative error to terminal
	double rel_err = absolute_error / norm_steady_state;
	if (rel_err < 0.005) {
		std::cout << "Steady state unit test passed: 1 dimensional heat equation with heat production \nagrees with analytical steady state solution. Rel error:" << rel_err << std::endl;
	}
}