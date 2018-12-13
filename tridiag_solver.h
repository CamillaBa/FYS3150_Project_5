#pragma once
#include "vector.h"

//================================================================================
// Tridiagonal matrix equation solver
//================================================================================

vector<double> sp_tri_ma_eq_so(double a, double b, double c, vector<double> d) {
	// get length of b and d
	int            n = d.n;
	vector<double> beta(n);
	vector<double> delta(n);

	// row reduction
	beta[0] = b;
	delta[0] = d[0];
	for (int i = 1; i < n; i++) {
		double epsilon_iminus1 = a / beta[i - 1];
		beta[i] = b - epsilon_iminus1 * c;
		delta[i] = d[i] - epsilon_iminus1 * delta[i - 1];
	}

	// solving linear set of equations
	vector<double> v(n);
	v[n - 1] = delta[n - 1] / beta[n - 1];
	for (int i = n - 2; i >= 0; i--) {
		v[i] = (delta[i] - c * v[i + 1]) / beta[i];
	}

	return v;
}
