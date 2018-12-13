#pragma once
#include <cmath>

//================================================================================
// Physics
//================================================================================

// the various functions used for heat production

double Q_pre_radioactive_enrichment(double x, double y, double t) {
	if (x > 0.26666666666) { return 0.05; }
	if (x < 0.13333333333) { return 1.4; }
	else return 0.35;
}

double Q_pre_radioactive_enrichment_1D(double x, double t) {
	return Q_pre_radioactive_enrichment(x, 0, t);
}

double Q_post_radioactive_enrichment(double x, double y, double t) {
	if (x > 0.26666666666) {
		return 0.05 + 0.2*exp(-0.154889*t) + 0.2*exp(-0.049454*t) + 0.1*exp(-0.553884*t);
	}
	if (x < 0.13333333333) { return 1.4; }
	else return 0.35;
}

// the steady state solution for the pre-radioactive enrichment problem

double steady_state_solution(double x) {
	double Q1, Q2, Q3, c1, c2, d1, d2, f1, f2;

	Q1 = 1.4;         Q2 = 0.35;     Q3 = 0.05;
	c1 = -71. / 180;   c2 = -1. / 1125;
	d1 = -229. / 900;  d2 = -23. / 2250;
	f1 = -157. / 900;  f2 = -47. / 2250;

	double output;
	if (x > 4. / 15) {
		output = -Q3 / 2 * x*x - f1 * x - f2;
	}
	else if (x < 2. / 15) {
		output = -Q1 / 2 * x*x - c1 * x - c2;
	}
	else  output = -Q2 / 2 * x*x - d1 * x - d2;
	return output;
}