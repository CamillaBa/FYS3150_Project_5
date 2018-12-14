#include "Header.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <string>
#include <omp.h>

int main() {
	unit_test_tridiagsolver();

	//============================================================================================
	// Problem 5 c)
	//============================================================================================
	
	{
		vector<double> dx_vector(2);
		dx_vector[0] = 0.1;
		dx_vector[1] = 0.01;

		for (int i = 0; i < 2; i++) {
			// set dx and dt values and initiate diffusion equation instance
			double dx = dx_vector[i];
			double dt = 0.01*dx*dx;
			int    n  = (int) 1 / dx + 1;
			vector<double> u_init(n); u_init[n - 1] = 1;

			one_dim_diffusion_eq eq1(u_init, dx, dt); // instance for explicit Euler
			one_dim_diffusion_eq eq2(u_init, dx, dt); // instance for implicit Euler
			one_dim_diffusion_eq eq3(u_init, dx, dt); // instance for Crank-Nicolson

			std::string filename1, filename2;
			std::ofstream myfile1, myfile2;

			// explicit Euler ======================================================

			// open files for explicit Euler
			filename1 = "explicit_euler_dx_"
				      + std::to_string(dx)
				      + "_t1.txt";
			myfile1.open(filename1);
			filename2 = "explicit_euler_dx_"
				      + std::to_string(dx)
				      + "_t2.txt";
			myfile2.open(filename2);

			// print u to file when t = 0.05
			for (int j = 0; j*dt <= 0.05; j++) { eq1.explicit_euler_update(); }
			print_vector_to_file(eq1.u, myfile1);

			// print u to file when t = 1
			for (int j = 0; j*dt <= 0.95; j++) { eq1.explicit_euler_update(); }
			print_vector_to_file(eq1.u, myfile2);

			// close files
			myfile1.close();
			myfile2.close();

			// implicit Euler ======================================================

			// open files for implicit Euler
			filename1 = "implicit_euler_dx_"
				      + std::to_string(dx)
				      + "_t1.txt";
			myfile1.open(filename1);
			filename2 = "implicit_euler_dx_"
				      + std::to_string(dx)
				      + "_t2.txt";
			myfile2.open(filename2);

			// print u to file when t = 0.05
			for (int j  = 0; j*dt <= 0.05; j++) { eq2.implicit_euler_update(); }
			print_vector_to_file(eq2.u, myfile1);

			// print u to file when t = 1
			for (int j = 0; j*dt <= 0.95; j++) { eq2.implicit_euler_update(); }
			print_vector_to_file(eq2.u, myfile2);

			// close files
			myfile1.close();
			myfile2.close();

			// Crank-Nicolson ======================================================

			// open files for Crank-Nicolson
			filename1 = "crank_nicolson_dx_"
				      + std::to_string(dx)
				      + "_t1.txt";
			myfile1.open(filename1);
			filename2 = "crank_nicolson_dx_"
				      + std::to_string(dx)
				      + "_t2.txt";
			myfile2.open(filename2);

			// print u to file when t = 0.05
			for (int j = 0; j*dt <= 0.05; j++) { eq3.crank_nicolson_update(); }
			print_vector_to_file(eq3.u, myfile1);

			// print u to file when t = 1
			for (int j = 0; j*dt <= 0.95; j++) { eq3.crank_nicolson_update(); }
			print_vector_to_file(eq3.u, myfile2);

			// close files
			myfile1.close();
			myfile2.close();
		}
	}

	//============================================================================================
	// Problem 5 f)
	//============================================================================================

	{
        omp_set_num_threads(3); // set the desired number of threads


		vector<double> h_vector(2);
		h_vector[0] = 0.1;
		h_vector[1] = 0.01;
		for (int j = 0; j < 2; j++) {
			std::ofstream myfile;
			int m, n;
			double h  = h_vector[j];
            double dt = 0.01*h*h;
			n = m = (int) 1.0 / h + 1;

			// set initial condition
			matrix<double> u_init(m, n);
			for (int i = 0; i < m; i++) {
				u_init[m - 1][i] = 1;
				u_init[i][0]     = i*h;
				u_init[i][n - 1] = i*h;
			}

			// initiate instace of two_dim_eq
			two_dim_diffusion_eq eq(u_init, h, dt);

			for (int i = 0; i*dt <= 0.05; i++) { eq.explicit_euler_update(); };
			myfile.open("two_dim_implicit_euler_h_"
				        + std::to_string(h)
				        + "_t1.txt");
			print_matrix_to_file(eq.u,myfile);
			myfile.close();

			for (int i = 0; i*dt <= 1; i++) { 
				eq.explicit_euler_update(); 
				if (i % 10000 == 0) {std::cout << "Completed iteration: "<< i << std::endl;}
			};
			myfile.open("two_dim_implicit_euler_h_"
				        + std::to_string(h)
				        + "_t2.txt");
			print_matrix_to_file(eq.u, myfile);
			myfile.close();
		}
    }

	//============================================================================================
	// Problem 5 g)
	//============================================================================================

	unit_test_steady_state();

	{
		// pre enrichment
		double dx = 0.01;
		double dt = 0.01*dx*dx;
		int     n = 0.8 / dx + 1;

		// initial/boundary conditions
		vector<double> T_init(n); 
		T_init[0]     = 8./9000;
		T_init[n - 1] = 1300./9000;

		one_dim_diffusion_eq_heat_prod eq(T_init, dx, dt, Q_pre_radioactive_enrichment_1D);

		std::string filename1, filename2;
		std::ofstream myfile1, myfile2;

		// explicit Euler =======================================================================

		// open files for explicit Euler
		filename1 = "pre_enr_explicit_euler_dx_"
			      + std::to_string(dx)
			      + "_t1.txt";
		myfile1.open(filename1);
		filename2 = "pre_enr_explicit_euler_dx_"
			      + std::to_string(dx)
			      + "_t2.txt";
		myfile2.open(filename2);

		// print u to file when t = 0.05
		for (int j = 0; j*dt <= 0.05; j++) { eq.explicit_euler_update(); }
		print_vector_to_file(eq.u, myfile1);

		// print u to file when t = 1
		for (int j = 0; j*dt <= 0.95; j++) { eq.explicit_euler_update(); }
		print_vector_to_file(eq.u, myfile2);

		// close files
		myfile1.close();
		myfile2.close();
	}

	{
		omp_set_num_threads(3); // set the desired number of threads

		// post enrichment
		double h  = 0.01;
		double dt = 0.01*h*h;
		int n     = 1.0 / h + 1;
		int m     = 0.8 / h + 1;

		// setup initial condition
		matrix<double> T_init(m, n);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				T_init[i][j] = steady_state_solution(i*h);
			}
		}

		{
			two_dim_diffusion_eq_heat_prod eq(T_init, h, dt, Q_post_radioactive_enrichment);
			eq.explicit_euler_update();

			// print matrix to file
			std::ofstream myfile;
			std::string filename;
			filename = "post_enr_explicit_euler_h_"
				     + std::to_string(h)
				     + "_after_one_dt.txt";
			myfile.open(filename);
			print_matrix_to_file(eq.u, myfile);
			myfile.close();
		}

		two_dim_diffusion_eq_heat_prod eq(T_init, h, dt, Q_post_radioactive_enrichment);
		for (int j = 0; j*dt <= 1.0; j++) { eq.explicit_euler_update(); }

		// print matrix to file
		std::ofstream myfile;
		std::string filename;
		filename = "post_enr_explicit_euler_h_"
			     + std::to_string(h)
			     + ".txt";
		myfile.open(filename);
		print_matrix_to_file(eq.u, myfile);
		myfile.close();
	}

	//============================================================================================
	// Success message
	//============================================================================================

	std::cout << "Success!" << std::endl;
	std::cin.get();
}