#pragma once
#include "vector.h"
#include "matrix.h"
#include <iostream>
#include <fstream>
#include <string>

//=============================================================================================

template<class T>
void print_vector_to_file(vector<T> &vec, std::ofstream &file) {
	/* Function to print the contents of an instance of vector<T>
	as a line in the given file, formated as CSV.
	*/
	int n = vec.n; // get length of vector
	for (int i = 0; i < n - 1; i++) {
		file << vec[i] << ',';
	}
	file << vec[n - 1] << "\n";
}

template<class T>
void print_vector_to_terminal(vector<T> &vec) {
	int n = vec.n; // get length of vector
	std::cout << "(";
	for (int i = 0; i < n - 1; i++) {
		std::cout << vec[i] << ',';
	}
	std::cout << vec[n - 1] << ")" << std::endl;
}

//=============================================================================================

template<class T>
void print_matrix_to_file(matrix<T> &mat, std::ofstream &file) {
	/* Function to print the contents of an instance of matrix<T>
	as a line in the given file, where each row is formated as a CSV.
	*/
	int m = mat.m;
	int n = mat.n;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n - 1; j++) {
			file << mat[i][j] << ",";
		}
		file << mat[i][n - 1];
		if (i < m - 1) {
			file << std::endl;
		}
	}
}

template<class T>
void print_matrix_to_terminal(matrix<T> &mat) {
	int m = mat.m; int n = mat.n; // get dimensions of matrix
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m - 1; j++) {
			std::cout << mat[i][j] << ", ";
		}
		std::cout << mat[i][n - 1] << std::endl;
	}
}