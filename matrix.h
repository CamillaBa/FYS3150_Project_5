#pragma once

//=============================================================================================
// simple matrix class
//=============================================================================================

template <class T>
class matrix {
private:
	// class variables
	T** entries;

public:
	// shape
	int m;
	int n;

	// default constructor
	matrix() {};

	// constructor
	matrix(int m, int n);

	// copy constructor
	matrix(matrix const & copy);

	// destructor
	~matrix();

	// access i'th row in matrix
	T* & operator[] (int i) {
		return entries[i];
	}

	// assignment operator
	matrix & operator = (const matrix & other)
	{
		if (this != &other) {

			for (int i = 0; i < m; i++) {
				delete[] entries[i];
			}
			delete[] entries;

			m = other.m;
			n = other.n;

			entries = new T *[m];
			for (int i = 0; i < m; i++) {
				entries[i] = new T[n]();
			}

			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					entries[i][j] = other.entries[i][j];
				}
			}
		}
		return *this;
	}
};

// constructor
template <class T>
matrix<T>::matrix(int m, int n){
	entries = new T *[m];
	for (int i = 0; i < m; i++) {
		entries[i] = new T[n]();
	}
	matrix::m = m;
	matrix::n = n;
}


// copy constructor
template <class T>
matrix<T>::matrix(matrix const & copy) {
	matrix::m = copy.m;
	matrix::n = copy.n;

	// create new entries
	entries = new T *[m];
	for (int i = 0; i < m; i++) {
		entries[i] = new T[n]();
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			entries[i][j] = copy.entries[i][j];
		}
	}
}

// destructor
template <class T>
matrix<T>::~matrix() {
	for (int i = 0; i < m; i++) {
		delete[] entries[i];
	}
	delete[] entries;
}