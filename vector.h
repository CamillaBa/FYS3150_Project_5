#pragma once
//=============================================================================================
// simple vector class
//=============================================================================================

template <class T>
class vector {
private:
	// class variables
	T* entries;

public:
	int n;

	// default constructor
	vector() {};

	// constructor
	vector(int n);

	// copy constructor
	vector(vector const & copy);

	// destructor
	~vector() {delete[] entries;}

	// access i'th entry in vector
	T & operator [] (int i) {
		return entries[i];
	}

	// assignment operator
	vector & operator = (const vector & other)
	{
		if (this != &other) {
			delete[] entries;
			n = other.n;
			entries = new T[n]{};
			for (int i = 0; i < n; i++) { entries[i] = other.entries[i]; };
		}
		return *this;
	}
};

// constructor
template <class T>
vector<T>::vector(int n) {
	entries = new T[n]{};
	vector::n = n;
}

// copy constructor
template <class T>
vector<T>::vector(vector const & copy) {
	vector::n = copy.n;

	// create new entries
	entries = new T[n]{};
	for (int i = 0; i < n; i++) { entries[i] = copy.entries[i]; };
}