#ifndef __QUANTUM_CPP_H__
#define __QUANTUM_CPP_H__

#include <algorithm>
#include <complex>
#include <cmath>
#include <functional>
#include <vector>

const std::complex<double> __quantum_sqrt_half = 1 / sqrt(2.0);

class Quantum {
    
    private:


    	size_t size;
    	std::vector< std::complex<double> > data, buffer;

    	void checkStatus(void);

    public:

    	void addQubits(size_t);
    	void addQubits(Quantum);

    	void Hadamard(size_t);
    	template<typename... Args> void Hadamard(size_t, Args...);
    	void Hadamard(std::vector<size_t>);

    	double get_probability(size_t);
    	size_t get_state(void);
    	std::complex<double> get_phase(size_t);

    	Quantum(size_t = 1u);
    	Quantum(std::vector< std::complex<double> >);
    	Quantum(size_t, std::vector< std::pair< size_t, std::complex<double> > >);

    	~Quantum(void);
};

// some functions that i can't put in other files
// it's because 'template' that makes this code less beautiful

// apply some Hadamard gates to some data
template<typename... Args> 
void Quantum::Hadamard(size_t idx, Args... args) {
	Hadamard(idx);
	Hadamard(args...);
}

#endif