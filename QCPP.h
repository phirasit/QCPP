#ifndef __QUANTUM_CPP_H__
#define __QUANTUM_CPP_H__

#include <algorithm>
#include <complex>
#include <cstdarg>
#include <functional>
#include <vector>

class Quantum {
    
    private:

    	size_t size;
    	std::vector< std::complex<double> > data;

    	void checkStatus(void);

    public:

    	void addQubits(size_t);
    	void addQubits(Quantum);

    	void Hadamard(size_t);
    	void Hadamard(size_t, size_t...);
    	void Hadamard(std::vector<size_t>);
    	
    	double get_probability(size_t);
    	size_t get_state(void);

    	Quantum(size_t);
    	Quantum(std::vector< std::complex<double> >);
    	Quantum(size_t, std::vector< std::pair< size_t, std::complex<double> > >);

};

#endif