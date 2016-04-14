#ifndef __QUANTUM_CPP_H__
#define __QUANTUM_CPP_H__

#include <algorithm>
#include <complex>
#include <cmath>
#include <cassert>
#include <functional>
#include <vector>

const std::complex<double> __quantum_sqrt_half = 1 / sqrt(2.0);

class Quantum {
    
    private:

    	size_t _size;
    	std::vector< std::complex<double> > data, buffer;

    	void checkStatus(void);

    	void CnotHashedIdx(size_t, size_t);
    	template<typename... Args> void CnotHashedIdx(size_t, size_t, Args...);

    public:

    	void addQubits(size_t);
    	void addQubits(Quantum);

    	void Hadamard(size_t);
    	template<typename... Args> void Hadamard(size_t, Args...);
    	void Hadamard(std::vector<size_t>);
    	void HadamardRange(size_t, size_t);

    	template<typename... Args> void Cnot(size_t, Args...);
    	void Cnot(std::vector<size_t>, size_t);
    	void Cnot(std::vector<size_t>);
    	void CnotRange(size_t, size_t, size_t);
    	template<typename... Args> void Toffoli(size_t, Args...);
    	void Toffoli(std::vector<size_t> , size_t);
    	void Toffoli(std::vector<size_t>);
    	void ToffoliRange(size_t, size_t, size_t);

    	void swap(size_t, size_t);

    	double getProbability(size_t);
    	size_t getState(void);
    	std::complex<double> getPhase(size_t);
    	size_t size(void);
    	
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

// apply C-NOT gates to all qubits in hash_val
template<typename... Args>
void Quantum::CnotHashedIdx(size_t hash_val, size_t idx, Args... args) {
	assert(0 <= idx and idx < _size);
	CnotHashedIdx(hash_val ^ (1 << idx), args...);
}
template<typename... Args>
void Quantum::Cnot(size_t idx, Args... args) {
	assert(0 <= idx and idx < _size);
	CnotHashedIdx(1u << idx, args...);
}
template<typename... Args>
void Quantum::Toffoli(size_t idx, Args... args) {
	assert(0 <= idx and idx < _size);
	CnotHashedIdx(1u << idx, args...);
}

#endif