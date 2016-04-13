#ifndef __QUANTUM_CPP_CPP__
#define __QUANTUM_CPP_CPP__

#include "QCPP.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

// check for errors
// the sum of all probabilities should be equal to 1
void Quantum::checkStatus(void) {
	assert((1 << size) == data.size());
	buffer.resize(data.size());
	double prob = 0.0; 
	for(size_t state = 0u;state < data.size();state++) {
		prob += get_probability(state);
	}
	assert(abs(prob - 1.0) < 1e-8);
}

// Initializer
Quantum::Quantum(size_t _size) {
	if(_size == 0u) return;
	size = _size;
	data.resize(1u << _size, 0.0);
	for(auto& state : data) {
		state = 0.0;
	}
	data[0] = 1.0;
	checkStatus();
}
Quantum::Quantum(std::vector< std::complex<double> > _data) {
	size = log(_data.size()) / log(2);
	assert((1 << size) == _data.size());
	data = _data;
	checkStatus();
}
Quantum::Quantum(size_t _size, std::vector< std::pair< size_t, std::complex<double> > > _data) {
	Quantum(size = _size);
	for(auto& state : _data) {
		data[state.first] = state.second; 
	}
	checkStatus();
}

// Deconstructor
Quantum::~Quantum(void) {
	data.clear();
}

// in case of needs for more qubits
// this functions should not be used (just in case of miscalculation)
void Quantum::addQubits(size_t _size) {
	assert(size > 0u);
	data.resize(1 << (size + _size), 0.0f);
	for(size_t new_state = (1 << size) + 1; new_state < (1 << (size + _size)); new_state++) {
		data[new_state] = 0.0f;
	}
	size += _size;
	checkStatus();
}
void Quantum::addQubits(Quantum qubits) {
	std::vector< std::complex<double> > new_data(1u << (size + qubits.size));
	for(size_t _idx = 0u; _idx < data.size(); _idx++) {
		for(size_t __idx = 0u; __idx < qubits.data.size(); __idx++) {
			new_data[(__idx << data.size()) ^ _idx] = data[_idx] * qubits.data[__idx];
		}
	}
	size += qubits.size;
	data = new_data;
	checkStatus();
}

// return probability of qubits collapse into _state
double Quantum::get_probability(size_t _state) {
	assert(0 <= _state and _state < data.size());
	return pow(abs(data[_state]), 2.0);
}

std::complex<double> Quantum::get_phase(size_t _state) {
	return 1.0;
}

// simulate the result if qubits are observed
size_t Quantum::get_state(void) {
	double rnd = (double) rand() / RAND_MAX;
	for(size_t state = 0u;state < data.size();state++) {
		if(rnd <= get_probability(state)) {
			return state;
		}else {
			rnd -= get_probability(state);
		}
	}
	assert(false);
}

// apply hadamard gate to a qubit
void Quantum::Hadamard(size_t idx) {
	assert(0 <= idx and idx < (size));
	for(size_t _state = 0u;_state < buffer.size();_state++) {
		buffer[_state] = 0.0;
	}
	for(size_t _state = 0u;_state < buffer.size();_state++) {
		buffer[_state] += __quantum_sqrt_half * data[_state];
		buffer[_state ^ (1u << idx)] += __quantum_sqrt_half * data[_state];
	}

	data = buffer;
	checkStatus();
}

// this function has been implemented in header file
// apply hadamard gates to some given qubits
/*
template<typename... Args> 
void Quantum::Hadamard(size_t idx, Args... args) {
	Hadamard(idx);
	Hadamard(args...);
}
*/

void Quantum::Hadamard(std::vector<size_t> idx_list) {
	for(size_t idx : idx_list) {
		Hadamard(idx);
	}
}

#endif