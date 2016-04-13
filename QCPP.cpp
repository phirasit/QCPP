#ifndef __QUANTUM_CPP_CPP__
#define __QUANTUM_CPP_CPP__

#include "QCPP.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

void Quantum::checkStatus(void) {
	assert((1 << size) == data.size());
	double prob = 0.0; 
	for(int state = 0;state < data.size();state++) {
		prob += get_probability(state);
	}
	assert(abs(prob - 1.0) < 1e-8);
}

void Quantum::addQubits(size_t _size) {
	data.resize(size + _size, 0.0f);
	for(size_t new_state = (1 << size) + 1; new_state < (1 << (size + _size)); new_state++) {
		data[new_state] = 0.0f;
	}
	size += _size;
}

double Quantum::get_probability(size_t _state) {
	assert(0 <= _state and _state < data.size());
	return pow(abs(data[_state]), 2.0);
}

size_t Quantum::get_state(void) {
	double rnd = (double) rand() / RAND_MAX;
	for(size_t state = 0;state < data.size();state++) {
		if(rnd <= get_probability(state)) {
			return state;
		}else {
			rnd -= get_probability(state);
		}
	}
	assert(false);
}

void Quantum::addQubits(Quantum qubits) {
	std::vector< std::complex<double> > new_data(size + qubits.size);
	for(size_t _idx = 0; _idx < data.size(); _idx++) {
		for(size_t __idx = 0; __idx < qubits.data.size(); __idx++) {
			new_data[(__idx << data.size()) ^ _idx] = data[_idx] * qubits.data[__idx];
		}
	}
	size += qubits.size;
	data = new_data;
	checkStatus();
}

Quantum::Quantum(size_t _size = 0) {
	size = _size;
	data.resize(_size, 0.0);
	for(auto& state : data) {
		state = 0.0;
	}
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

void Quantum::Hadamard(size_t) {
	
}
#endif