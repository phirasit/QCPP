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
		prob += getProbability(state);
	}
	assert(abs(prob - 1.0) < 1e-8);
}

// Initializer
Quantum::Quantum(size_t _size) {
	if(_size == 0u) return;
	size = _size;
	data.resize(1u << _size, 0.0);
	std::fill(data.begin(), data.end(), 0.0);
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
	std::fill(std::next(data.begin(), (1 << size) + 1), data.end(), 0);
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
double Quantum::getProbability(size_t _state) {
	assert(0 <= _state and _state < data.size());
	return pow(abs(data[_state]), 2.0);
}

std::complex<double> Quantum::getPhase(size_t _state) {
	return data[_state];
}

// simulate the result if qubits are observed
size_t Quantum::getState(void) {
	double rnd = (double) rand() / RAND_MAX;
	for(size_t state = 0u;state < data.size();state++) {
		if(rnd <= getProbability(state)) {
			return state;
		}else {
			rnd -= getProbability(state);
		}
	}
	assert(false);
}

// apply hadamard gate to a qubit
void Quantum::Hadamard(size_t idx) {
	assert(0 <= idx and idx < (size));
	std::fill(buffer.begin(), buffer.end(), 0.0);
	for(size_t _state = 0u;_state < buffer.size();_state++) {
		buffer[_state] += ((_state >> idx) & 1 ? -__quantum_sqrt_half : __quantum_sqrt_half) * data[_state];
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
// apply hadamard gates to every qubits in range [left, right]
void Quantum::HadamardRange(size_t left, size_t right) {
	for(size_t _idx = left;_idx <= right;_idx++) {
		Hadamard(_idx);
	}
}

// apply C-NOT to every qubits in hash_val and last_idx
void Quantum::Cnot(size_t hash_val, size_t last_idx) {
	assert(0 < hash_val and hash_val < data.size());
	for(size_t _idx = 0;_idx < data.size();_idx++) {
		if((__builtin_popcount(_idx & hash_val) & 1) and (_idx >> last_idx & 1)) {
			std::swap(data[_idx], data[_idx ^ (1u << last_idx)]);
		}
	}
}
// apply C-NOT to given qubits
void Quantum::Cnot(std::vector<size_t> idx_list, size_t last_idx) {
	size_t hash_val = 0u;
	for(size_t idx : idx_list) {
		assert(0 <= idx and idx < size);
		assert(idx != last_idx);
		hash_val ^= 1u << idx;
	}
	assert(0 <= last_idx and last_idx < data.size());
	Cnot(hash_val, last_idx);
}
void Quantum::Cnot(std::vector<size_t> idx_list) {
	assert(idx_list.size() > 0);
	size_t last_idx = idx_list.back();
	idx_list.pop_back();
	Cnot(idx_list, last_idx);
}
void Quantum::CnotRange(size_t left, size_t right, size_t last_idx) {
	assert(0 <= left and left <= right and right < size);
	assert(last_idx < left or right < last_idx);
	assert(0 <= last_idx and last_idx < size);
	size_t hash_val = (1u << (right + 1u)) - (1u << left);
	Cnot(hash_val, last_idx);
}
void Quantum::Toffoli(std::vector<size_t> idx_list, size_t last_idx) {
	Cnot(idx_list, last_idx);
}
void Quantum::Toffoli(std::vector<size_t> idx_list) {
	Cnot(idx_list);
}
void Quantum::ToffoliRange(size_t left, size_t right, size_t last_idx) {
	CnotRange(left, right, last_idx);
}

#endif