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


// check for errors
// the sum of all probabilities should be equal to 1
void Quantum::checkStatus(void) {
	assert((1 << _size) == data.size());
	buffer.resize(data.size());
	double prob = 0.0; 
	for(size_t state = 0u;state < data.size();state++) {
		prob += getProbability(state);
	}
	assert(abs(prob - 1.0) < 1e-8);
}
// return number of qubits
size_t Quantum::size(void) {
	return _size;
}

// Initializer
Quantum::Quantum(size_t __size) {
	if(__size == 0u) return;
	_size = __size;
	data.resize(1u << __size, 0.0);
	std::fill(data.begin(), data.end(), 0.0);
	data[0] = 1.0;
	checkStatus();
}
Quantum::Quantum(std::vector< std::complex<double> > _data) {
	_size = log(_data.size()) / log(2);
	assert((1 << _size) == _data.size());
	data = _data;
	checkStatus();
}
Quantum::Quantum(size_t __size, std::vector< std::pair< size_t, std::complex<double> > > _data) {
	Quantum(_size = __size);
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
void Quantum::addQubits(size_t __size) {
	assert(_size > 0u);
	data.resize(1 << (_size + __size), 0.0f);
	std::fill(std::next(data.begin(), (1 << _size) + 1), data.end(), 0);
	_size += __size;
	checkStatus();
}
void Quantum::addQubits(Quantum qubits) {
	std::vector< std::complex<double> > new_data(1u << (_size + qubits.size()));
	for(size_t _idx = 0u; _idx < data.size(); _idx++) {
		for(size_t __idx = 0u; __idx < qubits.data.size(); __idx++) {
			new_data[(__idx << data.size()) ^ _idx] = data[_idx] * qubits.data[__idx];
		}
	}
	_size += qubits.size();
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
	assert(0 <= idx and idx < (_size));
	std::fill(buffer.begin(), buffer.end(), 0.0);
	for(size_t _state = 0u;_state < buffer.size();_state++) {
		buffer[_state] += ((_state >> idx) & 1 ? -__quantum_sqrt_half : __quantum_sqrt_half) * data[_state];
		buffer[_state ^ (1u << idx)] += __quantum_sqrt_half * data[_state];
	}

	data = buffer;
	checkStatus();
}
// apply some Hadamard gates to some data
template<typename... Args> 
void Quantum::Hadamard(size_t idx, Args... args) {
	Hadamard(idx);
	Hadamard(args...);
}

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

// begin C-NOT and toffoli gates 
// apply C-NOT to every qubits in hash_val and last_idx
void Quantum::CnotHashedIdx(size_t hash_val, size_t last_idx) {
	assert(0 < hash_val and hash_val < data.size());
	assert(((hash_val >> last_idx) & 1) == 0);
	for(size_t _idx = 0;_idx < data.size();_idx++) {
		if((__builtin_popcount(_idx & hash_val) & 1) and (_idx >> last_idx & 1)) {
			std::swap(data[_idx], data[_idx ^ (1u << last_idx)]);
		}
	}
}
// apply C-NOT gates to all qubits in hash_val
template<typename... Args>
void Quantum::CnotHashedIdx(size_t hash_val, size_t idx, Args... args) {
	assert(0 <= idx and idx < _size);
	CnotHashedIdx(hash_val ^ (1 << idx), args...);
}
// apply C-NOT gates to some given qubits
void Quantum::Cnot(std::vector<size_t> idx_list, size_t last_idx) {
	size_t hash_val = 0u;
	for(size_t idx : idx_list) {
		assert(0 <= idx and idx < _size);
		assert(idx != last_idx);
		hash_val ^= 1u << idx;
	}
	assert(0 <= last_idx and last_idx < data.size());
	CnotHashedIdx(hash_val, last_idx);
}
template<typename... Args>
void Quantum::Cnot(size_t idx, Args... args) {
	assert(0 <= idx and idx < _size);
	CnotHashedIdx(1u << idx, args...);
}
void Quantum::Cnot(std::vector<size_t> idx_list) {
	assert(idx_list.size() > 0);
	size_t last_idx = idx_list.back();
	idx_list.pop_back();
	Cnot(idx_list, last_idx);
}
void Quantum::CnotRange(size_t left, size_t right, size_t last_idx) {
	assert(0 <= left and left <= right and right < _size);
	assert(last_idx < left or right < last_idx);
	assert(0 <= last_idx and last_idx < _size);
	size_t hash_val = (1u << (right + 1u)) - (1u << left);
	CnotHashedIdx(hash_val, last_idx);
}
template<typename... Args>
void Quantum::Toffoli(size_t idx, Args... args) {
	assert(0 <= idx and idx < _size);
	CnotHashedIdx(1u << idx, args...);
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
// end C-NOT and Toffoli gate

// swap 2 qubits' phases
// this function should be implement by using 3 c-not gates
void Quantum::swap(size_t idx1, size_t idx2) {
	assert(0 <= idx1 and idx1 < _size);
	assert(0 <= idx2 and idx2 < _size);
	if(idx1 == idx2) return;	
/*	
	// this is how it should be done	
	Cnot(idx1, idx2);
	Cnot(idx2, idx1);
	Cnot(idx1, idx2);
*/	
	for(size_t _idx = 0;_idx < data.size();_idx++) {
		if(((_idx >> idx1) & 1) == 0 and ((_idx >> idx2) & 1)) {
			std::swap(data[_idx], data[_idx ^ (1 << idx1) ^ (1 << idx2)]);
		}
	}
}

#endif