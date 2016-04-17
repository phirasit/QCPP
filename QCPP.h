#ifndef __QUANTUM_CPP_H__
#define __QUANTUM_CPP_H__

#include <algorithm>
#include <cassert>
#include <complex>
#include <cmath>
#include <functional>
#include <vector>

const std::complex<double> __quantum_sqrt_half__ = 1 / sqrt(2.0);

class Quantum {
    
    private:

    	class __qubit_index__ {
	    	public:
	    		size_t index;
	    		__qubit_index__(size_t _index) : index(_index) {}
    	};

    	size_t _size;
    	std::vector< std::complex<double> > data, buffer;

    	void checkStatus(void);

    	// -- start add qubits header
    	// add more qubits
    public:
    	void addQubits(size_t = 1);

    	// add qubits from a quantum state
    	void addQubits(Quantum);
    	// -- end add qubits

    	// -- start hadamard header
    	// hadamard function using with controlled function
    	// should not be called directly but passed as parameter
    public:
    	void hadamard(__qubit_index__, size_t);

    	// simple hadamard function
    	void hadamard(size_t);

    	// variadic hadamard function
    	template<typename... Args> void hadamard(size_t, Args...);

    	// hadamard function with list of indexed
    	void hadamard(std::vector<size_t>);

    	// hadamard function to every qubits in given range
    	void hadamardRange(size_t, size_t);
    	// -- end hadamard header

    	// start controlled function 
    	// controlled function consider the state of qubits before applying function
    private:
    	// controlled function with given hashed indice
    	template<typename... Args> void controlled(size_t, std::function<void(__qubit_index__, size_t, Args...)>);

    	// variadic controlled functions with given hashed indice
    	template<typename... Args> void controlled(size_t, std::function<void(__qubit_index__, size_t, Args...)>, size_t);
    	template<typename... Args1, typename... Args2> void controlled(size_t, std::function<void(__qubit_index__, size_t, Args1...)>, size_t, Args2...);

    public:
    	// simple controlled (variadic function)
    	template<typename... Args1, typename... Args2> void controlled(std::function<void(__qubit_index__, size_t, Args1...)>, size_t, Args2...);

    	// simple controlled with list of controlled qubits
    	template<typename... Args> void controlled(std::function<void(__qubit_index__, size_t, Args...)>, std::vector<size_t>);

    	// -- start not header
    	// not function swap current state of qubits
    	// function name has to be captial letter to avoid reversed word collision
    public:

    	// not funtion for controlled function
    	// void not(__qubit_index__, size_t);

    	// simple not function
    	void Not(size_t);

    	// -- end not header

    	// -- start phaseShift header
    	// phase shift shifts amplitude of |1> state
    	// note : phaseShift pi = phaseFlip
    public:
      	// phaseShift function that should be using with controlled 
    	void phaseShift(__qubit_index__, size_t, double = M_PI);

    	// simple phaseShift function
    	void phaseShift(size_t, double = M_PI);

    	// simple phaseShift function with list of shifted qubits
    	void phaseShift(std::vector<size_t>, double = M_PI);

    	// variadic phaseShift function
    	template<typename... Args> void phaseShift(size_t, Args..., double = M_PI);
    	// -- end phaseShift header

    	// -- start phaseFlip header
    	// phaseFlip inverses the phase of |1> state
    public:

    	// phaseFlip function is an alias of phaseShift
    	template<typename... Args> auto phaseFlip(Args&&... args) -> decltype(phaseShift(std::forward<Args>(args)...)) {
    		return phaseShift(std::forward<Args>(args)...);
    	}

    	// -- end phaseFlip header

    	// -- setPhase function header
    	// warning : this function should be used in oracle function only
    	// this function might corrupt the system 
    public:

    	// simple setPhase function 
    	void setPhase(size_t, std::complex<double>);

    	// setPhase with list of queries
    	void setPhase(std::vector< std::pair< size_t, std::complex<double> > >);

    	// variadic setPhase function
    	template<typename... Args> void setPhase(size_t, std::complex<double>, Args...);

    	// -- end setPhase header

    	// -- start swap function header
    	// swap function swaps phase of two qubits
	public:
		// swap function for controlled function
		// template<typename... Args> void swap(std::function<__qubit_index__, Args...>); -undone

		// simple swap function
    	void swap(size_t, size_t);
    	// -- end swap function header

    	// -- start response function header
    public:

    	// get probability of start collapsed to given state
    	double getProbability(size_t);

    	// simulate quantum machenic and output result
    	size_t getState(void);

    	// return phase of given state
    	std::complex<double> getPhase(size_t);

    	// return number of qubits
    	size_t size(void);

    	// -- end response function header

    	// -- start constructor funciton header
    public:

    	// simple constructor
    	Quantum(size_t = 1u);

    	// constructor with given state of every possibility
    	// this should use for phenomenones like 'quantum entanglement'
    	Quantum(std::vector< std::complex<double> >);

    	// constructor with some given state of evey possibility
    	Quantum(size_t, std::vector< std::pair< size_t, std::complex<double> > >);
  		// -- end constructor function header

    	// -- start destructor function header
  	public:

  		// simple destructor
    	~Quantum(void);
    	// -- end destructor function
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

// Initializers / constructors
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
Quantum::Quantum(size_t __size, std::vector< std::pair< size_t, std::complex<double> > > _data) : Quantum(__size) {

	data[0] = 0.0;
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

// return phase of given state
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

// -- start hadamard function
// apply hadamard gate to a qubit
void Quantum::hadamard(size_t idx) {
	assert(0 <= idx and idx < (_size));
	std::fill(buffer.begin(), buffer.end(), 0.0);
	for(size_t _state = 0u;_state < buffer.size();_state++) {
		buffer[_state] += ((_state >> idx) & 1 ? -__quantum_sqrt_half__ : __quantum_sqrt_half__) * data[_state];
		buffer[_state ^ (1u << idx)] += __quantum_sqrt_half__ * data[_state];
	}

	data = buffer;
	checkStatus();
}
// apply some hadamard gates to some data
template<typename... Args> 
void Quantum::hadamard(size_t idx, Args... args) {
	hadamard(idx);
	hadamard(args...);
}

void Quantum::hadamard(std::vector<size_t> idx_list) {
	for(size_t idx : idx_list) {
		hadamard(idx);
	}
}
// apply hadamard gates to every qubits in range [left, right]
void Quantum::hadamardRange(size_t left, size_t right) {
	for(size_t _idx = left;_idx <= right;_idx++) {
		hadamard(_idx);
	}
}
// -- end hadamard function

// -- start swap function
// swap 2 qubits' phases
// this function should be implement by using 3 c-not gates
void Quantum::swap(size_t idx1, size_t idx2) {
	assert(0 <= idx1 and idx1 < _size);
	assert(0 <= idx2 and idx2 < _size);
	if(idx1 == idx2) return;	
/*	
	// this is how it should be done	
	cnot(idx1, idx2);
	cnot(idx2, idx1);
	cnot(idx1, idx2);
*/	
	for(size_t _idx = 0;_idx < data.size();_idx++) {
		if(((_idx >> idx1) & 1) == 0 and ((_idx >> idx2) & 1)) {
			std::swap(data[_idx], data[_idx ^ (1 << idx1) ^ (1 << idx2)]);
		}
	}
}
// -- end swap function

// -- start phase shift and phase flip functions 
void Quantum::phaseShift(size_t qubit_idx, double ang) {
	assert(0 <= qubit_idx and qubit_idx < _size);
	for(size_t _idx = 0;_idx < data.size();_idx++) {
		if((_idx >> qubit_idx) & 1) {
			data[_idx] *= std::exp((0, 1) * ang);
		}
	}
}
void Quantum::phaseShift(std::vector<size_t> qubit_idx_list, double ang) {
	for(size_t qubit_idx : qubit_idx_list) {
		phaseShift(qubit_idx, ang);
	}
}
template<typename... Args> 
void Quantum::phaseShift(size_t qubit_idx, Args... args, double ang) {
	phaseShift(qubit_idx, ang);
	phaseShift(args..., ang);
}
// -- end of rotate phase and flip phase functions

// -- start setPhase
// this function is pretty dangerous
// setPhase will override phase of some state without checking
// this might violate some contrains in the begining
// please, ensure that before doing anything else the constrain is hold
void Quantum::setPhase(size_t idx, std::complex<double> value) {
	assert(0 <= idx and idx < data.size());
	data[idx] = value;
}
void Quantum::setPhase(std::vector< std::pair< size_t, std::complex<double> > > list) {
	for(auto it : list) {
		setPhase(it.first, it.second);
	}
}
// -- end setPhase

// -- start controlled function
// controlled gate with hashed index
template<typename... Args> 
void Quantum::controlled(size_t hashed, std::function<void(__qubit_index__, size_t, Args...)> func) {
	assert(0 <= hashed and hashed < data.size());
	for(size_t _index = 0; _index < data.size();_index++) {
		if((_index & hashed) == hashed) {
			func(__qubit_index__(_index));
		}
	}
}
template<typename... Args> 
void Quantum::controlled(size_t hashed, std::function<void(__qubit_index__, size_t, Args...)> func, size_t _index) {
	assert(0 <= _index and _index < _size);
	controlled(hashed ^ (1 << _index), func);
}

// variadic version of controlled gate with hashed index
template<typename... Args1, typename... Args2> 
void Quantum::controlled(size_t hashed, std::function<void(__qubit_index__, size_t, Args1...)> func, size_t _index, Args2... args) {
	assert(0 <= _index and _index < _size);
	controlled(hashed ^ (1 << _index), func, args...);
}


// controlled gate
template<typename... Args1, typename... Args2>
void Quantum::controlled(std::function<void(__qubit_index__, size_t, Args1...)> func, size_t _index, Args2... args) {
	controlled(0, func, _index, args...);
}

// controlled gate with list of controlled qubits
template<typename... Args> 
void Quantum::controlled(std::function<void(__qubit_index__, size_t, Args...)> func, std::vector<size_t> index_list) {

	size_t hashed = 0;
	for(size_t _index : index_list) {
		assert(0 <= _index and _index < _size);
		hashed ^= 1u << _index;
	}

	controlled(hashed, func);
}
// -- end controlled function

// -- start not function
// not function change qubit from |0> to |1>
void Quantum::Not(size_t idx) {
	assert(0 <= idx and idx < _size);
	for(size_t _index = 0;_index < data.size();_index++) {
		if((_index >> idx) & 1) {
			std::swap(data[_index], data[_index ^ (1 << idx)]);
		} 
	}
}
// -- end not function

#endif