/*
example of using cnot gate in QCPP library
*/
#include <iostream>
#include <cmath>
#include <QCPP.h>

int main() {

	// init state
	double prob = 1 / sqrt(2.0 * 4);
	Quantum qubits(3, {{1, prob}, {3, prob}, {5, prob}, {7, prob}});

	std::cout << "Before apply cnot gate" << std::endl;
	for(int i = 0;i < (1 << qubits.size());i++) {
		std::cout << "Probability of being " << i << " : " << qubits.getProbability(i) << std::endl;
	}

	qubits.controlled([&](Quantum::__qubit_index__ idx1, size_t idx2) {
		qubits.Not(idx1, idx2);
	}, 2, 0);

	std::cout << "After apply cnot gate" << std::endl;
	for(int i = 0;i < (1 << qubits.size());i++) {
		std::cout << "Probability of being " << i << " : " << qubits.getProbability(i) << std::endl;
	}

	return 0;
}