#include <iostream>
#include <QCPP.h>

int main() {

	Quantum qubits(5);

	qubits.hadamard(0, 1, 2);

	// check probabilities
	std::cout << " ==== Before Swaps ==== " << std::endl;
	for(int i = 0;i < (1 << qubits.size());i++) {
		std::cout << "Probability of being "  << i << " = " << qubits.getProbability(i) << std::endl;
	}

	// swap qubits
	qubits.swap(0, 3);
	qubits.swap(1, 4);

	// check again
	std::cout << " ==== After Swaps ==== " << std::endl;
	for(int i = 0;i < (1 << qubits.size());i++) {
		std::cout << "Probability of being "  << i << " = " << qubits.getProbability(i) << std::endl;
	}

	return 0;
}