#include <iostream>
#include "QCPP.h"
#include <ctime>

int main() {

	srand(time(NULL));

	Quantum qubits(5);

	// turns all qubits into every state with same probability
	qubits.Hadamard(0u, 1u, 2u, 3u, 4u);
	// or
	// for(int i = 0;i < 5;i++) qubits.Hadamard(i);
	// or
	// for(int i = 0;i < 5;i++) 


	for(int i = 0;i < 10;i++) {
		std::cout << "State of qubits : " << qubits.get_state() << std::endl;
	}

	return 0;
}