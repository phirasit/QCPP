#include <iostream>
#include <ctime>
#include <QCPP.h>

int main() {

	// make it more non-deterministic
	srand(time(NULL));

	// initialize Quantum state with 5 qubits
	Quantum qubits(5);

	// turns all qubits into every state with same probabilities
	qubits.hadamard(0, 1, 2, 3, 4);
	// or
	// for(int i = 0;i < qubits.size();i++) qubits.Hadamard(i);
	// or
	// qubits.hadamard({0, 1, 2, 3, 4});
	// would result the same

	// the possibility should be the same
	for(int i = 0;i < (1 << qubits.size());i++) {
		std::cout << "Probabilities of being " << i << " : " << qubits.getProbability(i) << std::endl;
	}

	// try to observe
	for(int i = 0;i < 10;i++) {
		std::cout << "State of qubits : " << qubits.getState() << std::endl;
	}

	return 0;
}