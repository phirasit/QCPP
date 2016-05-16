/*
example of how to use QFT
*/

#include <iostream>
#include <cmath>
#include <QCPP.h>

using namespace std;

int main() {

	Quantum qubits(5);

	// direct implementation
	/*
	for(int i = 0;i < qubits.size();i++) {
		qubits.hadamard(i);
		for(int j = i+1;j < qubits.size();j++) {
			qubits.controlled(qubits.phaseShiftFunc(i, 2 * M_PI / pow(2, j-i+1)), j);
		}
	}
	
	for(int i = 0, j = qubits.size()-1;i < j;i++, j--) {
		qubits.swap(i, j);
	}
	*/

	qubits.QFT();
	
	for(int i = 0;i < (1 << qubits.size());i++) {
		cout << "Possibility of being " << i << " = " << qubits.getPhase(i) << std::endl;
	}

	return 0;
}