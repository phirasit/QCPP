/*
Superdense Coding
=================

This is a new way to send information by using quantum mechanic.
The algorithm can transmit data of two conventional bits with 
the requirement of two qubits, but only one of them will be transmitted
through the channel (it saves lots of bandwidth). Combined this algorithm 
with no-cloning theorem and no-delete theorem, the method is very
effective in terms of package size and security.

*/

#include <iostream>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

#include <QCPP.h>

int main() {

	/*
	In the begining There are two qubits entangled to each other.
		= 1/sqrt(2) * (|00> + |11>)
	Both Alice and Bob will keep each qubit seperately.
	*/

	// create quantum entanglement (a bell pair)
	std::complex<double> _sqrt_half_ = 1.0 / sqrt(2.0);
	Quantum qubits({_sqrt_half_, 0.0f, 0.0f, _sqrt_half_});	

	/*
	When Alice try to send information of Bob,
	she will manipulate her qubit according to the text
	she want to send.
	*/

	// -- begin Alice's part
	int text;

	// check for validation of text
	do {
		std::cout << "input number [0-3] : ";
		std::cin >> text;
	}while(text < 0 or 3 < text);

	// Alice will manipulate only the first qubit
	if(text == 0) {
		// do nothing
	}else if(text == 1) {
		// apply not to first qubit
		qubits.Not(0);
	}else if(text == 2) {
		// apply phaseFlip
		qubits.phaseFlip(0);
	}else if(text == 3) {
		// apply Not and phaseFlip
		qubits.Not(0);
		qubits.phaseFlip(0);
	}else {
		assert(false);
	}

	// check probabilities -- the status of the system
	std::cout << " == After Alice == " << std::endl;
	for(int i = 0;i < (1 << qubits.size());i++) {
		std::cout << "Probability of being : " << i << " = " << qubits.getProbability(i) << qubits.getPhase(i) << std::endl;
	}

	// -- end Alice's part

	/*
	Alice send her qubit to Bob, and Bob will have both qubits in his hand.
	He will decrypt the text by apply c-not and hadamard respectively.
	*/

	// -- begin Bob's part

	// perform c-not
	qubits.controlled(qubits.NotFunc(1), 0);
	// perform hadamard to first qubit
	qubits.hadamard(0);

	// check probabilities
	std::cout << " == After Bob == " << std::endl;
	for(int i = 0;i < (1 << qubits.size());i++) {
		std::cout << "Probability of being : " << i << " = " << qubits.getProbability(i) << qubits.getPhase(i) << std::endl;
	}

	std::vector<int> table = {0, 2, 1, 3};

	// observe result
	int result = qubits.getState();
	std::cout << "Observed value = " << result << std::endl;
	// use look-up table to gain the actual message
	std::cout << "Real value = " << table[result] << std::endl;
	// -- end Bob's part

	// check for correctness
	assert(table[result] == text);
	
	return 0;

}