/*
Deutsch-Jozsa Problem
=====================

Given a list of number(0 / 1) of length n which 
have one of the following condtions

	- all number are equal to each other
	- number of 0s is equal to number of 1s

The task is to determine which condition the list has.

By using conventional compution in worst case scenario. 
At least n/2+1 queries have to be made in order to ensure 
that the answer is correct.

By using randomized algorithm, more than half probability of
achieving the correct answer can be done by using only 2 queries.

However, in quantum computing, this problem can be solved by 
using only 1 quantum query which will be faster than other computation.

input :
first line : n (size of list) Note: this number should be power of 2
next n lines : ai (values in the list) Note: 0 <= ai <= 1 

*/

#include <QCPP.h>
#include <iostream>
#include <cmath>
#include <cassert>
#include <ctime>
#include <vector>

std::vector<int> arr;
int n;

int main() {

	srand(time(NULL));
	std::cin >> n;

	// calcuate number of required qubits
	int bit = log(n) / log(2);

	// n should be power of 2
	assert((1 << bit) == n);

	// initialize qubits
	Quantum qubits(bit + 1);

	// enumerate every possibilities with equal probability
	qubits.hadamardRange(0, bit);

	// this part should be oracle black box
	// set value according to algorithm
	double prob = 1.0f / sqrt(2.0f * n);
	arr.resize(n);
	for(int i = 0;i < n;i++) {
		int &val = arr[i];
		std::cin >> val;
		assert(val == 0 or val == 1);

		if(val == 0) {
			qubits.setPhase(2 * i, prob);
			qubits.setPhase(2 * i + 1, -prob);
		}else  {
			qubits.setPhase(2 * i, -prob);
			qubits.setPhase(2 * i + 1, prob);
		}
	}

	// pass qubits through hadamard gates again
	qubits.hadamardRange(1, bit);

	for(int i = 0;i < (1 << qubits.size());i++) {
		std::cout << "Probability of being " << i << " : " << qubits.getPhase(i) << std::endl;
	}

	// there are 2 possible results 
	// condition will check for states of qubits 1 to bit
	if((qubits.getState() & (1 << (bit+1) - (1 << 1))) == 0) {
		// all data are equal
		std::cout << "all numbers are equal" << std::endl;
	}else {
		// num of 0s and 1s are equal
		std::cout << "number of 0s is equal to number of 1s" << std::endl;
	}

	return 0;
}