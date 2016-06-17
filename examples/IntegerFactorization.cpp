/*
=====================
Integer Factorization
=====================

Factor an integer by using Shor's algorithm.

*/

#include <iostream>
#include <algorithm>
#include <cmath>

#include <QCPP.h>

int main() {

	int N;

	const int LIMIT = 50;

	std :: cout << " Input your number [1-" << LIMIT << "] : ";
	std :: cin >> N;

	if(N > LIMIT) {
		std::cout << N << " is too large" << std::endl;
		return 0;
	}

	for ( int i = 2 ; i < N ; i++) {
		for( int val = i * i , cnt = 2 ; val <= N ; val *= i , cnt++ ) {
			if ( val == N ) {
				std :: cout << N << " can be expressed as ";
				for ( int j = 0 ; j < cnt ; j++ ) {
					std :: cout << i << "\nx"[ j + 1 < cnt ];
				}
				return 0;
			}
		}
	}

	if(N == 1) {
		std :: cout << " Why do you want to factorize 1 ? " << std :: endl;
		return 0;
	}

	if(N % 2 == 0) {
		int factor = 2;
		std :: cout << " Solution found : " << factor << " x " << N / factor << std :: endl;
		return 0;
	}

	int q = log2(2 * N * N);

	Quantum qubits(2 * q);

	std :: cerr << " Number of qubits = " << qubits.size() << std :: endl;

	// x should be chosen randomly, but the program will terminate too early
	// int x = rand() % (N-1) + 1;
	int x = N - 1;

	// apply QFT to the first half ( equivalent to hadamard in this case ) 
	qubits.QFT( 0, q-1 );

	// Quantum Oracle
	double prob = 1.0f / sqrt((double) (1 << q));
	for(int state = 0;state < (1 << qubits.size());state++) {
		qubits.setPhase(state, 0.0);
	}
	for(int state = 0, mul = 1;state < (1 << q);state++) {
		// change |0>|0> to |s>|f(s)>, where f(s) = (x ^ s) % N
		int idx = (mul << q) ^ state;
		qubits.setPhase(idx, prob);
		mul = (mul * x) % N;
	}

	// check for error ( should be done after using setPhase() )
	qubits.checkStatus();

	int result = qubits.observeRange( q, 2*q-1 );

	// mask the value 
	result >>= q;

	// std::cout << qubits.size() << std::endl;
	// for(int state = 0;state < (1 << qubits.size());state++) {
	// 	if(qubits.getProbability(state) > 0) {
	// 		std::cout << "Probability of being " << state << " = " << qubits.getProbability(state) << std::endl;
	// 	}
	// }

	// apply QFT to the first half again ( this time isn't equivalent to hadamard )
	qubits.QFT( 0, q-1 );

	// std::cout << result << std::endl;

	const int cnt_limit = 100;

	for ( int try_cnt = 0 ; try_cnt < cnt_limit ; try_cnt++ ) {

		// get result
		int c = qubits.getState();

		// mask the result
		c = c & ( (1 << q) - 1u);

		if ( c == 0 ) {
			std :: cerr << " Measured zero , try again. " << std :: endl;
			continue;
		}

		int w = 1 << q;

		std :: cerr << " Measured " << c << " (" << c / w << ")" << std :: endl;

		int _gcd = std :: __gcd ( c ,  w );

		c /= _gcd;
		w /= _gcd;

		std :: cerr << " Fractional Approximation is " << c << " / " << w <<  std :: endl;

		if ( ( w % 2 ) && ( 2 * w < ( 1 << q ) ) ) {
			std :: cerr << " Odd denominator , trying to expand by 2 . " << std :: endl;
			w *= 2;
		}

		if ( w % 2 ) {
			std :: cerr << " Odd period , try agian. " << std :: endl;
			continue;
		}

		std :: cout << " Possible Period : " << w << std :: endl;

		int a = ( ( int ) pow ( x , w / 2 ) % N + N + 1 ) % N;
		int b = ( ( int ) pow ( x , w / 2 ) % N + N - 1 ) % N;

		std :: cout << " Candidate coprime ( " << a << " , " << b << " ) " << std :: endl;

		a = std :: __gcd ( N , a );
		b = std :: __gcd ( N , b );

		int factor = std :: max ( a , b );

		if ( 1 < factor and factor < N ) {
			std :: cout << " Solution found : " << factor << " x " << N / factor << std :: endl;
			return 0;
		}
	}

	std :: cout << N << " is a prime " << std :: endl;

	return 0;
}