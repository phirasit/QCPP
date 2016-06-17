/*
====
QCPP
====
created by : Phirasit Charoenchitseriwong

QCPP is a C++ library for quantum computer simulation
it can be used to simulate quantum algorithms by using 
a conventional method. Hence, the performance might not
be good. However, the purpose of this projects is to
be a tool that can show the correctness of quantum algorithms
such as Shor's algorithm and Grover's algorithm.

*/

// check for duplicate inclusion
#ifndef __QUANTUM_CPP_H__
#define __QUANTUM_CPP_H__

// list of header 
#include <algorithm>
#include <cassert>
#include <complex>
#include <cmath>
#include <functional>
#include <vector>

#ifndef __QUANTUM_CPP_CONSTS__
#define __QUANTUM_CPP_CONSTS__

const std :: complex < double > __quantum_sqrt_half__ = 1.0 / sqrt ( 2.0 ) ;

#endif

// function for controlled
#ifndef __QUANTUM_CPP_TYPEDEF__
#define __QUANTUM_CPP_TYPEDEF__

typedef std :: function < void ( size_t ) > ControlledFunc;

#endif

class Quantum {
    
    // -- begin helper function headers
    
    private:

    	// contians number of qubits
    	size_t _size;

    	// represent all possibility of every state 
    	std :: vector < std :: complex < double > > data, buffer;

    public:
    	// check for validity of state
    	void checkStatus ( void );
    
    // -- end helper function headers

	// -- begin add qubits headers
	// add more qubits
    
    public:

    	void addQubits ( size_t = 1u );

    	// add qubits from a quantum state
    	void addQubits( Quantum );

	// -- end add qubits

	// -- begin hadamard headers
	// hadamard function using with controlled function
	// should not be called directly but passed as parameter
    
    public:

    	// simple hadamard functions
    	void hadamard ( void );
    	void hadamard ( size_t );

    	// variadic hadamard function
    	template < typename... Args > void hadamard ( size_t , Args... );

    	// hadamard function for controlled gate
    	ControlledFunc hadamardFunc ( size_t );

    	// hadamard function with list of indices
    	void hadamard ( std :: vector < size_t > );

    	// variadic function with list of indices
    	template < typename... Args > void hadamard ( std :: vector < size_t > , Args... );

    	// hadamard function to every qubits in given range
    	void hadamardRange ( size_t , size_t );

	// -- end hadamard headers

	// -- begin controlled function headers
	// controlled function consider the state of qubits before applying function

    private:

    	// controlled function with given hashed indices
    	void controlled ( size_t , ControlledFunc );

    	// variadic controlled function with given hashed indice
    	void controlled ( size_t , ControlledFunc , size_t );
    	template < typename... Args > void controlled ( size_t , ControlledFunc , size_t , Args... );

    public:

    	// variadic controlled function without given hashed indice 
    	void controlled ( ControlledFunc , size_t );
    	template < typename... Args > void controlled (  ControlledFunc , size_t, Args... );

    	// simple controlled with list of controlled qubits
    	void controlled ( ControlledFunc , std :: vector < size_t > );

    // -- end controlled function headers
	// -- begin not headers
	// not function swap current state of qubits
	// function name has to be captial letter to avoid reversed word collision
    public:

    	// return not funtion for controlled function
    	ControlledFunc NotFunc ( size_t );

    	// simple not function
    	void Not ( size_t );

    	// apply not to every qubits
    	void Not ( void );

    	// variadic not function
    	template < typename... Args > void Not ( size_t , Args... );

	// -- end not headers

	// -- begin phaseShift headers
	// phase shift shifts amplitude of |1> state
	// note : phaseShift pi = phaseFlip
    public:

    	// phaseShift function for controlled
    	ControlledFunc phaseShiftFunc ( size_t , double );
    	
    	// simple phaseShift function
    	void phaseShift ( double );
    	void phaseShift ( size_t , double );

    	// simple phaseShift function with list of shifted qubits
    	void phaseShift ( std :: vector < size_t >, double );

    	// variadic phaseShift function
    	template < typename... Args > void phaseShift ( size_t, Args..., double );
	// -- end phaseShift headers

	// -- begin phaseFlip headers
	// phaseFlip inverses the phase of |1> state
	public:

    	// phaseShift function for controlled
    	ControlledFunc phaseFlipFunc ( size_t );

    	// simple phaseFlip function
    	void phaseFlip ( void );
    	void phaseFlip ( size_t );

    	// simple phaseFlip function with list of shifted qubits
    	void phaseFlip ( std :: vector < size_t > );

    	// variadic phaseFlip function
    	template < typename... Args > void phaseFlip ( size_t , Args... );

	// -- end phaseFlip headers

	// -- begin setPhase function headers
	// warning : this function should be used in oracle function only
	// this function might corrupt the system 
    public:

    	// simple setPhase function 
    	void setPhase ( size_t, std :: complex < double > );

    	// setPhase with list of queries
    	void setPhase ( std :: vector < std :: pair < size_t , std :: complex < double > > > );

    	// variadic setPhase function
    	template < typename... Args > void setPhase ( size_t , std :: complex < double > , Args... );

	// -- end setPhase headers

	// -- begin swap function headers
	// swap function swaps phase of two qubits
	public:
		// swap function for controlled function
		ControlledFunc swapFunc ( size_t , size_t );

		// simple swap function
    	void swap ( size_t , size_t );
	// -- end swap function headers

	// -- begin QFT function headers
    public:

    	// apply QFT to every qubits
    	void QFT ( void );

    	// using QFT with given index vector
    	void QFT( std :: vector < size_t > );

    	// apply QFT to every qubits in given range
    	void QFT( size_t , size_t );

	// -- end QFT function headers

	// -- begin helper function headers
    public:

    	// get probability of qubits collapsed to given state
    	double getProbability ( size_t );

    	// simulate quantum mechanic and output result
    	size_t getState ( void );


    	// return phase of given state
    	std :: complex < double > getPhase ( size_t );

    	// return number of qubits
    	size_t size ( void );

    // -- begin observe function headers

    private:
    	// observe with given hashed indices
    	size_t observeHashed ( size_t );

    public:
		// observe quantum mechanic
    	size_t observe ( void );

    	// observe only determined qubits
    	size_t observe ( std :: vector < size_t > );
    	size_t observeRange ( size_t , size_t );

    // -- end observe function headers

	// -- end helper function headers

	// -- begin constructor funciton headers
    public:

    	// simple constructor
    	Quantum ( size_t = 1u );

    	// constructor with given state of every possibility
    	// this should use for phenomenones like 'quantum entanglement'
    	Quantum ( std :: vector < std :: complex < double > > );

    	// constructor with some given possibilities
    	Quantum ( size_t ,  std :: vector < std :: pair < size_t, std :: complex < double > > > );
	// -- end constructor function headers

	// -- begin destructor function headers
  	public:

  		// simple destructor
    	~Quantum ( void );
	// -- end destructor function headers
};

// -- end list of headers

// -- begin helper functions
// check status function
// check for errors
//  - sum of all probabilities should be equal to 1.0
//  - each probability must be in range of 0.0 - 1.0
void Quantum::checkStatus ( void ) {
	assert ( ( 1u << _size ) == data.size() );
	buffer.resize ( data.size() );
	double prob = 0.0; 
	for ( size_t state = 0u ; state < ( 1u << _size ) ; state++ ) {
		double pstate = getProbability ( state );
		assert ( 0.0 <= pstate and pstate <= 1.0 );
		prob += pstate;
	}
	assert ( abs ( prob - 1.0 ) < 1e-8 );
}

// return number of qubits
size_t Quantum::size ( void ) {
	return _size;
}

// return probability of qubits collapse into state
double Quantum::getProbability ( size_t state ) {
	assert ( 0 <= state and state < data.size() );
	return pow ( abs ( data [ state ] ) , 2.0 );
}

// return phase of given state
std :: complex < double > Quantum::getPhase ( size_t state ) {
	return data [ state ] ;
}

// simulate the result if qubits are observed
size_t Quantum::getState ( void ) {
	double rnd = ( double ) rand () / RAND_MAX;
	for ( size_t state = 0u ; state < data.size() ; state++ ) {
		if( rnd <= getProbability ( state ) ) {
			return state;
		} else {
			rnd -= getProbability ( state );
		}
	}
	assert ( false );
}

// -- begin observe functions

// observe with hashed indices 
size_t Quantum::observeHashed ( size_t hashed ) {

	size_t result = getState(); 
	
	std :: fill ( buffer.begin() , buffer.end() , 0.0 );
	
	double sum_prob = 0.0;
	for ( size_t state = 0 ; state < ( 1u << size() ) ; state++ ) {
		if( ( state & hashed ) == ( result & hashed ) ) {
			sum_prob += getProbability ( state );
		}
	}

	double sqrt_prob = sqrt ( sum_prob );
	for ( size_t state = 0 ; state < ( 1u << size() ) ; state++ ) {
		if( ( state & hashed ) == ( result & hashed ) ) {
			data [ state ] /= sqrt_prob;
		} else {
			data [ state ] = 0.0;
		}
	}

	checkStatus(); 

	return result;
}

// observe the whole system
size_t Quantum::observe ( void ) {
	return observeRange(0, size() - 1u );
}
// observe with given list of indices
size_t Quantum::observe ( std :: vector < size_t > indices ) {
	size_t hashed = 0u;
	for ( size_t idx : indices ) {
		hashed ^= 1u << idx;
	}
	return observeHashed ( hashed );
}
// observe with given range of indices
size_t Quantum::observeRange ( size_t left, size_t right ) {
	size_t hashed = ( 1u << ( right + 1 ) ) - ( 1u << left );
 	return observeHashed ( hashed );
}

// -- end observe functions
// -- end helper functions

// -- begin Initializers / constructors
Quantum::Quantum ( size_t __size ) {
	if ( __size == 0u ) return;
	_size = __size;
	data.resize ( 1u << __size , 0.0 );
	std :: fill ( data.begin() , data.end() , 0.0 );
	data [ 0 ] = 1.0;
	checkStatus(); 
}
Quantum::Quantum( std :: vector < std :: complex < double > > _data ) {
	_size = log2(_data.size() );
	assert ( _size > 0 and ( 1u << _size ) == _data.size() );
	data.resize ( _data.size() );
	std :: copy ( _data.begin() , _data.end() , data.begin() );
	checkStatus(); 
}
Quantum::Quantum( size_t __size , std :: vector < std :: pair < size_t, std :: complex < double > > > _data ) : Quantum ( __size ) {
	std :: fill ( data.begin() , data.end() , 0.0 );
	for ( auto& state : _data ) {
		assert ( 0 <= state.first and state.first < data.size() );
		data [ state.first ] = state.second; 
	}
	checkStatus(); 
}
// -- end constructors

// -- begin destructor
Quantum::~Quantum ( void ) {}
// -- end destructor

// -- begin addQubits
// in case of needs for more qubits
// this functions should not be used (just in case of miscalculation)
void Quantum::addQubits ( size_t __size ) {
	assert ( _size > 0u );
	data.resize ( 1 << ( _size + __size ) , 0.0f );
	std :: fill ( std :: next ( data.begin() , ( 1 << _size ) + 1 ) , data.end() , 0.0 );
	_size += __size;

	// check for correctness
	checkStatus(); 
}
void Quantum::addQubits ( Quantum qubits ) {
	 std :: vector < std :: complex < double > > new_data(1u << (_size + qubits.size() ), 0.0);
	for ( size_t _idx = 0u ; _idx < data.size() ;  _idx++ ) {
		for ( size_t __idx = 0u ; __idx < qubits.data.size() ;  __idx++ ) {
			new_data [ ( __idx << data.size() ) ^ _idx ] = data [ _idx ] * qubits.data [ __idx ];
		}
	}
	_size += qubits.size(); 
	data.assign ( new_data.begin() , new_data.end() );

	// check for correctness
	checkStatus(); 
}
// -- end addQubits

// -- begin hadamard function
// apply hadamard gate to a qubit
void Quantum::hadamard ( size_t idx ) {

	assert ( 0 <= idx and idx < ( _size ) );

	std :: fill ( buffer.begin() , buffer.end() , 0.0 );
	for ( size_t state = 0u ; state < buffer.size() ; state++ ) {
		buffer [ state ] += ( ( state >> idx ) & 1 ? -1.0 : 1.0 ) * __quantum_sqrt_half__ * data [ state ];
		buffer [ state ^ ( 1u << idx ) ] += __quantum_sqrt_half__ * data [ state ];
	}
	std :: copy ( buffer.begin() , buffer.end() , data.begin() );

	// check for correctness
	checkStatus(); 
}

// apply some hadamard gates to some data
template < typename... Args > 
void Quantum::hadamard ( size_t idx , Args... args ) {
	hadamard ( idx );
	hadamard ( args... );
}

// apply hadamard to every qubits in indices
void Quantum::hadamard ( std :: vector < size_t > indices ) {
	for ( size_t idx : indices ) {
		hadamard ( idx );
	}
}

template < typename... Args >
void Quantum::hadamard ( std :: vector < size_t > indices , Args... args ) {
	for ( size_t index : indices ) {
		hadamard ( index );
	}
	hadamard ( args... );
}

// apply hadamard gates to every qubits in range [left, right]
void Quantum::hadamardRange ( size_t left, size_t right ) {
	for ( size_t _idx = left ; _idx <= right ; _idx++ ) {
		hadamard ( _idx );
	}
}

// apply hadamard to every qubits
void Quantum::hadamard ( void ) {
	hadamardRange ( 0, size() -1u );
}

// hadamard function for controlled gate
ControlledFunc Quantum::hadamardFunc ( size_t idx1 ) {
	return [ & , idx1 ] ( size_t idx ) {
		buffer [ idx ] += ( ( idx >> idx1 ) & 1 ? -1.0 : 1.0 ) * __quantum_sqrt_half__ * data [ idx ];
		buffer [ idx ^ ( 1u << idx1 ) ] += __quantum_sqrt_half__ * data [ idx ];
	};
}

// -- end hadamard function

// -- begin swap function
// swap 2 qubits' phases
// this function should be implemented by using 3 c-not gates
void Quantum::swap ( size_t idx1 , size_t idx2 ) {
	assert ( 0 <= idx1 and idx1 < size() );
	assert ( 0 <= idx2 and idx2 < size() );
	if ( idx1 == idx2 ) {
		return; 
	}

	// this is how it should be done	
	// controlled ( NotFunc ( idx1 ) , idx2 );
	// controlled ( NotFunc ( idx2 ) , idx1 );
	// controlled ( NotFunc ( idx1 ) , idx2 );
	
	// but, for a better performance
	std :: fill ( buffer.begin() , buffer.end() , 0.0 );
	for ( size_t _idx = 0 ; _idx < data.size() ; _idx++ ) {
		size_t bit1 = ( _idx >> idx1 ) bitand 1u;
		size_t bit2 = ( _idx >> idx2 ) bitand 1u;
		size_t nwidx = _idx xor ( bit1 << idx1 ) xor ( bit2 << idx2 );
		buffer [ nwidx xor ( bit1 << idx2 ) xor ( bit2 << idx1 ) ] = data [ _idx ];
	}
	std :: copy ( buffer.begin() , buffer.end() , data.begin() );
}

// swap function for controlled gate
ControlledFunc Quantum::swapFunc ( size_t idx1 , size_t idx2 ) {
	return [ = , idx1 , idx2 ] ( size_t idx ) {
		size_t bit1 = ( idx >> idx1 ) bitand 1u;
		size_t bit2 = ( idx >> idx2 ) bitand 1u;
		size_t nwidx = idx xor ( bit1 << idx1 ) xor ( bit2 << idx2 );
		buffer [ nwidx xor ( bit1 << idx2 ) xor ( bit2 << idx1 ) ] = data [ idx ];		
	};
}
// -- end swap function

// -- begin phase shift and phase flip functions 
void Quantum::phaseShift ( double ang ) {
	for ( size_t i = 0 ; i < size() ; i++ ) {
		phaseShift ( i , ang );
	}	
}
void Quantum::phaseShift ( size_t qubit_idx , double ang ) {

	assert ( 0 <= qubit_idx and qubit_idx < size() );

	std :: fill ( buffer.begin() , buffer.end() , 0 );
	for ( size_t _idx = 0 ; _idx < data.size() ; _idx++ ) {
		if( ( _idx >> qubit_idx ) & 1 ) {
			buffer [ _idx ] = data [ _idx ] * std :: exp ( std :: complex < double > ( 0 , 1 ) * ang );
		} else {
			buffer [ _idx ] = data [ _idx ];
		}
	}
	std :: copy ( buffer.begin() , buffer.end() , data.begin() );
}
void Quantum::phaseShift ( std :: vector < size_t > qubit_indices , double ang ) {
	for ( size_t qubit_idx : qubit_indices ) {
		phaseShift ( qubit_idx , ang );
	}
}
template < typename... Args > 
void Quantum::phaseShift ( size_t qubit_idx , Args... args , double ang ) {
	phaseShift ( qubit_idx , ang );
	phaseShift ( args... , ang );
}

// phaseShift for controlled function
ControlledFunc Quantum::phaseShiftFunc ( size_t idx , double ang ) {
	
	assert ( 0 <= idx and idx < size() );

	return [ = , idx , ang ] ( size_t index ) {
		if( ( index >> idx ) & 1 ) {
			buffer [ index ] = data [ index ] * std :: exp ( std :: complex < double > ( 0 , 1 ) * ang );
		} else {
			buffer [ index ] = data [ index ];
		}
	};
}

// phaseFlip is a simple form of phaseShift
void Quantum::phaseFlip ( void ) {
	phaseShift ( M_PI );
}
void Quantum::phaseFlip ( size_t qubit_idx ) {
	phaseShift ( qubit_idx , M_PI );	
}
void Quantum::phaseFlip ( std :: vector < size_t > qubit_indices ) {
	phaseShift ( qubit_indices , M_PI );
}
template < typename... Args > 
void Quantum::phaseFlip ( size_t qubit_idx , Args... args ) {
	phaseShift ( qubit_idx , args... , M_PI );
}

// -- end of rotate phase and flip phase functions

// -- begin setPhase
// this function is pretty dangerous
// setPhase will override phase of some state without checking
// this might violate some contrains in the begining
// please, ensure that before doing anything else the constrain is hold
void Quantum::setPhase ( size_t idx , std :: complex < double > value ) {
	assert ( 0 <= idx and idx < data.size() );
	data [ idx ] = value;
}
void Quantum::setPhase ( std :: vector < std :: pair < size_t, std :: complex < double > > > list ) {
	for ( auto it : list ) {
		setPhase ( it.first , it.second );
	}
}
// -- end setPhase

// -- begin controlled function
// controlled gate with hashed index
void Quantum::controlled ( size_t hashed , ControlledFunc func ) {

	assert ( 0 <= hashed and hashed < data.size() );

	std :: fill ( buffer.begin() , buffer.end() , 0 );
	for ( size_t index = 0 ; index < data.size() ; index++ ) {
		if ( ( index & hashed ) == hashed ) {
			func ( index );
		} else {
			buffer [ index ] += data [ index ];
		}
	}
	std :: copy ( buffer.begin() , buffer.end() , data.begin() );
}

// variadic version of controlled gate with hashed index
void Quantum::controlled ( size_t hashed , ControlledFunc func , size_t index ) {
	controlled ( hashed ^ ( 1u << index ) , func );
}
template < typename... Args > 
void Quantum::controlled ( size_t hashed , ControlledFunc func , size_t index, Args... args ) {
	assert ( 0 <= index and index < size() );
	controlled ( hashed ^ ( 1u << index ) , func , args... );
}


// variadic controlled gate functions
void Quantum::controlled ( ControlledFunc func , size_t idx ) {
	controlled ( 0 , func , idx );
}
template < typename... Args >
void Quantum::controlled ( ControlledFunc func , size_t idx , Args... args ) {
	controlled ( 0 , func , idx , args... );
}

// controlled gate with list of controlled qubits
void Quantum::controlled ( ControlledFunc func , std :: vector < size_t > index_list ) {

	size_t _idx = index_list.back(); 
	index_list.pop_back(); 

	size_t hashed = 0;
	for ( size_t index : index_list ) {
		assert ( 0 <= index and index < size() );
		hashed ^= 1u << index;
	}

	controlled ( hashed , func , _idx );
}
// -- end controlled function

// -- begin not function
// not function change qubit from |0> to |1>
void Quantum::Not ( size_t idx ) {

	assert ( 0 <= idx and idx < size() );

	std :: fill ( buffer.begin() , buffer.end() , 0 );
	for ( size_t index = 0 ; index < data.size() ; index++ ) {
		buffer [ index ^ ( 1u << idx ) ] = data [ index ];
	}

	std :: copy ( buffer.begin() , buffer.end() , data.begin() );
}

// variadic not function
template < typename... Args >
void Quantum::Not ( size_t idx , Args... args ) {
	Not ( idx );
	Not ( args... );
}

// not function for every qubits
void Quantum::Not ( void )  {
	for ( size_t i = 0 ; i < size() ; i++ ) {
		Not ( i );
	}
}
// -- end not function

// not function for controlled function
ControlledFunc Quantum::NotFunc ( size_t idx ) {

	assert ( ( 0 <= idx and idx < size() , " index out of bound " ) );
	
	return [ = , idx ] ( size_t index ) {
		buffer [ index ^ ( 1 << idx ) ] = data [ index ];
	};
}
// -- end not function

// -- begin QFT

// normal QFT (apply QFT to every qubits)
void Quantum::QFT ( void ) {
	
	for ( size_t i = 0 ; i < size() ; i++ ) {
		hadamard ( i );
		for ( size_t j = i + 1 ; j < size() ; j++ ) {
			controlled ( phaseShiftFunc ( i , 2.0f * M_PI / pow ( 2.0f , j-i+1 ) ) , j );
		}
	}
	
	for ( size_t i = 0 , j = size() - 1u ; i < j ; i++, j-- ) {
		swap ( i , j );
	}

}

// QFT with given index vector
void Quantum::QFT ( std :: vector < size_t > indices ) {

	for ( size_t id : indices ) {
		assert ( (0 <= id and id < size() , " index out of bound " ) );
	}
	
	for ( size_t i = 0 ; i < indices.size() ; i++ ) {
		hadamard ( indices[i] ) ;
		for ( size_t j = i+1 ; j < indices.size() ; j++ ) {
			controlled ( phaseShiftFunc ( indices[i] , 2.0f * M_PI / pow ( 2.0f , j-i+1 ) ) , indices[j] );
		}
	}

	for ( size_t i = 0 , j = indices.size() - 1u ; i < j ; i++ , j-- ) {
		swap ( indices[i] , indices[j] );
	}

}

// QFT with a given range of indices
void Quantum::QFT ( size_t l , size_t r ) {

	assert ( ( 0 <= l and l <= r and r < size() , " index out of bound " ) );

	for ( size_t i = l ; i <= r ; i++ ) {
		hadamard ( i );
		for ( size_t j = i + 1 ; j <= r ; j++ ) {
			controlled ( phaseShiftFunc ( i , 2.0f * M_PI / pow ( 2.0f , j-i+1 ) ) , j );
		}
	}

	for ( size_t i = l , j = r ; i < j ; i++ , j-- ) {
		swap ( i , j );
	}
}
// -- end QFT

#endif

// -- end QCPP.h 