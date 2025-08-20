// QuantumCircuit implementation

#include "QuantumCircuit.hpp"
#include <cmath>
#include <iostream>
#include <bitset>

#define C std::complex<double>

#define E    2.7182818f
#define PI   3.1415926f
#define PI_4 0.7853981f

using namespace qcf;

QuantumCircuit::QuantumCircuit(int numQubits)
	: numQubits(numQubits),
	stateVector((size_t(1) << numQubits), 1, C(0.0, 0.0)) {
	stateVector(0, 0) = C(1.0, 0.0);
}



void QuantumCircuit::H(int qi) {

	Matrix<C> H(2, 2, C(1.0, 0.0));
	double invSqrt2 = 1.0 / sqrt(2.0);
	H = H * invSqrt2;
	H(1, 1) = -invSqrt2;

	Matrix I1 = Matrix<C>::identity(size_t(1) << qi);
	Matrix I2 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));

	Matrix fullGate = I1.tensorProduct(H).tensorProduct(I2);

	stateVector = fullGate * stateVector;
}

void QuantumCircuit::X(int qi) {

	Matrix<C> X(2, 2, C(0.0, 0.0));
	X(0, 1) = C(1.0, 0.0);
	X(1, 0) = C(1.0, 0.0);

	Matrix I1 = Matrix<C>::identity(size_t(1) << qi);
	Matrix I2 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));

	Matrix fullGate = I1.tensorProduct(X).tensorProduct(I2);

	stateVector = fullGate * stateVector;
}

void QuantumCircuit::Y(int qi) {
	Matrix<C> Y(2, 2, C(0.0, 0.0));
	Y(0, 1) = C(0.0, 1.0);
	Y(1, 0) = C(0.0, 1.0);

	Matrix I1 = Matrix<C>::identity(size_t(1) << qi);
	Matrix I2 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));

	Matrix fullGate = I1.tensorProduct(Y).tensorProduct(I2);

	stateVector = fullGate * stateVector;
}

void QuantumCircuit::Z(int qi) {
	Matrix<C> Z(2, 2, C(0.0, 0.0));
	Z(0, 0) = C(1.0, 0.0);
	Z(1, 1) = C(-1.0, 0.0);

	Matrix I1 = Matrix<C>::identity(size_t(1) << qi);
	Matrix I2 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));

	Matrix fullGate = I1.tensorProduct(Z).tensorProduct(I2);

	stateVector = fullGate * stateVector;
}


void QuantumCircuit::S(int qi) {
	Matrix<C> S(2, 2, C(0.0, 0.0));
	S(0, 0) = C(1.0, 0.0);
	S(1, 1) = C(0.0, 1.0);

	Matrix I1 = Matrix<C>::identity(size_t(1) << qi);
	Matrix I2 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));

	Matrix fullGate = I1.tensorProduct(S).tensorProduct(I2);

	stateVector = fullGate * stateVector;
}

void QuantumCircuit::T(int qi) {
	Matrix<C> T(2, 2, C(0.0, 0.0));
	T(0, 0) = C(1.0, 0.0);
	T(1, 1) = std::exp(C(0.0, PI_4));

	Matrix I1 = Matrix<C>::identity(size_t(1) << qi);
	Matrix I2 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));

	Matrix fullGate = I1.tensorProduct(T).tensorProduct(I2);

	stateVector = fullGate * stateVector;
}



void QuantumCircuit::printState() const {
	stateVector.print();
}

void QuantumCircuit::printProb() const {
	for (size_t i = 0; i < stateVector.rows; i++) {
		if (stateVector(i, 0) != C(0.0, 0.0)) {
			double prob = std::norm(stateVector(i, 0));
			std::cout << "|";
			for (size_t j = 0; j < numQubits; j++)
				std::cout << ((i >> (numQubits - 1 - j)) & 1);
			std::cout << "> : " << prob << std::endl;
		}
	}
}