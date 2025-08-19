// QuantumCircuit implementation

#include "QuantumCircuit.hpp"
#include <cmath>
#include <iostream>
#include <bitset>

#define C std::complex<double>

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

void QuantumCircuit::NOT(int qi) {

	Matrix<C> NOT(2, 2, C(0.0, 0.0));
	NOT(0, 1) = C(1.0, 0.0);
	NOT(1, 0) = C(1.0, 0.0);

	Matrix I1 = Matrix<C>::identity(size_t(1) << qi);
	Matrix I2 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));

	Matrix fullGate = I1.tensorProduct(NOT).tensorProduct(I2);

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