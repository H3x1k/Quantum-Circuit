#include <vector>
#include <complex>
#include <cmath>
#include <iostream>

#include "QuantumCircuit.hpp"
#include "Matrix.hpp"
//#include "Gate.hpp"

double static mag(std::complex<double> z) {
	return z.real() * z.real() + z.imag() * z.imag();
}

int main() {

	Matrix qubitState(2, 1, std::complex<double>(0.0, 0.0));
	qubitState(0, 0) = { 1.0, 0.0 };

	Matrix H(2, 2, std::complex<double>(1.0, 0.0));
	H = H * (1 / sqrt(2));
	H(1, 1) = -(1 / sqrt(2));

	qubitState.print();

	qubitState = H * qubitState;

	qubitState.print();

	double prob0 = mag(qubitState(0, 0));
	double prob1 = mag(qubitState(1, 0));
	
	std::cout << "Probability of measuring |0> : " << prob0 << std::endl;
	std::cout << "Probability of measuring |1> : " << prob1 << std::endl;

	while (1);

	return 0;
}