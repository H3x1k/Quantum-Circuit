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

	/*
	
	QuantumCircuit qc(qn (number of qubits), ? cn (number of cbits))

	qc.<single qubit gate symbol>( qi (qubit index) or qi[] (multiple target qubits) )
	// H = Hadamard

	qc.<controlled gate symbol>( qi[] (indices of control qubits), qi (qubit index) or qi[] (multiple target qubits) )

	qc.<multi-qubit gate symbol>( qi[] (multiple target qubits) )

	qc.addGate() // multiple options


	*/

	/*

	Matrix qubitState1(2, 1, std::complex<double>(0.0, 0.0));
	qubitState1(0, 0) = { 1.0, 0.0 };
	Matrix qubitState2(2, 1, std::complex<double>(0.0, 0.0));
	qubitState2(0, 0) = { 1.0, 0.0 };

	Matrix H(2, 2, std::complex<double>(1.0, 0.0));
	H = H * (1 / sqrt(2));
	H(1, 1) = -(1 / sqrt(2));

	qubitState1.print();

	qubitState1 = H * qubitState1;

	qubitState1.print();

	double prob0 = mag(qubitState1(0, 0));
	double prob1 = mag(qubitState1(1, 0));
	
	std::cout << "Probability of measuring |0> : " << prob0 << std::endl;
	std::cout << "Probability of measuring |1> : " << prob1 << std::endl;

	qubitState1 = qubitState1.tensorProduct(qubitState1);

	qubitState1 = qubitState1.tensorProduct(qubitState1);

	qubitState1.print();

	*/

	int nq = 5;

	qcf::QuantumCircuit qc(nq);

	for (int i = 0; i < nq; i++) {
		qc.H(i);
	}

	qc.printProb();

	while (1);

	return 0;
}