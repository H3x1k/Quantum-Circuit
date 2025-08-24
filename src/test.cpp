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

	int nq = 2;
	qcf::QuantumCircuit qc(nq);

	qc.H(0);
	qc.CNOT(0, 1);

	qc.printProb();

	qcf::Measurement m1 = qc.measure(0);
	m1.print();
	std::cout << std::endl;

	qc.printProb();

	qcf::Measurement m2 = qc.measure(1);
	m2.print();

	while (1);

	return 0;
}