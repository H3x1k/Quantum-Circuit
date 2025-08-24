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

	int nq = 1;
	qcf::QuantumCircuit qc(nq);

	qc.H(0);

	qc.printProb();

	qcf::MeasurementBatch mb = qc.measure_batch({ 0 }, 1024);
	mb.print();

	while (1);

	return 0;
}