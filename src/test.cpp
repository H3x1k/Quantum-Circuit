#include <vector>
#include <complex>
#include <cmath>
#include <iostream>

#include "QuantumCircuit.hpp"
#include "Matrix.hpp"
//#include "Gate.hpp"

int main() {
	int nq = 2;
	qcf::QuantumCircuit qc(nq);

	qc.H(0);
	qc.H(0);

	qc.printProb();

	qcf::MeasurementBatch mb = qc.measure_batch({ 0 }, 1024);
	mb.print();

	qc.printDiagram();

	while (1);

	return 0;
}