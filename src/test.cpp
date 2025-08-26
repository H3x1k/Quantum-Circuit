#include "QuantumCircuit.hpp"

int main() {
	int nq = 2;
	qcf::QuantumCircuit qc(nq);

	qc.H(0);
	qc.CNOT(0, 1);

	qc.printProb();

	qcf::MeasurementBatch mb = qc.measure_batch({ 0 }, 1024);
	mb.print();

	qc.printDiagram();

	while (1);

	return 0;
}