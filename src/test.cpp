#include "QuantumCircuit.hpp"

int main() {
	int nq = 4;
	qcf::QuantumCircuit qc(nq);

	qc.H(0);
	qc.RX(2, 1.57);
	qc.CNOT(1, 3);

	qc.printProb();

	qcf::MeasurementBatch mb = qc.measure_batch({ 0 }, 1024);
	mb.print();

	qc.printDiagram();

	while (1);

	return 0;
}