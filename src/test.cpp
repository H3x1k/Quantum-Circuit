#include "QuantumCircuit.hpp"

int main() {
	int nq = 2;
	qcf::QuantumCircuit qc(nq);

	// add index class
	// will use single, multi, and all index types
	// for easy use of applying gates to qubits

	qc.H({ {0, 1} });

	//qc.RX(2, Angle::degrees(90));
	//qc.CZ(1, 2);
	//qc.CNOT(1, 3);

	qc.printProb();

	qcf::MeasurementBatch mb = qc.measure_batch({ 0 }, 1024);
	mb.print();

	qc.printDiagram();

	while (1);

	return 0;
}