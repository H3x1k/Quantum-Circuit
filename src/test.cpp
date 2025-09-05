#include "QuantumCircuit.hpp"

int main() {
	int nq = 8;
	qcf::QuantumCircuit qc(nq);

	int r = 0.785398 * sqrt(1 << nq) + 1; // pi/4 * sqrt(N)   N=2^n

	std::vector<size_t> all;
	for (int i = 0; i < nq; i++)
		all.push_back(i);

	// uniform state
	qc.H(all);

	for (int i = 0; i < r; i++) {
		// oracle
		qc.stateVector(101, 0) *= -1.0;

		// diffusion
		qc.H(all);
		qc.X(all);
		qc.CZ(all);
		qc.X(all);
		qc.H(all);
	}

	//qc.RX(2, Angle::degrees(90));
	//qc.CZ(1, 2);
	//qc.CNOT(1, 3);

	qc.printProb();

	qcf::Measurement m = qc.measure_all();
	m.print();

	//qc.printDiagram();

	while (1);

	return 0;
}