#include "QuantumCircuit.hpp"

int main() {
	int nq = 2;
	qcf::QuantumCircuit qc(nq);

	qc.X(1);
	qc.QFT({ 0,1 });

	qc.printState();

	/*
	int r = 0.785398 * sqrt(1 << nq) + 1; // pi/4 * sqrt(N)   N=2^n

	std::vector<size_t> all;
	for (int i = 0; i < nq; i++)
		all.push_back(i);

	qc.H(all); // uniform state
	for (int i = 0; i < r; i++) {
		qc.stateVector(101, 0) *= -1.0; // oracle
		// diffusion
		qc.H(all);
		qc.X(all);
		qc.CZ(all);
		qc.X(all);
		qc.H(all);
	}

	qc.printProb();

	qcf::Measurement m = qc.measure_all();
	m.print();
	*/

	//qc.printDiagram();

	while (1);

	return 0;
}