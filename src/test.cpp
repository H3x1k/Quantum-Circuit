#include "QuantumCircuit.hpp"

// TODO:
// add display for multi-qubit gates and multi control qubits
// add time component for operations so that multiple gates can be displayed on one line

static Matrix<std::complex<double>> modular_mult_matrix(int a, int n, int q) {
	using C = std::complex<double>;
	int dim = 1 << q;
	Matrix<C> m(dim, dim, C(0.0, 0.0));
	for (int y = 0; y < dim; ++y) {
		int result = (y * a) % n;
		m(result, y) = C(1.0, 0.0);
	}
	return m;
}

int main() {
	const int n = 21; // doesnt work
	const int a = 2;
	const int t = 5; // number of counting qubits
	const int q = 5; // number of work qubits
	int nq = t + q;

	qcf::QuantumCircuit qc(nq);

	std::vector<size_t> counting;
	for (int i = 0; i < t; ++i)
		counting.push_back(i);
	qc.H(counting);

	qc.X({ t }); // set the first work qubit to |1>

	for (size_t i = 0; i < t; ++i) {
		int power = 1 << i;
		int a_exp = 1;
		for (int j = 0; j < power; ++j)
			a_exp = (a_exp * a) % n;
		auto modmult = modular_mult_matrix(a_exp, n, q);
		qcf::Gate modmult_gate(modmult, Index::range(t, t + q));
		qc.apply_controlled(modmult_gate, Index({ i }));
	}

	qc.IQFT(Index::range(0, t));

	for (int i = 0; i < 100; i++) {
		qcf::Measurement m = qc.measure_all(false);
		m.print();
	}

	qc.printState();
	qc.printProb();

	while (1);

	return 0;
}