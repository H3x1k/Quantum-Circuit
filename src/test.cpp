#include "QuantumCircuit.hpp"

#include <iomanip>

// TODO:
// add display for multi-qubit gates and multi control qubits
// add time component for operations so that multiple gates can be displayed on one line

// add multi-qubit measurement to finish shor's algorithm
// find the highest occurring number then use continued fractions to approximate the computed fraction from c/2^w
// then use r to compute factors of N using gcd(a^r/2 +- 1, N)


static Matrix<std::complex<double>> modular_mult_matrix(int a, int n, int q) {
	using C = std::complex<double>;
	int dim = 1 << q;
	Matrix<C> m(dim, dim, C(0.0, 0.0));
	for (int y = 0; y < dim; ++y) {
		int result;
		if (y < n)
			result = (y * a) % n;
		else
			result = y;
		//m(y, result) = C(1.0, 0.0);
		m(result, y) = C(1.0, 0.0);
	}
	return m;
}

static void printCountingMarginals(const std::map<std::string, double>& fullDist, int t, int q) {
	std::map<int, double> marginals;

	for (auto& kv : fullDist) {
		const std::string& state = kv.first;   // e.g. "0000100000"
		double prob = kv.second;

		// take first t bits as counting register
		//std::string countBits = state.substr(0, t);
		std::string countBits = state.substr(t, t + q);
		int countVal = std::stoi(countBits, nullptr, 2);

		marginals[countVal] += prob;
	}

	std::cout << "Marginal distribution over counting register:\n";
	for (auto& kv : marginals) {
		std::cout << std::setw(2) << kv.first << " : " << kv.second << "\n";
	}
}


int main() {

	const int n = 21;
	const int a = 8;
	const int t = 6; // number of counting qubits
	const int q = 6; // number of work qubits
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
		qc.apply_controlled_test(modmult_gate, Index({ i }));
	}

	qc.IQFT(Index::range(0, t));

	//for (int i = 0; i < 100; i++) {
	//	qcf::Measurement m = qc.measure_all(false);
	//	m.print();
	//}

	//qc.printState();
	//qc.printProb();

	printCountingMarginals(qc.probabilityDistribution(), t, q);

	while (1);

	return 0;
}