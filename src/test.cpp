#include <vector>
#include <complex>
#include <cmath>

#include "QuantumCircuit.hpp"
#include "Gate.hpp"

int main() {

	std::vector<std::vector<std::complex<double>>> HMat = {
		{1 / sqrt(2), 1 / sqrt(2)},
		{1 / sqrt(2), -1 / sqrt(2)},
	};

	qcf::Gate H(HMat, 0);

	return 0;
}