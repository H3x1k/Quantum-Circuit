// QuantumCircuit header

#include <vector>
#include <complex>

#include "Gate.hpp"

namespace qcf {

	class QuantumCircuit {
	public:

		size_t numQubits;
		std::vector<std::complex<double>> stateVector;

		QuantumCircuit(int numQubits);
		
		void addGate(const Gate& gate);
	};

}