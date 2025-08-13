// QuantumCircuit header

#include <vector>
#include <complex>

#include "Gate.hpp"

namespace qcf {

	class QuantumCircuit {
	public:

		size_t numQubits;
		std::vector<std::complex<double>> stateVector;
		std::vector<Gate> operations;

		QuantumCircuit(int numQubits);
		
		void addGate(const Gate& gate);
	};

}