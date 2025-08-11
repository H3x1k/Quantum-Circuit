// QuantumCircuit header

#include <vector>
#include <complex>

#include "Gate.hpp"

namespace qcf {

	class QuantumCircuit {
	public:

		QuantumCircuit(int numQubits);
		
		void addGate(const Gate& gate);
	};

}