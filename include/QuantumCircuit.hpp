// QuantumCircuit header

#include <vector>
#include <complex>
#include <map>

#include "Matrix.hpp"

namespace qcf {

	struct Measurement {
		std::vector<uint8_t> bits{};
		double probability = 0.0;

		void print() const;
	};

	struct MeasurementBatch {
		std::vector<size_t> qubits{};
		std::map<std::string, size_t> counts{};
	};

	class QuantumCircuit {
	public:

		size_t numQubits;
		Matrix<std::complex<double>> stateVector;

		QuantumCircuit(int numQubits);
		
		// Hadamard Gate
		void H(int qi);
		// Pauli Gates
		void X(int qi);
		void Y(int qi);
		void Z(int qi);
		// Phase Gates
		void S(int qi); void Sdag(int qi);
		void T(int qi); void Tdag(int qi);
		// Rotation Gates
		void RX(int qi, double angle);
		void RY(int qi, double angle);
		void RZ(int qi, double angle);
		// Controlled Gates
		void CNOT(int ci, int ti);

		// Measurement
		Measurement measure(int qi, bool collapse = true);
		Measurement measure(const std::vector<size_t>& qi, bool collapse = true);
		Measurement measure_all(bool collapse = true);
		MeasurementBatch measure_batch(const std::vector<size_t>& qi, int shots = 100);

		// Display
		void printState() const;
		void printProb() const;
		void printDiagram();
	};

}