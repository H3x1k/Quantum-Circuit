// QuantumCircuit header

#include <vector>
#include <complex>
#include <map>

#include "Matrix.hpp"
#include "Angle.hpp"

namespace qcf {

	enum class OperationType {
		H, X, Y, Z, S, Sdag, T, Tdag, RX, RY, RZ, CNOT, CZ, Measure
	};

	struct Measurement {
		std::vector<uint8_t> bits{};
		double probability = 0.0;

		void print() const;
	};

	struct MeasurementBatch {
		std::vector<size_t> qubits{};
		std::map<std::string, size_t> counts{};
		int shotCount;

		void print(bool fraction = false) const;
	};

	struct Operation {
		OperationType type;
		std::vector<int> qubits{};
		double parameter = 0.0;
	};

	class QuantumCircuit {
	public:

		size_t numQubits;
		Matrix<std::complex<double>> stateVector;
		std::vector<Operation> operations;

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
		void RX(int qi, Angle angle);
		void RY(int qi, Angle angle);
		void RZ(int qi, Angle angle);
		// Controlled Gates
		void CNOT(int ci, int ti);
		void CZ(int ci, int ti);

		// Measurement
		Measurement measure(int qi, bool collapse = true, bool saveOp = true);
		Measurement measure(const std::vector<size_t>& qi, bool collapse = true, bool saveOp = true);
		Measurement measure_all(bool collapse = true, bool saveOp = true);
		MeasurementBatch measure_batch(const std::vector<size_t>& qi, int shots = 100);

		// Display
		void printState() const;
		void printProb() const;
		void printDiagram() const;
	};

}