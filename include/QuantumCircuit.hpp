// QuantumCircuit header

#include <vector>
#include <complex>
#include <map>

#include "Matrix.hpp"
#include "Angle.hpp"
#include "Index.hpp"
#include "Gate.hpp"

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
		void H(Index qi);
		// Pauli Gates
		void X(Index qi);
		void Y(Index qi);
		void Z(Index qi);
		// Phase Gates
		void S(Index qi); void Sdag(Index qi);
		void T(Index qi); void Tdag(Index qi);
		// Rotation Gates
		void RX(Index qi, Angle angle);
		void RY(Index qi, Angle angle);
		void RZ(Index qi, Angle angle);
		// Controlled Gates
		void CNOT(int ci, int ti);
		void CZ(Index ci);
		void Rm(int ci, int ti, int m); void Rmdag(int ci, int ti, int m);
		void SWAP(int q1, int q2); // only takes two qubits
		// Composite Gates
		void QFT(Index qi); void IQFT(Index qi);
		// Custom Gates
		void apply(Gate g);
		void apply_controlled(Gate g, Index ci);
		void apply_controlled_test(const Gate& g, const Index& ci);

		// Measurement
		Measurement measure(int qi, bool collapse = true, bool saveOp = true);
		Measurement measure(const std::vector<size_t>& qi, bool collapse = true, bool saveOp = true);
		Measurement measure_all(bool collapse = true, bool saveOp = true);
		MeasurementBatch measure_batch(const std::vector<size_t>& qi, int shots = 100); // mulitple doesn't work

		// Display
		void printState() const;
		void printProb() const;
		void printDiagram() const;
	};

}