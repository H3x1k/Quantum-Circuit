// QuantumCircuit header

#include <vector>
#include <complex>

#include "Matrix.hpp"

namespace qcf {

	struct Measurement {
		std::vector<uint8_t> bits{};
		double probability = 0.0;
	};

	struct MeasurementBatch {

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
		Measurement measure(int qi);
		Measurement measure(const std::vector<size_t>& qi);
		Measurement measure_all();
		MeasurementBatch measure_batch(int shots);


		void printState() const;
		void printProb() const;
		//void printDiagram();
	};

}