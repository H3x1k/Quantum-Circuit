// QuantumCircuit header

#include <vector>
#include <complex>

#include "Matrix.hpp"

namespace qcf {

	class QuantumCircuit {
	public:

		size_t numQubits;
		Matrix<std::complex<double>> stateVector;

		QuantumCircuit(int numQubits);
		
		void H(int qi);


		//const Matrix<std::complex<double>>& getState() const;
		void printState() const;
		//void printDiagram();
	};

}