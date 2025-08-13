// QuantumCircuit implementation

#include "QuantumCircuit.hpp"
#include <cmath>

using namespace qcf;

QuantumCircuit::QuantumCircuit(int numQubits) : numQubits(numQubits), stateVector(size_t(1) << numQubits, {0.0, 0.0}) {
	stateVector[0] = { 1.0, 0.0 };
}

void QuantumCircuit::addGate(const Gate& gate) {
	operations.push_back(gate);
}