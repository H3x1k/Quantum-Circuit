#include "QuantumCircuit.hpp"


// Articles
// https://www.nature.com/articles/s41534-022-00583-7
// https://www.nature.com/articles/s41567-022-01539-6?fromPaywallRec=false
// https://en.wikipedia.org/wiki/Quantum_logic_gate

// TODO:
// add display for multi-qubit gates and multi control qubits
// add time component for operations so that multiple gates can be displayed on one line

// add multi-qubit measurement to finish shor's algorithm
// find the highest occurring number then use continued fractions to approximate the computed fraction from c/2^w
// then use r to compute factors of N using gcd(a^r/2 +- 1, N)

// optimize state vector operations
// optimize matrix operations
// add more measurement or simulations
// add noise models and realistic hardware simulation
// different simulation modes (statevector, stabilizer (Clifford circuits), density matrix, tensor network (scale to more qubits))
// integrate with standard formats (OpenQASM, Cirq, Qiskit)
// multi-threading, gqu acceleration (CUDA)
// high-level api

// Prioritize QASM support, noise models, and optimization


int main() {

	qcf::QuantumCircuit qc(2);

	qc.H(0);

	qc.printProb();
	qc.printState();

	while (1);
	return 0;
}