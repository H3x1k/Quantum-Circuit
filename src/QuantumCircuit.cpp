// QuantumCircuit implementation

#include "QuantumCircuit.hpp"
#include <cmath>
#include <iostream>
#include <bitset>
#include <random>

#define C std::complex<double>

#define E        2.7182818f

#define PI       3.1415926f
#define PI_4     0.7853981f

#define SQRT2    1.4142135f
#define INVSQRT2 0.7071067f


static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<double> dist(0.0, 1.0);


using namespace qcf;

void Measurement::print() const {
	std::cout << "Measured: ";
	for (int i = bits.size() - 1; i >= 0; i--)
		std::cout << int(bits[i]);
	std::cout << std::endl << "Probability: " << probability << std::endl;
}

void MeasurementBatch::print(bool fraction) const {
	std::cout << "Shot Count : " << shotCount << std::endl;
	std::cout << "Outcomes (" << (fraction ? "fraction" : "count") << ")" << std::endl;
	for (const auto& [key, value] : counts) {
		std::cout << key << " : " << (fraction ? (float)value / shotCount : value) << std::endl;
	}
	std::cout << std::endl;
}


QuantumCircuit::QuantumCircuit(int numQubits)
	: numQubits(numQubits),
	stateVector((size_t(1) << numQubits), 1, C(0.0, 0.0)) {
	stateVector(0, 0) = C(1.0, 0.0);
	operations = {};
}



void QuantumCircuit::H(int qi) {

	Matrix<C> H(2, 2);
	H(0, 0) = INVSQRT2; H(0, 1) = INVSQRT2;
	H(1, 0) = INVSQRT2; H(1, 1) = -INVSQRT2;
	
	Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));
	Matrix I2 = Matrix<C>::identity(size_t(1) << qi);

	Matrix fullGate = I1.tensorProduct(H).tensorProduct(I2);

	stateVector = fullGate * stateVector;

	operations.push_back({ OperationType::H, {qi} });
}

void QuantumCircuit::X(int qi) {

	Matrix<C> X(2, 2, C(0.0, 0.0));
	X(0, 1) = C(1.0, 0.0);
	X(1, 0) = C(1.0, 0.0);

	Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));
	Matrix I2 = Matrix<C>::identity(size_t(1) << qi);

	Matrix fullGate = I1.tensorProduct(X).tensorProduct(I2);

	stateVector = fullGate * stateVector;

	operations.push_back({ OperationType::X, {qi} });
}

void QuantumCircuit::Y(int qi) {
	Matrix<C> Y(2, 2, C(0.0, 0.0));
	Y(0, 1) = C(0.0, 1.0);
	Y(1, 0) = C(0.0, 1.0);

	Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));
	Matrix I2 = Matrix<C>::identity(size_t(1) << qi);

	Matrix fullGate = I1.tensorProduct(Y).tensorProduct(I2);

	stateVector = fullGate * stateVector;

	operations.push_back({ OperationType::Y, {qi} });
}

void QuantumCircuit::Z(int qi) {
	Matrix<C> Z(2, 2, C(0.0, 0.0));
	Z(0, 0) = C(1.0, 0.0);
	Z(1, 1) = C(-1.0, 0.0);

	Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));
	Matrix I2 = Matrix<C>::identity(size_t(1) << qi);

	Matrix fullGate = I1.tensorProduct(Z).tensorProduct(I2);

	stateVector = fullGate * stateVector;

	operations.push_back({ OperationType::Z, {qi} });
}


void QuantumCircuit::S(int qi) {
	Matrix<C> S(2, 2, C(0.0, 0.0));
	S(0, 0) = C(1.0, 0.0);
	S(1, 1) = C(0.0, 1.0);

	Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));
	Matrix I2 = Matrix<C>::identity(size_t(1) << qi);

	Matrix fullGate = I1.tensorProduct(S).tensorProduct(I2);

	stateVector = fullGate * stateVector;

	operations.push_back({ OperationType::S, {qi} });
}

void QuantumCircuit::Sdag(int qi) {
	Matrix<C> Sdag(2, 2, C(0.0, 0.0));
	Sdag(0, 0) = C(1.0,  0.0);
	Sdag(1, 1) = C(0.0, -1.0);

	Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));
	Matrix I2 = Matrix<C>::identity(size_t(1) << qi);

	Matrix fullGate = I1.tensorProduct(Sdag).tensorProduct(I2);

	stateVector = fullGate * stateVector;

	operations.push_back({ OperationType::Sdag, {qi} });
}

void QuantumCircuit::T(int qi) {
	Matrix<C> T(2, 2, C(0.0, 0.0));
	T(0, 0) = C(1.0, 0.0);
	T(1, 1) = C(INVSQRT2, INVSQRT2);

	Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));
	Matrix I2 = Matrix<C>::identity(size_t(1) << qi);

	Matrix fullGate = I1.tensorProduct(T).tensorProduct(I2);

	stateVector = fullGate * stateVector;

	operations.push_back({ OperationType::T, {qi} });
}

void QuantumCircuit::Tdag(int qi) {
	Matrix<C> Tdag(2, 2, C(0.0, 0.0));
	Tdag(0, 0) = C(1.0, 0.0);
	Tdag(1, 1) = C(INVSQRT2, -INVSQRT2);

	Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));
	Matrix I2 = Matrix<C>::identity(size_t(1) << qi);

	Matrix fullGate = I1.tensorProduct(Tdag).tensorProduct(I2);

	stateVector = fullGate * stateVector;

	operations.push_back({ OperationType::Tdag, {qi} });
}


void QuantumCircuit::RX(int qi, double angle) {
	Matrix<C> RX(2, 2, C(0.0, 0.0));
	double theta = angle / 2;
	double costheta = cos(theta);
	double sintheta = sin(theta);
	RX(0, 0) = C(costheta,  0.0);
	RX(0, 1) = C(0.0, -sintheta);
	RX(1, 0) = C(0.0, -sintheta);
	RX(1, 1) = C(costheta,  0.0);

	Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));
	Matrix I2 = Matrix<C>::identity(size_t(1) << qi);

	Matrix fullGate = I1.tensorProduct(RX).tensorProduct(I2);

	stateVector = fullGate * stateVector;

	operations.push_back({ OperationType::RX, {qi}, angle });
}

void QuantumCircuit::RY(int qi, double angle) {
	Matrix<C> RY(2, 2, C(0.0, 0.0));
	double theta = angle / 2;
	double costheta = cos(theta);
	double sintheta = sin(theta);
	RY(0, 0) = C(costheta,  0.0);
	RY(0, 1) = C(-sintheta, 0.0);
	RY(1, 0) = C(-sintheta, 0.0);
	RY(1, 1) = C(costheta,  0.0);
	
	Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));
	Matrix I2 = Matrix<C>::identity(size_t(1) << qi);

	Matrix fullGate = I1.tensorProduct(RY).tensorProduct(I2);

	stateVector = fullGate * stateVector;

	operations.push_back({ OperationType::RY, {qi}, angle });
}

void QuantumCircuit::RZ(int qi, double angle) {
	Matrix<C> RZ(2, 2, C(0.0, 0.0));
	C expVal = std::exp(C(0.0, angle * 0.5));
	RZ(0, 0) = 1.0 / expVal;
	RZ(1, 1) = expVal;

	Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - qi - 1));
	Matrix I2 = Matrix<C>::identity(size_t(1) << qi);

	Matrix fullGate = I1.tensorProduct(RZ).tensorProduct(I2);

	stateVector = fullGate * stateVector;

	operations.push_back({ OperationType::RZ, {qi}, angle });
}


void QuantumCircuit::CNOT(int ci, int ti) {
	for (size_t i = 0; i < stateVector.rows; i++) {
		if (((i >> ci) & 1) && (((i >> ti) & 1) == 0)) {
			size_t j = i ^ (1ull << ti);
			C temp = stateVector(i, 0);
			stateVector(i, 0) = stateVector(j, 0);
			stateVector(j, 0) = temp;
		}
	}

	operations.push_back({ OperationType::CNOT, {ci, ti} });
}



Measurement QuantumCircuit::measure(int qi, bool collapse, bool saveOp) {
	double prob0 = 0.0;

	for (size_t i = 0; i < stateVector.rows; i++) {
		if (((i >> qi) & 1) == 0) {
			prob0 += std::norm(stateVector(i, 0));
		}
	}

	double r = dist(gen);
	Measurement m{};
	if (r < prob0) {
		m.bits = { 0 };
		m.probability = prob0;
	} else {
		m.bits = { 1 };
		m.probability = 1.0 - prob0;
	}

	if (collapse) {
		double norm_factor = 1.0 / std::sqrt(m.probability);
		for (size_t i = 0; i < stateVector.rows; i++) {
			if (((i >> qi) & 1) == m.bits[0]) {
				stateVector(i, 0) *= norm_factor;
			}
			else {
				stateVector(i, 0) = C(0.0, 0.0);
			}
		}
	}

	if (saveOp)
		operations.push_back({ OperationType::Measure, {qi} });

	return m;
}

Measurement QuantumCircuit::measure(const std::vector<size_t>& qi, bool collapse, bool saveOp) {
	const size_t k = qi.size();
	std::vector<double> probs(numQubits, 0.0);

	for (size_t i = 0; i < stateVector.rows; i++) {
		size_t outcome = 0;
		for (size_t j = 0; j < k; j++) {
			size_t bit = (i >> qi[j]) & 1;
			outcome |= (bit << j);
		}
		probs[outcome] += std::norm(stateVector(i, 0));
	}

	double r = dist(gen);
	double cumulativeProb = 0.0;
	size_t index = 0;
	for (size_t i = 0; i < probs.size(); i++) {
		cumulativeProb += probs[i];
		if (r < cumulativeProb) { index = i; break; }
	}

	if (collapse) {
		double norm_factor = 1.0 / std::sqrt(probs[index]);
		for (size_t i = 0; i < stateVector.rows; i++) {
			size_t outcome = 0;
			for (size_t j = 0; j < k; j++) {
				size_t bit = (i >> qi[j]) & 1;
				outcome |= (bit << j);
			}
			if (outcome == index) {
				stateVector(i, 0) *= norm_factor;
			}
			else {
				stateVector(i, 0) = C(0.0, 0.0);
			}
		}
	}

	Measurement m{};
	m.bits.resize(k);
	for (size_t i = 0; i < k; i++) {
		m.bits[i] = (index >> i) & 1;
	}
	m.probability = probs[index];

	return m;
}

Measurement QuantumCircuit::measure_all(bool collapse, bool saveOp) {
	double r = dist(gen);
	double cumulativeProb = 0.0;
	size_t index = 0;
	double indexProb = 0.0;
	for (size_t i = 0; i < stateVector.rows; i++) {
		indexProb = std::norm(stateVector(i, 0));
		cumulativeProb += indexProb;
		if (r < cumulativeProb) { index = i; break; }
	}

	if (collapse) {
		stateVector = Matrix<C>(stateVector.rows, 1, C(0.0, 0.0));
		stateVector(index, 0) = C(1.0, 0.0);
	}

	Measurement m{};
	m.probability = indexProb;
	m.bits.resize(numQubits);
	for (size_t i = 0; i < numQubits; i++) {
		m.bits[i] = (index >> i) & 1;
	}
	return m;
}

MeasurementBatch QuantumCircuit::measure_batch(const std::vector<size_t>& qi, int shots) {
	MeasurementBatch batch{};
	if (qi.size() == 1) {
		size_t index = qi[0];
		batch.shotCount = shots;
		batch.qubits = { index };
		for (int i = 0; i < shots; i++) {
			Measurement m = this->measure(index, false, false);
			batch.counts[(m.bits[0] == 1 ? "1" : "0")]++;
		}
	}

	return batch;
}



void QuantumCircuit::printState() const {
	stateVector.print();
}

void QuantumCircuit::printProb() const {
	for (size_t i = 0; i < stateVector.rows; i++) {
		if (stateVector(i, 0) != C(0.0, 0.0)) {
			double prob = std::norm(stateVector(i, 0));
			std::cout << "|";
			for (size_t j = 0; j < numQubits; j++)
				std::cout << ((i >> (numQubits - 1 - j)) & 1);
			std::cout << "> : " << prob << std::endl;
		}
	}
}

void QuantumCircuit::printDiagram() const {
	const int numOps = operations.size();
	const int rows = numQubits * 6 - 1;
	const int cols = numOps * 7 + 5 * (numOps - 1) + 5;
	std::vector<std::vector<char>> canvas(rows, std::vector<char>(cols, ' '));
	
	for (int i = 0; i < numQubits; i++) {
		int y = 6 * i + 3;
		canvas[y][0] = 'q';
		canvas[y][1] = '0' + i;
		canvas[y][2] = ' ';
		for (int x = 3; x < cols; x++)
			canvas[y][x] = (char)196;
	}

	


	for (int y = 0; y < rows; y++) {
		for (int x = 0; x < cols; x++)
			std::cout << canvas[y][x];
		std::cout << std::endl;
	}
}