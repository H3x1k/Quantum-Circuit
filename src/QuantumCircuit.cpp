// QuantumCircuit implementation

#include "QuantumCircuit.hpp"
#include <cmath>
#include <iostream>
#include <bitset>
#include <random>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <map>

#define C std::complex<double>

#define E        2.7182818f

#define PI       3.1415926f
#define PI_4     0.7853981f

#define SQRT2    1.4142135f
#define INVSQRT2 0.707106781187f


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



void QuantumCircuit::H(Index qi) {
	for (int i = 0; i < qi.i.size(); i++) {
		int index = qi.i[i];
		Matrix<C> H(2, 2);
		H(0, 0) = INVSQRT2; H(0, 1) = INVSQRT2;
		H(1, 0) = INVSQRT2; H(1, 1) = -INVSQRT2;

		Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - index - 1));
		Matrix I2 = Matrix<C>::identity(size_t(1) << index);

		Matrix fullGate = I1.tensorProduct(H).tensorProduct(I2);

		stateVector = fullGate * stateVector;

		operations.push_back({ OperationType::H, {index} });
	}
	//operations.push_back({ OperationType::H, qi });
}

void QuantumCircuit::X(Index qi) {
	for (int i = 0; i < qi.i.size(); i++) {
		int index = qi.i[i];
		Matrix<C> X(2, 2, C(0.0, 0.0));
		X(0, 1) = C(1.0, 0.0);
		X(1, 0) = C(1.0, 0.0);

		Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - index - 1));
		Matrix I2 = Matrix<C>::identity(size_t(1) << index);

		Matrix fullGate = I1.tensorProduct(X).tensorProduct(I2);

		stateVector = fullGate * stateVector;

		operations.push_back({ OperationType::X, {index} });
	}
	//operations.push_back({ OperationType::X, {qi} });
}

void QuantumCircuit::Y(Index qi) {
	for (int i = 0; i < qi.i.size(); i++) {
		int index = qi.i[i];
		Matrix<C> Y(2, 2, C(0.0, 0.0));
		Y(0, 1) = C(0.0, 1.0);
		Y(1, 0) = C(0.0, 1.0);

		Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - index - 1));
		Matrix I2 = Matrix<C>::identity(size_t(1) << index);

		Matrix fullGate = I1.tensorProduct(Y).tensorProduct(I2);

		stateVector = fullGate * stateVector;

		operations.push_back({ OperationType::Y, {index} });
	}
}

void QuantumCircuit::Z(Index qi) {
	for (int i = 0; i < qi.i.size(); i++) {
		int index = qi.i[i];
		Matrix<C> Z(2, 2, C(0.0, 0.0));
		Z(0, 0) = C(1.0, 0.0);
		Z(1, 1) = C(-1.0, 0.0);

		Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - index - 1));
		Matrix I2 = Matrix<C>::identity(size_t(1) << index);

		Matrix fullGate = I1.tensorProduct(Z).tensorProduct(I2);

		stateVector = fullGate * stateVector;

		operations.push_back({ OperationType::Z, {index} });
	}
}


void QuantumCircuit::S(Index qi) {
	for (int i = 0; i < qi.i.size(); i++) {
		int index = qi.i[i];
		Matrix<C> S(2, 2, C(0.0, 0.0));
		S(0, 0) = C(1.0, 0.0);
		S(1, 1) = C(0.0, 1.0);

		Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - index - 1));
		Matrix I2 = Matrix<C>::identity(size_t(1) << index);

		Matrix fullGate = I1.tensorProduct(S).tensorProduct(I2);

		stateVector = fullGate * stateVector;

		operations.push_back({ OperationType::S, {index} }); \
	}
}

void QuantumCircuit::Sdag(Index qi) {
	for (int i = 0; i < qi.i.size(); i++) {
		int index = qi.i[i];
		Matrix<C> Sdag(2, 2, C(0.0, 0.0));
		Sdag(0, 0) = C(1.0, 0.0);
		Sdag(1, 1) = C(0.0, -1.0);

		Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - index - 1));
		Matrix I2 = Matrix<C>::identity(size_t(1) << index);

		Matrix fullGate = I1.tensorProduct(Sdag).tensorProduct(I2);

		stateVector = fullGate * stateVector;

		operations.push_back({ OperationType::Sdag, {index} });
	}
}

void QuantumCircuit::T(Index qi) {
	for (int i = 0; i < qi.i.size(); i++) {
		int index = qi.i[i];
		Matrix<C> T(2, 2, C(0.0, 0.0));
		T(0, 0) = C(1.0, 0.0);
		T(1, 1) = C(INVSQRT2, INVSQRT2);

		Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - index - 1));
		Matrix I2 = Matrix<C>::identity(size_t(1) << index);

		Matrix fullGate = I1.tensorProduct(T).tensorProduct(I2);

		stateVector = fullGate * stateVector;

		operations.push_back({ OperationType::T, {index} });
	}
}

void QuantumCircuit::Tdag(Index qi) {
	for (int i = 0; i < qi.i.size(); i++) {
		int index = qi.i[i];
		Matrix<C> Tdag(2, 2, C(0.0, 0.0));
		Tdag(0, 0) = C(1.0, 0.0);
		Tdag(1, 1) = C(INVSQRT2, -INVSQRT2);

		Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - index - 1));
		Matrix I2 = Matrix<C>::identity(size_t(1) << index);

		Matrix fullGate = I1.tensorProduct(Tdag).tensorProduct(I2);

		stateVector = fullGate * stateVector;

		operations.push_back({ OperationType::Tdag, {index} });
	}
}


void QuantumCircuit::RX(Index qi, Angle angle) {
	for (int i = 0; i < qi.i.size(); i++) {
		int index = qi.i[i];
		Matrix<C> RX(2, 2, C(0.0, 0.0));
		double theta = angle.get() / 2;
		double costheta = cos(theta);
		double sintheta = sin(theta);
		RX(0, 0) = C(costheta, 0.0);
		RX(0, 1) = C(0.0, -sintheta);
		RX(1, 0) = C(0.0, -sintheta);
		RX(1, 1) = C(costheta, 0.0);

		Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - index - 1));
		Matrix I2 = Matrix<C>::identity(size_t(1) << index);

		Matrix fullGate = I1.tensorProduct(RX).tensorProduct(I2);

		stateVector = fullGate * stateVector;

		operations.push_back({ OperationType::RX, {index}, angle.get() });
	}
}

void QuantumCircuit::RY(Index qi, Angle angle) {
	for (int i = 0; i < qi.i.size(); i++) {
		int index = qi.i[i];
		Matrix<C> RY(2, 2, C(0.0, 0.0));
		double theta = angle.get() / 2;
		double costheta = cos(theta);
		double sintheta = sin(theta);
		RY(0, 0) = C(costheta, 0.0);
		RY(0, 1) = C(-sintheta, 0.0);
		RY(1, 0) = C(-sintheta, 0.0);
		RY(1, 1) = C(costheta, 0.0);

		Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - index - 1));
		Matrix I2 = Matrix<C>::identity(size_t(1) << index);

		Matrix fullGate = I1.tensorProduct(RY).tensorProduct(I2);

		stateVector = fullGate * stateVector;

		operations.push_back({ OperationType::RY, {index}, angle.get() });
	}
}

void QuantumCircuit::RZ(Index qi, Angle angle) {
	for (int i = 0; i < qi.i.size(); i++) {
		int index = qi.i[i];
		Matrix<C> RZ(2, 2, C(0.0, 0.0));
		C expVal = std::exp(C(0.0, angle.get() * 0.5));
		RZ(0, 0) = 1.0 / expVal;
		RZ(1, 1) = expVal;

		Matrix I1 = Matrix<C>::identity(size_t(1) << (numQubits - index - 1));
		Matrix I2 = Matrix<C>::identity(size_t(1) << index);

		Matrix fullGate = I1.tensorProduct(RZ).tensorProduct(I2);

		stateVector = fullGate * stateVector;

		operations.push_back({ OperationType::RZ, {index}, angle.get() });
	}
}


void QuantumCircuit::CNOT(int ci, int ti) {
	for (size_t i = 0; i < stateVector.rows; i++) {
		if (((i >> ci) & 1) == 1) {
			size_t j = i ^ (1ull << ti);
			if (i < j) {
				std::swap(stateVector(i, 0), stateVector(j, 0));
			}
		}
	}

	operations.push_back({ OperationType::CNOT, {ci, ti} });
}

void QuantumCircuit::CZ(Index ci) {
	for (size_t i = 0; i < stateVector.rows; i++) {
		bool allOne = true;
		for (size_t c : ci.i) {
			if (!((i >> c) & 1)) {
				allOne = false;
				break;
			}
		}
		if (allOne) {
			stateVector(i, 0) *= -1.0;
		}
	}

	std::vector<int> i(ci.i.begin(), ci.i.end()); // temp
	operations.push_back({ OperationType::CZ, i });
}

void QuantumCircuit::Rm(int ci, int ti, int m) {
	double theta = 2.0 * PI / (1ull << m);
	C phase = std::exp(C(0, theta));
	for (size_t i = 0; i < stateVector.rows; i++)
		if (((i >> ci) & 1) && ((i >> ti) & 1))
			stateVector(i, 0) *= phase;
}

void QuantumCircuit::Rmdag(int ci, int ti, int m) {
	double theta = -2.0 * PI / (1ull << m);
	C phase = std::exp(C(0, theta));
	for (size_t i = 0; i < stateVector.rows; i++)
		if (((i >> ci) & 1) && ((i >> ti) & 1))
			stateVector(i, 0) *= phase;
}

void QuantumCircuit::SWAP(int q1, int q2) {
	if (q1 == q2) return;

	for (size_t i = 0; i < stateVector.rows; i++) {
		if (((i >> q1) & 1) != ((i >> q2) & 1)) {
			size_t j = i ^ ((1ull << q1) | (1ull << q2));
			if (i < j) { // to avoid double swapping
				std::swap(stateVector(i, 0), stateVector(j, 0));
			}
		}
	}
}


void QuantumCircuit::QFT(Index qi) {
	size_t n = qi.i.size();
	for (size_t i = 0; i < n; i++) {
		int index = qi.i[i];
		H(index);
		for (size_t j = i + 1; j < n; j++) {
			int jindex = qi.i[j];
			Rm(jindex, index, j - i);
		}
	}
	for (size_t i = 0; i < n / 2; i++) {
		SWAP(qi.i[i], qi.i[n - i - 1]);
	}
}

void QuantumCircuit::IQFT(Index qi) {
	size_t n = qi.i.size();
	for (size_t i = 0; i < n / 2; i++) {
		SWAP(qi.i[i], qi.i[n - i - 1]);
	}
	for (int i = n - 1; i >= 0; i--) {
		int index = qi.i[i];
		for (size_t j = n - 1; j > i; j--) {
			int jindex = qi.i[j];
			Rmdag(jindex, index, j - i + 1);
		}
		H(index);
	}
}


Matrix<C> constructPermutationMatrix(int numQubits, const std::vector<size_t>& targets) {
	int dim = 1 << numQubits; // 2^n
	int k = targets.size();

	// Record which qubits are targets
	std::vector<int> targetPos(numQubits, -1);
	for (int i = 0; i < k; i++) {
		int q = targets[i];
		if (q < 0 || q >= numQubits)
			throw std::invalid_argument("Target qubit index out of range");
		targetPos[q] = i; // position among targets
	}

	Matrix<C> P(dim, dim, static_cast<C>(0));

	for (int basis = 0; basis < dim; basis++) {
		int targetBits = 0;
		int nonTargetBits = 0;
		int nonTargetShift = 0;

		// extract bits from basis
		for (int q = 0; q < numQubits; q++) {
			int bit = (basis >> q) & 1;
			if (targetPos[q] != -1) {
				// pack target qubits into lower k bits
				targetBits |= (bit << targetPos[q]);
			}
			else {
				// pack non-target qubits above them
				nonTargetBits |= (bit << nonTargetShift);
				nonTargetShift++;
			}
		}

		// new index = non-targets in higher bits | targets in lower bits
		int newIndex = (nonTargetBits << k) | targetBits;

		// permutation matrix: column = old index, row = new index
		P(newIndex, basis) = static_cast<C>(1);
	}

	return P;
}


void QuantumCircuit::apply(Gate g) { // needs optimization
	int n = g.target.i.size();
	Matrix<C> P = constructPermutationMatrix(numQubits, g.target.i);
	Matrix<C> I = Matrix<C>::identity(numQubits - n);
	Matrix<C> G = I.tensorProduct(g.matrix);
	Matrix<C> U = P * G * P.inverse();
	stateVector = U * stateVector;
	// add operation to operations
}

Matrix<C> permuteMatrix(const Matrix<C>& U, const std::vector<size_t>& perm) {
	int dim = U.rows;
	Matrix<C> result(dim, dim);

	assert(perm.size() == U.rows);  // for rows
	for (size_t k = 0; k < perm.size(); ++k)
		assert(perm[k] < U.rows);


	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			result(perm[i], perm[j]) = U(i, j);

	return result;
}

void QuantumCircuit::apply_controlled(Gate g, Index ci) { // only single qubit
	int n = g.target.i.size();
	int d = 1 << n; // dimension of target gate

	// Gather control + targets
	std::vector<size_t> allIndices = g.target.i;
	allIndices.insert(allIndices.begin(), ci.i[0]);

	// Permutation to bring control + targets to front
	Matrix<C> P = constructPermutationMatrix(numQubits, allIndices);

	// Build controlled matrix (2d x 2d)
	Matrix<C> CU(2 * d, 2 * d, C(0));

	// Top-left block = identity
	for (int i = 0; i < d; i++)
		CU(i, i) = C(1);

	// Bottom-right block = g.matrix
	for (int i = 0; i < d; i++)
		for (int j = 0; j < d; j++)
			CU(d + i, d + j) = g.matrix(i, j);

	// Expand with identities on the remaining qubits
	Matrix<C> I_rest = Matrix<C>::identity(1 << (numQubits - (n + 1)));
	Matrix<C> U = I_rest.tensorProduct(CU);

	// Permute back
	U = permuteMatrix(U, allIndices);
	//U = P * U * P.inverse();
	stateVector = U * stateVector;

	// add to operations history
}

void QuantumCircuit::apply_controlled_test(const Gate& g, const Index& ci) {
	int n = static_cast<int>(g.target.i.size());  // number of target qubits
	if (n == 0) return;

	size_t numStates = 1ULL << numQubits;       // total state vector size
	size_t d = 1ULL << n;                       // size of target gate
	size_t controlPos = ci.i[0];

	// Sanity check: control qubit not in target
	for (size_t t : g.target.i)
		if (t == controlPos)
			throw std::invalid_argument("Control qubit cannot be a target qubit");

	// Precompute bit masks for target qubits
	std::vector<size_t> targetMasks(n);
	for (int k = 0; k < n; ++k)
		targetMasks[k] = 1ULL << g.target.i[k];

	// Temporary storage for amplitudes of a target block
	std::vector<C> in(d), out(d);

	// Loop over all basis states
	for (size_t idx = 0; idx < numStates; ++idx) {
		// Skip states where control qubit is 0
		if (((idx >> controlPos) & 1ULL) == 0) continue;

		// Compute target sub-index
		size_t targetSub = 0;
		for (int k = 0; k < n; ++k)
			if (idx & targetMasks[k]) targetSub |= 1ULL << k;

		// Only process each block once (targetSub == 0)
		if (targetSub != 0) continue;

		// Gather amplitudes for this block
		for (size_t j = 0; j < d; ++j) {
			size_t idxj = idx;
			for (int k = 0; k < n; ++k)
				if (j & (1ULL << k)) idxj |= targetMasks[k];
			in[j] = stateVector(idxj, 0);
		}

		// Apply the gate: out = g.matrix * in
		for (size_t r = 0; r < d; ++r) {
			C acc = 0;
			for (size_t c = 0; c < d; ++c)
				acc += g.matrix(r, c) * in[c];
			out[r] = acc;
		}

		// Write back results
		for (size_t j = 0; j < d; ++j) {
			size_t idxj = idx;
			for (int k = 0; k < n; ++k)
				if (j & (1ULL << k)) idxj |= targetMasks[k];
			stateVector(idxj, 0) = out[j];
		}
	}
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



std::map<std::string, double> QuantumCircuit::probabilityDistribution() {
	std::map<std::string, double> dist;
	for (size_t i = 0; i < stateVector.rows; i++) {
		std::string bitString;
		double prob = std::norm(stateVector(i, 0));
		for (size_t j = 0; j < numQubits; j++) {
			char bit = ((i >> (numQubits - 1 - j)) & 1) ? '1' : '0';
			bitString += bit;
		}
		dist.insert({ bitString, prob });
	}
	return dist;
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

static void drawControlledGate(std::vector<std::string>& canvas, int ci, int ti, std::string c, int numQubits) {
	int l = c.length() + 4;
	int ll = (l % 2 == 0) ? std::max(0, (l / 2) - 1) : l / 2;
	int lr = l / 2;
	if (ci < ti) {
		for (int i = 0; i < numQubits; i++) {
			int y = i * 3 + 1;
			if (i == ci) {
				canvas[y - 1] += std::string(l, ' ');
				canvas[y]     += std::string(ll, char(196)) + std::string(1, char(254)) + std::string(lr, char(196));
				canvas[y + 1] += std::string(ll, ' ') + std::string(1, char(179)) + std::string(lr, ' ');
			} else if (i == ti) {
				canvas[y - 1] += std::string(1, char(218)) + std::string(std::max(0, ll - 1), char(196)) + std::string(1, char(193)) + std::string(std::max(0, lr - 1), char(196)) + std::string(1, char(191));
				canvas[y] += std::string(1, char(180)) + " " + c + " " + std::string(1, char(195));
				canvas[y + 1] += std::string(1, char(192)) + std::string(l - 2, char(196)) + std::string(1, char(217));
			} else if (i > ci && i < ti) {
				canvas[y - 1] += std::string(ll, ' ') + std::string(1, char(179)) + std::string(lr, ' ');
				canvas[y] += std::string(ll, char(196)) + std::string(1, char(197)) + std::string(lr, char(196));
				canvas[y + 1] += std::string(ll, ' ') + std::string(1, char(179)) + std::string(lr, ' ');
			} else {
				canvas[y - 1] += std::string(l, ' ');
				canvas[y] += std::string(l, char(196));
				canvas[y + 1] += std::string(l, ' ');
			}
		}
	} else {
		for (int i = 0; i < numQubits; i++) {
			int y = i * 3 + 1;
			if (i == ci) {
				canvas[y - 1] += std::string(ll, ' ') + std::string(1, char(179)) + std::string(lr, ' ');
				canvas[y] += std::string(ll, char(196)) + std::string(1, char(254)) + std::string(lr, char(196));
				canvas[y + 1] += std::string(l, ' ');
			}
			else if (i == ti) {
				canvas[y - 1] += std::string(1, char(218)) + std::string(l - 2, char(196)) + std::string(1, char(191));
				canvas[y] += std::string(1, char(180)) + " " + c + " " + std::string(1, char(195));
				canvas[y + 1] += std::string(1, char(192)) + std::string(std::max(0, ll - 1), char(196)) + std::string(1, char(194)) + std::string(std::max(0, lr - 1), char(196)) + std::string(1, char(217));
			}
			else if (i < ci && i > ti) {
				canvas[y - 1] += std::string(ll, ' ') + std::string(1, char(179)) + std::string(lr, ' ');
				canvas[y] += std::string(ll, char(196)) + std::string(1, char(197)) + std::string(lr, char(196));
				canvas[y + 1] += std::string(ll, ' ') + std::string(1, char(179)) + std::string(lr, ' ');
			}
			else {
				canvas[y - 1] += std::string(l, ' ');
				canvas[y] += std::string(l, char(196));
				canvas[y + 1] += std::string(l, ' ');
			}
		}
	}
}


static void drawBoxGate(std::vector<std::string>& canvas, int ti, std::string c, int numQubits) {
	size_t l = c.length();
	for (int i = 0; i < numQubits; i++) {
		int y = i * 3 + 1;
		if (i == ti) {
			canvas[y - 1] += std::string(1, char(218)) + std::string(l + 2, char(196)) + std::string(1, char(191));
			canvas[y]     += std::string(1, char(180)) + " " + c + " " + std::string(1, char(195));
			canvas[y + 1] += std::string(1, char(192)) + std::string(l + 2, char(196)) + std::string(1, char(217));
		} else {
			canvas[y - 1] += std::string(l + 4, ' ');
			canvas[y]     += std::string(l + 4, char(196));
			canvas[y + 1] += std::string(l + 4, ' ');
		}
	}
}

static std::string doubleToString(double value, int precision) {
	std::ostringstream out;
	out << std::fixed << std::setprecision(precision) << value;
	std::string s = out.str();
	s.erase(s.find_last_not_of('0') + 1, std::string::npos);
	if (!s.empty() && s.back() == '.') s.pop_back();
	return s;
}

void QuantumCircuit::printDiagram() const {
	const int numOps = operations.size();
	const int rows = numQubits * 3;

	std::vector<std::string> canvas(rows, "");

	for (int i = 0; i < rows; i++) {
		if ((i - 1) % 3 == 0) {
			int qn = (i - 1) / 3;
			canvas[i] += "q" + std::to_string(qn) + " " + char(196);
		} else {
			canvas[i] += "    ";
		}
	}

	for (int i = 0; i < numOps; i++) {
		switch (operations[i].type) {
		case OperationType::H: {
			drawBoxGate(canvas, operations[i].qubits[0], "H", numQubits);
			break;
		}
		case OperationType::X: {
			drawBoxGate(canvas, operations[i].qubits[0], "X", numQubits);
			break;
		}
		case OperationType::Y: {
			drawBoxGate(canvas, operations[i].qubits[0], "Y", numQubits);
			break;
		}
		case OperationType::Z: {
			drawBoxGate(canvas, operations[i].qubits[0], "Z", numQubits);
			break;
		}
		case OperationType::RX: {
			std::string rparam = "RX(" + doubleToString(operations[i].parameter, 3) + ")";
			drawBoxGate(canvas, operations[i].qubits[0], rparam, numQubits);
			break;
		}
		case OperationType::RY: {
			std::string rparam = "RY(" + doubleToString(operations[i].parameter, 3) + ")";
			drawBoxGate(canvas, operations[i].qubits[0], rparam, numQubits);
			break;
		}
		case OperationType::RZ: {
			std::string rparam = "RZ(" + doubleToString(operations[i].parameter, 3) + ")";
			drawBoxGate(canvas, operations[i].qubits[0], rparam, numQubits);
			break;
		}
		case OperationType::CNOT: {
			int ci = operations[i].qubits[0];
			int ti = operations[i].qubits[1];
			drawControlledGate(canvas, ci, ti, "X", numQubits);
			break;
		}
		case OperationType::CZ: {
			int ci = operations[i].qubits[0];
			int ti = operations[i].qubits[1];
			drawControlledGate(canvas, ci, ti, "Z", numQubits);
			break;
		}
		default: {
			break;
		}
		}
	}

	for (int i = 0; i < numQubits; i++) {
		int y = i * 3 + 1;
		canvas[y - 1] += " "; canvas[y] += char(196); canvas[y + 1] += " ";
	}

	for (int i = 0; i < canvas.size(); i++) {
		std::cout << canvas[i] << std::endl;
	}
}