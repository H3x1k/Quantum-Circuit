// Gate implementation

#include "Gate.hpp"

using namespace qcf;

Gate::Gate(const Matrix<std::complex<double>>& mat, Index target) : matrix(mat), target(target) {}