// Gate implementation

#include "Gate.hpp"

using namespace qcf;

Gate::Gate(const std::vector<std::vector<std::complex<double>>>& mat, int target) : matrix(mat), target(target) {

}