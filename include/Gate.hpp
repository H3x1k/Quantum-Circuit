// Gate header

#include <string>
#include <vector>
#include <complex>

#include "Matrix.hpp"
#include "Index.hpp"

namespace qcf {
	class Gate {
	public:
		Matrix<std::complex<double>> matrix;
		Index target;
		Gate(const Matrix<std::complex<double>>& mat, Index target);
	};
}