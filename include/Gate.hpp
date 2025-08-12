// Gate header

#include <string>
#include <vector>
#include <complex>

namespace qcf {
	class Gate {
	public:
		std::vector<std::vector<std::complex<double>>> matrix;
		size_t target;
		Gate(const std::vector<std::vector<std::complex<double>>>& mat, int target);
	};
}