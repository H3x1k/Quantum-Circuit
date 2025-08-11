// Gate header

#include <string>
#include <vector>
#include <complex>

namespace qcf {
	class Gate {
	public:
		std::string name;
		std::vector<std::complex<double>> matrix;
		Gate();
	};
}