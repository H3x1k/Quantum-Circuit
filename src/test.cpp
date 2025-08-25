#include <vector>
#include <complex>
#include <cmath>
#include <iostream>

#include "QuantumCircuit.hpp"
#include "Matrix.hpp"
//#include "Gate.hpp"

int main() {
	int nq = 1;
	qcf::QuantumCircuit qc(nq);

	qc.H(0);

	qc.printProb();

	qcf::MeasurementBatch mb = qc.measure_batch({ 0 }, 1024);
	mb.print();



	char he = (char)196;
	char ve = (char)179;
	char tlc = (char)218;
	char trc = (char)191;
	char blc = (char)192;
	char brc = (char)217;
	char vel = (char)180;
	char ver = (char)195;

	std::cout << "    " << tlc << he << he << he << trc << " " << std::endl;
	std::cout << "q0 " << he << vel << " " << "H" << " " << ver << he << std::endl;
	std::cout << "    " << blc << he << he << he << brc << " " << std::endl;



	while (1);

	return 0;
}