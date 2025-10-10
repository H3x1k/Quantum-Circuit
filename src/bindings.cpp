#include <emscripten/bind.h>
#include <emscripten/val.h>
#include "QuantumCircuit.hpp"
#include "Index.hpp"

using namespace emscripten;
using namespace qcf;

val probabilityDistributionToJS(QuantumCircuit& qc) {
    auto dist = qc.probabilityDistribution();
    val obj = val::object();

    for (const auto& pair : dist) {
        obj.set(pair.first, pair.second);
    }

    return obj;
}

EMSCRIPTEN_BINDINGS(quantum_module) {

    register_map<std::string, double>("MapStringDouble");
    register_vector<size_t>("VectorSizeT");

    class_<QuantumCircuit>("QuantumCircuit")
        .constructor<int>()
        .function("H", &QuantumCircuit::H)
        .function("X", &QuantumCircuit::X)
        .function("Y", &QuantumCircuit::Y)
        .function("Z", &QuantumCircuit::Z)
        .function("printState", &QuantumCircuit::printState)
        .function("printProb", &QuantumCircuit::printProb)
        .function("probabilityDistribution", &probabilityDistributionToJS);

    class_<Gate>("Gate")
        .constructor<const Matrix<std::complex<double>>&, Index>();

    class_<Index>("Index")
        .constructor<size_t>()
        .constructor<std::vector<size_t>>()
        .class_function("range", &Index::range);

    // Register std::vector<size_t> for gate operations
    register_vector<size_t>("VectorSizeT");
}