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
        .function("S", &QuantumCircuit::S)
        .function("Sdag", &QuantumCircuit::Sdag)
        .function("T", &QuantumCircuit::T)
        .function("Tdag", &QuantumCircuit::Tdag)
        .function("RX", &QuantumCircuit::RX)
        .function("RY", &QuantumCircuit::RY)
        .function("RZ", &QuantumCircuit::RZ)
        .function("CNOT", &QuantumCircuit::CNOT)
        .function("CZ", &QuantumCircuit::T)
        .function("Rm", &QuantumCircuit::Rm)
        .function("Rmdag", &QuantumCircuit::Rmdag)
        .function("SWAP", &QuantumCircuit::SWAP)
        .function("QFT", &QuantumCircuit::QFT)
        .function("IQFT", &QuantumCircuit::IQFT)
        //.function("apply", &QuantumCircuit::apply) // untested
        //.function("apply_controlled", &QuantumCircuit::apply_controlled) // unfinished
        //.function("measure", &QuantumCircuit::measure) // unfinished
        .function("printState", &QuantumCircuit::printState)
        .function("printProb", &QuantumCircuit::printProb)
        .function("probabilityDistribution", &probabilityDistributionToJS)
        .function("toQASM", &QuantumCircuit::toQASM);

    class_<Gate>("Gate")
        .constructor<const Matrix<std::complex<double>>&, Index>();

    class_<Index>("Index")
        .constructor<size_t>()
        .constructor<std::vector<size_t>>()
        .class_function("range", &Index::range);

    // Register std::vector<size_t> for gate operations
    register_vector<size_t>("VectorSizeT");
}