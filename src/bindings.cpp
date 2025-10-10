#include <emscripten/bind.h>
#include <emscripten/val.h>
#include "QuantumCircuit.hpp"
#include "Index.hpp"

using namespace emscripten;
using namespace qcf;

EMSCRIPTEN_BINDINGS(quantum_module) {
    class_<QuantumCircuit>("QuantumCircuit")
        .constructor<int>()

        // Single qubit gates - use proper select_overload syntax
        .function("H", &QuantumCircuit::H)
        .function("X", &QuantumCircuit::X)
        .function("Y", &QuantumCircuit::Y)
        .function("Z", &QuantumCircuit::Z)

        // QFT
        .function("IQFT", &QuantumCircuit::IQFT)

        // Measurement
        .function("measure_all", &QuantumCircuit::measure_all)
        .function("probabilityDistribution", &QuantumCircuit::probabilityDistribution)

        // Custom gate application
        .function("apply_controlled_test", &QuantumCircuit::apply_controlled_test)

        // Print functions
        .function("printState", &QuantumCircuit::printState)
        .function("printProb", &QuantumCircuit::printProb);

    emscripten::class_<Gate>("Gate")
        .constructor<>();

    class_<Index>("Index")
        .class_function("range", &Index::range);

    // Register std::vector<size_t> for gate operations
    register_vector<size_t>("VectorSizeT");
}