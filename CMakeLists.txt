cmake_minimum_required(VERSION 3.15)

project(QuantumFramework VERSION 0.1 LANGUAGES CXX)

# Use C++17 standard (or change to your preferred version)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories (headers)
include_directories(include)

# Gather all source files
file(GLOB SOURCES "src/*.cpp")

# Create the library (static or shared)
add_library(quantum_framework STATIC ${SOURCES})

# (Optional) Enable testing if you add tests later
enable_testing()
