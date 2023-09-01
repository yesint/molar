#!/bin/bash

# To be run in the root of Gromacs source directory
# wrapper.hpp have to be copied or linked there

bindgen --no-layout-tests --allowlist-type "t_topology t_state t_inputrec" --allowlist-function read_tpr \
wrapper.hpp -o gromacs_bindings.rs -- -x c++ -std=c++17 \
-I ./src \
-I ./src/gromacs/utility/include \
-I ./src/gromacs/math/include \
-I ./src/gromacs/topology/include \
-I ./api/legacy/include \
-I ./build/api/legacy/include \
-I ./src/external \
