#!/bin/bash

bindgen --no-layout-tests --allowlist-type t_topology \
bindings.hpp -o bindings.rs -- \
-I ./src \
-I ./src/gromacs/utility/include \
-I ./src/gromacs/math/include \
-I ./src/gromacs/topology/include \
-I ./api/legacy/include \
-I ./build/api/legacy/include \
-I ./src/external \
