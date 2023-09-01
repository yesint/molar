#pragma once

#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/idef.h"

extern void read_tpr(const char* fname, t_topology* top, t_state* state, t_inputrec* ir);
