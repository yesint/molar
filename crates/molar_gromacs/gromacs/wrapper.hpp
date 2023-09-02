#pragma once

#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/idef.h"

// Helper class which incapsulates Gromacs internals
// that can't be translated to rust
class TprHelper {
public:
    TprHelper(const char* fname);
    //~TprHelper();
    t_topology* get_top();
    float* get_atom_xyz(int ind);
    float *get_box();
private:
    t_topology top;
    t_state state;
    t_inputrec ir;
};
