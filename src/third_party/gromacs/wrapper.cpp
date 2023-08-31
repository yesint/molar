#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/idef.h"

t_topology read_tpr() {
    t_inputrec ir;    
    gmx_mtop_t mtop;
    //t_topology top;
    t_state state;

    read_tpx_state(fname.c_str(), &ir, &state, &mtop);

    return( gmx_mtop_t_to_t_topology(&mtop,false) );
}
