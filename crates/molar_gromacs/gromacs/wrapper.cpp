#include "wrapper.hpp"

void read_tpr(const char* fname, t_topology* top, t_state* state, t_inputrec* ir) {
    gmx_mtop_t mtop;
    read_tpx_state(fname, ir, state, &mtop);
    *top = gmx_mtop_t_to_t_topology(&mtop,false);
}
