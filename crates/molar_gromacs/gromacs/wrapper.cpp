#include "wrapper.hpp"

TprHelper::TprHelper(const char *fname)
{
    gmx_mtop_t mtop;
    read_tpx_state(fname, &ir, &state, &mtop);
    top = gmx_mtop_t_to_t_topology(&mtop,false);
}

//TprHelper::~TprHelper()
//{
//}

t_topology *TprHelper::get_top()
{
    return &top;
}

float *TprHelper::get_atom_xyz(int ind)
{
    return (float*)&state.x[ind];
}

float *TprHelper::get_box()
{
    return (float*)&state.box;
}
