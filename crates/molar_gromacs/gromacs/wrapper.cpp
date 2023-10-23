#include "wrapper.hpp"

TprHelper::TprHelper(const char *fname)
{
    read_tpx_state(fname, &ir, &state, &mtop);
    top = gmx_mtop_t_to_t_topology(&mtop,true);
}

TprHelper::~TprHelper()
{
    done_top_mtop(&top,&mtop);
}

t_topology *TprHelper::get_top()
{
    return &top;
}

size_t TprHelper::get_natoms()
{
    return top.atoms.nr;
}

float *TprHelper::get_atom_xyz(size_t ind)
{
    return (float*)&state.x[ind];
}

float *TprHelper::get_box()
{
    return (float*)&state.box;
}

char *TprHelper::get_atomname(size_t ind)
{
    return *(top.atoms.atomname[ind]);
}
