#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"

#include "wrapper.hpp"

/* Internal structure — Gromacs types are fully hidden from the C ABI. */
struct TprHandle {
    gmx_mtop_t               mtop;
    t_topology               top;
    t_state                  state;
    t_inputrec               ir;
    std::vector<TprBond>     bonds;
    std::vector<TprMolecule> molecules;
};

static thread_local std::string s_last_error;

extern "C" {

TprHandle* tpr_open(const char* path)
{
    try {
        TprHandle* h = new TprHandle();

        read_tpx_state(path, &h->ir, &h->state, &h->mtop);
        h->top = gmx_mtop_t_to_t_topology(&h->mtop, true);

        /* --- Parse bonds from idef ----------------------------------------
         * idef.il is indexed by function type (0..F_NRE).
         * Each il[ftype].iatoms contains tuples:
         *   [param_idx, atom1, atom2]        for 2-body interactions
         *   [param_idx, atom1, atom2, atom3] for F_SETTLE
         */
        for (int ftype = 0; ftype < F_NRE; ++ftype) {
            const t_ilist& il = h->top.idef.il[ftype];
            if (il.nr == 0) continue;

            switch (ftype) {
                case F_BONDS:
                case F_G96BONDS:
                case F_HARMONIC:
                case F_FENEBONDS:
                case F_CUBICBONDS:
                case F_CONSTR:
                case F_CONSTRNC:
                    for (int i = 0; i < il.nr; i += 3) {
                        h->bonds.push_back({
                            (uint32_t)il.iatoms[i + 1],
                            (uint32_t)il.iatoms[i + 2]
                        });
                    }
                    break;

                case F_SETTLE:
                    /* Each SETTLE entry covers O-H1-H2; represents 2 bonds. */
                    for (int i = 0; i < il.nr; i += 4) {
                        h->bonds.push_back({(uint32_t)il.iatoms[i+1], (uint32_t)il.iatoms[i+2]});
                        h->bonds.push_back({(uint32_t)il.iatoms[i+1], (uint32_t)il.iatoms[i+3]});
                    }
                    break;

                default:
                    break;
            }
        }

        /* --- Parse molecules from mols block --------------------------------
         * t_block.nr  = number of molecules
         * t_block.index[i] = first atom of molecule i
         * t_block.index[i+1] = first atom of molecule i+1
         *
         * We replicate the Rust logic: read mols.nr entries and chunk by 2,
         * treating them as (start, end_exclusive) pairs.
         */
        {
            int  mol_nr    = h->top.mols.nr;
            int* mol_index = h->top.mols.index;
            for (int i = 0; i + 1 < mol_nr; i += 2) {
                h->molecules.push_back({
                    (uint32_t)mol_index[i],
                    (uint32_t)(mol_index[i + 1] - 1)
                });
            }
        }

        return h;
    }
    catch (const std::exception& e) {
        s_last_error = e.what();
        return nullptr;
    }
    catch (...) {
        s_last_error = "unknown error in tpr_open";
        return nullptr;
    }
}

void tpr_close(TprHandle* h)
{
    if (!h) return;
    done_top_mtop(&h->top, &h->mtop);
    delete h;
}

const char* tpr_last_error(void)
{
    return s_last_error.c_str();
}

size_t tpr_natoms(TprHandle* h)      { return (size_t)h->top.atoms.nr; }
size_t tpr_nbonds(TprHandle* h)      { return h->bonds.size(); }
size_t tpr_nmolecules(TprHandle* h)  { return h->molecules.size(); }

void tpr_fill_atoms(TprHandle* h, TprAtom* out)
{
    int   natoms    = h->top.atoms.nr;
    auto* atoms     = h->top.atoms.atom;
    auto* atomnames = h->top.atoms.atomname;
    auto* resinfo   = h->top.atoms.resinfo;
    auto* atomtypes = h->top.atoms.atomtype;
    auto* pdbinfo   = h->top.atoms.pdbinfo;

    for (int i = 0; i < natoms; ++i) {
        TprAtom& a = out[i];
        memset(&a, 0, sizeof(TprAtom));

        if (atomnames[i] && *atomnames[i])
            strncpy(a.name, *atomnames[i], 7);

        int resi = atoms[i].resind;

        if (resinfo[resi].name && *resinfo[resi].name)
            strncpy(a.resname, *resinfo[resi].name, 7);

        char chain = resinfo[resi].chainid;
        a.chain = (chain == '\0') ? ' ' : chain;

        if (atomtypes[i] && *atomtypes[i])
            strncpy(a.type_name, *atomtypes[i], 7);

        a.resind        = (uint32_t)resi;
        a.type_id       = (uint32_t)atoms[i].type;
        a.atomic_number = (uint32_t)atoms[i].atomnumber;
        a.charge        = atoms[i].q;
        a.mass          = atoms[i].m;

        if (pdbinfo) {
            a.occupancy = pdbinfo[i].occup;
            a.bfactor   = pdbinfo[i].bfac;
        }
    }
}

void tpr_fill_bonds(TprHandle* h, TprBond* out)
{
    memcpy(out, h->bonds.data(), h->bonds.size() * sizeof(TprBond));
}

void tpr_fill_molecules(TprHandle* h, TprMolecule* out)
{
    memcpy(out, h->molecules.data(), h->molecules.size() * sizeof(TprMolecule));
}

void tpr_fill_coords(TprHandle* h, float* out)
{
    int natoms = h->top.atoms.nr;
    for (int i = 0; i < natoms; ++i) {
        const float* xyz = reinterpret_cast<const float*>(&h->state.x[i]);
        out[i * 3 + 0]   = xyz[0];
        out[i * 3 + 1]   = xyz[1];
        out[i * 3 + 2]   = xyz[2];
    }
}

void tpr_fill_box(TprHandle* h, float* out9)
{
    /* state.box is float[3][3] in row-major order — copy the raw bytes.
     * Rust side calls Matrix3::from_column_slice on these 9 floats. */
    memcpy(out9, h->state.box, 9 * sizeof(float));
}

} /* extern "C" */
