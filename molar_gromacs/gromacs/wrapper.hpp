#pragma once

#include <stddef.h>
#include <stdint.h>

/* Plain C-compatible data structures for the molar TPR plugin.
   No Gromacs types appear here — this header is the stable ABI boundary. */

typedef struct {
    char     name[8];
    char     resname[8];
    char     type_name[8];
    char     chain;
    /* 3 bytes implicit padding before first uint32_t */
    uint32_t resind;
    uint32_t type_id;
    uint32_t atomic_number;
    float    charge;
    float    mass;
    float    occupancy;   /* 0.0 if TPR has no pdbinfo */
    float    bfactor;     /* 0.0 if TPR has no pdbinfo */
} TprAtom;

typedef struct {
    uint32_t atom1;
    uint32_t atom2;
} TprBond;

typedef struct {
    uint32_t start;  /* first atom index, inclusive */
    uint32_t end;    /* last atom index, inclusive */
} TprMolecule;

/* Opaque handles — internals defined only in wrapper.cpp */
typedef struct TprHandle TprHandle;
typedef struct CptHandle CptHandle;

#ifdef __cplusplus
extern "C" {
#endif

/* Open a TPR file; returns NULL on error (call tpr_last_error for message). */
TprHandle*  tpr_open(const char* path);
void        tpr_close(TprHandle* handle);
/* Error message from the last failed tpr_open() or cpt_open(); valid until next call. */
const char* tpr_last_error(void);

/* Counts — computed once at open, O(1) thereafter. */
size_t tpr_natoms(TprHandle* handle);
size_t tpr_nbonds(TprHandle* handle);
size_t tpr_nmolecules(TprHandle* handle);

/* Fill caller-allocated arrays.
   Rust side: allocate Vec<TprX> of the correct count, pass .as_mut_ptr(). */
void tpr_fill_atoms(TprHandle* handle, TprAtom* out);
void tpr_fill_bonds(TprHandle* handle, TprBond* out);
void tpr_fill_molecules(TprHandle* handle, TprMolecule* out);
/* out must hold natoms*3 floats (row-major XYZ). */
void tpr_fill_coords(TprHandle* handle, float* out);
/* out must hold 9 floats (box matrix flattened row-major). */
void tpr_fill_box(TprHandle* handle, float* out9);

/* Open a CPT (checkpoint) file; returns NULL on error (call tpr_last_error for message).
   Reads and caches coordinates + box only; velocities and forces are ignored. */
CptHandle* cpt_open(const char* path);
void       cpt_close(CptHandle* handle);
size_t     cpt_natoms(CptHandle* handle);
float      cpt_time(CptHandle* handle);
int64_t    cpt_step(CptHandle* handle);
/* out must hold natoms*3 floats (row-major XYZ, nm). */
void       cpt_fill_coords(CptHandle* handle, float* out);
/* out must hold 9 floats (box matrix flattened row-major, nm). */
void       cpt_fill_box(CptHandle* handle, float* out9);

#ifdef __cplusplus
}
#endif
