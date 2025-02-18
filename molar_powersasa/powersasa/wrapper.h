#pragma once
#include "stddef.h"

void run_powersasa(
    float *area_per_atom,             // Returns areas per atom in this array. Has to point to properly sized array!
    float *volume_per_atom,           // Returns volume per atom in this array. Has to point to properly sized array!
    float *(*cb_crd)(void *, size_t), // callback for coordinates
    float (*cb_vdw)(void *, size_t),  // callback VdW radii. !!! They must be VdW+probe_r !!!
    void *context_crd,                // Context pointer for coordinates
    void *context_vdw,                // Context pointer for VdW
    size_t num                        // Number of atoms
);
