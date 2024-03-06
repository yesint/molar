#pragma once
#include "stddef.h"

void powersasa(    
    float* area_per_atom, // Returns areas per atom in this array. Has to point to properly sized array!    
    float* volume_per_atom, // Returns volume per atom in this array. Has to point to properly sized array!
    float* (*cb)(size_t,void*), // callback for coordinates
    float (*cb_vdw)(size_t,void*), // callback VdW radii. !!! They must be VdW+probe_r !!!
    size_t num, // Numebr of atoms
    void* context
    ); 
