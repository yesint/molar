#include <Eigen/Core>
#include "power_sasa.h"
#include "wrapper.hpp"

// Minimal iterator-like class for coordinates
class SasaCoordIter {
public:
    using value_type = Eigen::Vector3f ;
    using difference_type = int;
    using pointer = value_type*;
    using reference = value_type& ;
    using iterator_category = std::forward_iterator_tag;

    SasaCoordIter(float* (*cb)(void*,size_t), void* context): callback(cb), pos(0), context(context) {}

    SasaCoordIter operator++(int junk) {auto tmp = *this; ++pos; return tmp;}
    SasaCoordIter& operator++() { ++pos; return *this; }

    reference operator*() const {
        // Take a raw pointer from callback and reinterpret as *Vector3f
        return *reinterpret_cast<pointer>(callback(context,pos));
    }    
private:
    float* (*callback)(void*,size_t);
    size_t pos;
    void* context;
};

// Minimal iterator-like class for VdW
class SasaVDWIter {
public:
    using value_type = float;
    using difference_type = int;
    using pointer = value_type*;
    using reference = value_type& ;
    using iterator_category = std::forward_iterator_tag;

    SasaVDWIter(float (*cb)(void*,size_t), void* context): callback(cb), pos(0), context(context) {}

    SasaVDWIter operator++(int junk) {auto tmp = *this; ++pos; return tmp;}
    SasaVDWIter& operator++() { ++pos; return *this; }

    float operator*() const {
        // Just return a value from callback
        return callback(context,pos);
    }    
private:
    float (*callback)(void*,size_t);
    size_t pos;
    void* context;
};


// This is a minimal Container-like class, which takes a Rust callback that provides elements
// PowerSasa uses only begin() and size() methods, so no iterator comparison logic is needed
class SasaCoordContainer {
public:
    SasaCoordContainer(float* (*cb)(void*,size_t), size_t sz, void* context): callback(cb), sz(sz), context(context) {}    
    SasaCoordIter begin() const { return SasaCoordIter(callback,context); }
    int size() const {return sz;}
private:
    float* (*callback)(void*,size_t);
    size_t sz;
    void* context;
};

// Container for VdW radii
class SasaVDWContainer {
public:
    SasaVDWContainer(float (*cb)(void*,size_t), size_t sz, void* context): callback(cb), sz(sz), context(context) {}    
    SasaVDWIter begin() const { return SasaVDWIter(callback,context); }
    int size() const {return sz;}
private:
    float (*callback)(void*,size_t);
    size_t sz;
    void* context;
};


void powersasa(    
    float* area_per_atom, // Returns areas per atom in this array. Has to point to properly sized array!    
    float* volume_per_atom, // Returns volume per atom in this array. Has to point to properly sized array!
    float* (*cb_crd)(void*,size_t), // callback for coordinates
    float (*cb_vdw)(void*,size_t), // callback VdW radii. !!! They must be VdW+probe_r !!!
    void* context_crd, // Context pointer for coordinates
    void* context_vdw, // Context pointer for VdW
    size_t num // Number of atoms
    )
{        
    SasaCoordContainer cont_crd(cb_crd,num,context_crd);
    SasaVDWContainer cont_vdw(cb_vdw,num,context_vdw);

    
    // Call POWERSASA
    POWERSASA::PowerSasa<float,Eigen::Vector3f> ps(cont_crd, cont_vdw, 1, 0, 1, 0);
    ps.calc_sasa_all();

    float v,surf=0.0;

    if(volume_per_atom){        
        for(int i = 0; i < num; ++i){
            volume_per_atom[i] = ps.getVol()[i];
        }
    }

    if(area_per_atom){
        for(int i = 0; i < num; ++i){                    
            area_per_atom[i] = ps.getSasa()[i];
        }
    }
    
}
