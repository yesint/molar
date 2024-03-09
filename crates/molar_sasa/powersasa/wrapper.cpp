#include <molar_sasa/Eigen/Dense>
#include "molar_sasa/powersasa/power_sasa.h"
#include "molar_sasa/powersasa/wrapper.h"
#include <memory>
#include "molar_sasa/src/main.rs.h"

// Minimal iterator-like class for coordinates
class SasaCoordIter {
public:
    SasaCoordIter(RustCoordProvider& provider): provider(provider) {}
    SasaCoordIter& operator++() { provider.advance(); return *this; }
    Eigen::Vector3f& operator*() const {
        // Take a raw pointer from rust provider and reinterpret as *Vector3f
        return *reinterpret_cast<Eigen::Vector3f*>(provider.get());
    }    
private:
    RustCoordProvider& provider;
};

// Minimal iterator-like class for VdW
class SasaVDWIter {
public:
    SasaVDWIter(RustVdwProvider& provider): provider(provider) {}
    SasaVDWIter& operator++() { provider.advance(); return *this; }
    float operator*() const {
        // Just return a value from rust provider
        return provider.get();
    }    
private:
    RustVdwProvider& provider;
};


// This is a minimal Container-like class, which takes a Rust callback that provides elements
// PowerSasa uses only begin() and size() methods, so no iterator comparison logic is needed
class SasaCoordContainer {
public:
    SasaCoordContainer(RustCoordProvider& provider, size_t sz): provider(provider), sz(sz) {}    
    SasaCoordIter begin() const { return SasaCoordIter(provider); }
    int size() const {return sz;}
private:
    RustCoordProvider& provider;
    size_t sz;
};

// Container for VdW radii
class SasaVDWContainer {
public:
    SasaVDWContainer(RustVdwProvider& provider, size_t sz): provider(provider), sz(sz) {}    
    SasaVDWIter begin() const { return SasaVDWIter(provider); }
    int size() const {return sz;}
private:
    RustVdwProvider& provider;
    size_t sz;
};

/*
void run_powersasa(    
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
*/
//----------------------------------------

std::unique_ptr<PowerSasaWrapper> new_power_sasa(
    RustCoordProvider& crd_prov, 
    RustVdwProvider& vdw_prov,
    size_t sz
) {
    return std::make_unique<PowerSasaWrapper>(crd_prov,vdw_prov,sz);
}

void PowerSasaWrapper::get_next_coord() const {
    std::cout << "crd_prov:" << crd.get()[0] << "\n";
    std::cout << "vsw_prov" << vdw.get() << "\n";
    vdw.advance();
}

PowerSasaWrapper::PowerSasaWrapper(RustCoordProvider &crd_prov, RustVdwProvider &vdw_prov, size_t sz):
    crd(crd_prov),vdw(vdw_prov)
{
    SasaCoordContainer cont_crd(crd,sz);
    SasaVDWContainer cont_vdw(vdw,sz);

    std::cout << "Computing SASA...\n";
    // Call POWERSASA
    POWERSASA::PowerSasa<float,Eigen::Vector3f> ps(cont_crd, cont_vdw, 1, 0, 1, 0);
    ps.calc_sasa_all();
    std::cout << "First area: " << ps.getSasa()[0] << "\n";
}