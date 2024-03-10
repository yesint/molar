#include <molar_sasa/Eigen/Dense>
#include "molar_sasa/powersasa/power_sasa.h"
#include "molar_sasa/powersasa/wrapper.h"
#include <memory>
#include "molar_sasa/src/lib.rs.h"

// Minimal iterator-like class for coordinates
class SasaCoordIter {
public:
    SasaCoordIter(RustCoordProvider& provider): provider(provider) {}
    SasaCoordIter& operator++() { provider.advance(); return *this; }
    Eigen::Vector3f const& operator*() const {
        // Take a raw pointer from rust provider and reinterpret as *Vector3f
        return *reinterpret_cast<Eigen::Vector3f const*>(provider.get());
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


// A minimal Container-like class wrapping a Rust coord provider
class SasaCoordContainer {
public:
    SasaCoordContainer(RustCoordProvider& provider, size_t sz): provider(provider), sz(sz) {}    
    SasaCoordIter begin() const {
        provider.reset();
        return SasaCoordIter(provider);
    }
    int size() const {return sz;}
private:
    RustCoordProvider& provider;
    size_t sz;
};

// A minimal Container-like class wrapping a Rust coord provider
class SasaVDWContainer {
public:
    SasaVDWContainer(RustVdwProvider& provider, size_t sz): provider(provider), sz(sz) {}    
    SasaVDWIter begin() const {
        provider.reset();
        return SasaVDWIter(provider);
    }
    int size() const {return sz;}
private:
    RustVdwProvider& provider;
    size_t sz;
};

// Constructor function
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