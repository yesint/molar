#pragma once
#include "stddef.h"
#include <memory>

struct RustCoordProvider;
struct RustVdwProvider;

// A type exposed to Rust
class PowerSasaWrapper
{
public:
    PowerSasaWrapper(RustCoordProvider &crd_prov, RustVdwProvider &vdw_prov, size_t sz);

    void get_next_coord() const;
private:
    RustCoordProvider& crd;
    RustVdwProvider& vdw;
};

std::unique_ptr<PowerSasaWrapper> new_power_sasa(
    RustCoordProvider &crd_prov,
    RustVdwProvider &vdw_prov,
    size_t sz);
