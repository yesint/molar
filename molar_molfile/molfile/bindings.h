#pragma once

#include "molfile_plugin.h"
//struct molfile_plugin_t;

extern molfile_plugin_t* pdb_get_plugin_ptr();
extern molfile_plugin_t* xyz_get_plugin_ptr();
extern molfile_plugin_t* dcd_get_plugin_ptr();
