#pragma once

#include "xdrfile.h"
#include <stdbool.h>

int xdr_xtc_get_natoms(XDRFILE* handle, int* natoms);
int xtc_get_current_frame_number(XDRFILE* handle, int natoms, bool* bOK);
int xtc_get_next_frame_number(XDRFILE* handle, int natoms);
float xtc_get_current_frame_time(XDRFILE* handle, int natoms, bool* bOK);
float xtc_get_next_frame_time(XDRFILE* handle, int natoms, bool* bOK);
float xdr_xtc_estimate_dt(XDRFILE* handle, int natoms, bool* bOK);
int xdr_xtc_get_last_frame_number(XDRFILE* handle, int natoms, bool* bOK);
float xdr_xtc_get_last_frame_time(XDRFILE* handle, int natoms, bool* bOK);
int xdr_xtc_seek_frame(int frame, XDRFILE* handle, int natoms);
int xdr_xtc_seek_time(float time, XDRFILE* handle, int natoms, bool bSeekForwardOnly);
int check_trr_content(XDRFILE* handle, int* natoms, int* xsz, int* vsz, int* fsz);
