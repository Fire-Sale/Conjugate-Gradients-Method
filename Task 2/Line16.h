#pragma once

#include "Parameters.h"

void KernelLine16(float (&x)[XDIM][YDIM][ZDIM], float (&y)[XDIM][YDIM][ZDIM], float (&z)[XDIM][YDIM][ZDIM], const float scale1, const float scale2);