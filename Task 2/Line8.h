#pragma once

#include "Parameters.h"

float KernelLine8(const float (&x)[XDIM][YDIM][ZDIM], const float (&y)[XDIM][YDIM][ZDIM], float (&z)[XDIM][YDIM][ZDIM], const float scale);

