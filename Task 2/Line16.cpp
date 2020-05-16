#include "Line16.h"
#include <algorithm>


void KernelLine16(float (&x)[XDIM][YDIM][ZDIM], float (&y)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM],
    const float scale1, const float scale2)
{
    // Should we use OpenMP parallel for here?
// #pragma omp parallel for    
    for (int i = 1; i < XDIM-1; i++)
    for (int j = 1; j < YDIM-1; j++)
    for (int k = 1; k < ZDIM-1; k++)
        y[i][j][k] = x[i][j][k] * scale1 + y[i][j][k];

    for (int i = 1; i < XDIM-1; i++)
    for (int j = 1; j < YDIM-1; j++)
    for (int k = 1; k < ZDIM-1; k++)
        x[i][j][k] = x[i][j][k] * scale2+ z[i][j][k];

}