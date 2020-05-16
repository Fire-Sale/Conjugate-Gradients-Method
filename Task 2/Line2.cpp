#include "Line2.h"
#include <algorithm>


float KernelLine2(float (&x)[XDIM][YDIM][ZDIM], float (&y)[XDIM][YDIM][ZDIM],
    const float (&z)[XDIM][YDIM][ZDIM], float (&r)[XDIM][YDIM][ZDIM])
{
    #pragma omp parallel for
    for (int i = 1; i < XDIM-1; i++)
    for (int j = 1; j < YDIM-1; j++)
    for (int k = 1; k < ZDIM-1; k++)
        y[i][j][k] =
            -6 * x[i][j][k]
            + x[i+1][j][k]
            + x[i-1][j][k]
            + x[i][j+1][k]
            + x[i][j-1][k]
            + x[i][j][k+1]
            + x[i][j][k-1];

    float result = 0.;

#pragma omp parallel for reduction(max:result)
    for (int i = 1; i < XDIM-1; i++)
    for (int j = 1; j < YDIM-1; j++)
    for (int k = 1; k < ZDIM-1; k++)
    {
        r[i][j][k] = -y[i][j][k] + z[i][j][k];
        result = std::max(result, std::abs(r[i][j][k]));
    }
        
    return result;

}