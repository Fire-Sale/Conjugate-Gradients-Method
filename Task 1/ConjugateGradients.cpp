#include "Laplacian.h"
#include "Parameters.h"
#include "PointwiseOps.h"
#include "Reductions.h"
#include "Utilities.h"
#include "Timer.h"

#include <iostream>

extern Timer timerTotal;
extern Timer timerLaplacian;
extern Timer timerInnerProduct;
extern Timer timerSaxpy;
extern Timer timerNorm;
extern Timer timerCopy;

void ConjugateGradients(
    float (&x)[XDIM][YDIM][ZDIM],
    const float (&f)[XDIM][YDIM][ZDIM],
    float (&p)[XDIM][YDIM][ZDIM],
    float (&r)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM],
    const bool writeIterations)
{
    timerTotal.Restart();

    // Algorithm : Line 2
    timerLaplacian.Restart();
    ComputeLaplacian(x, z);
    timerLaplacian.Pause();

    timerSaxpy.Restart();
    Saxpy(z, f, r, -1);
    timerSaxpy.Pause();

    timerNorm.Restart();
    float nu = Norm(r);
    timerNorm.Pause();

    // Algorithm : Line 3
    if (nu < nuMax) return;
        
    // Algorithm : Line 4
    timerCopy.Restart();
    Copy(r, p);
    timerCopy.Pause();
    timerInnerProduct.Restart();
    float rho=InnerProduct(p, r);
    timerInnerProduct.Pause();
        
    // Beginning of loop from Line 5
    for(int k=0;;k++)
    {
        std::cout << "Residual norm (nu) after " << k << " iterations = " << nu << std::endl;

        // Algorithm : Line 6
        timerLaplacian.Restart();
        ComputeLaplacian(p, z); 
        timerLaplacian.Pause();

        timerInnerProduct.Restart();
        float sigma=InnerProduct(p, z);
        timerInnerProduct.Pause();

        // Algorithm : Line 7
        float alpha=rho/sigma;

        // Algorithm : Line 8
        timerSaxpy.Restart();
        Saxpy(z, r, r, -alpha);
        timerSaxpy.Pause();

        timerNorm.Restart();
        nu=Norm(r);
        timerNorm.Pause();

        // Algorithm : Lines 9-12
        if (nu < nuMax || k == kMax) {
            timerSaxpy.Restart();
            Saxpy(p, x, x, alpha);
            timerSaxpy.Pause();
            timerTotal.Pause();
            std::cout << "Conjugate Gradients terminated after " << k << " iterations; residual norm (nu) = " << nu << std::endl;
            if (writeIterations) WriteAsImage("x", x, k, 0, 127);
            return;
        }
            
        // Algorithm : Line 13
        timerCopy.Restart();
        Copy(r, z);
        timerCopy.Pause();
        timerInnerProduct.Restart();
        float rho_new = InnerProduct(z, r);
        timerInnerProduct.Pause();

        // Algorithm : Line 14
        float beta = rho_new/rho;

        // Algorithm : Line 15
        rho=rho_new;

        // Algorithm : Line 16
        timerSaxpy.Restart();
        Saxpy(p, x, x, alpha);
        Saxpy(p, r, p, beta);
        timerSaxpy.Pause();

        if (writeIterations) WriteAsImage("x", x, k, 0, 127);
    }
    
}
