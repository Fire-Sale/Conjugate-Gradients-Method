#include "ConjugateGradients.h"
#include "Timer.h"
#include "Utilities.h"

Timer timerTotal;
Timer timerLaplacian;
Timer timerInnerProduct;
Timer timerSaxpy;
Timer timerNorm;
Timer timerCopy;

int main(int argc, char *argv[])
{
    using array_t = float (&) [XDIM][YDIM][ZDIM];

    float *xRaw = new float [XDIM*YDIM*ZDIM];
    float *fRaw = new float [XDIM*YDIM*ZDIM];
    float *pRaw = new float [XDIM*YDIM*ZDIM];
    float *rRaw = new float [XDIM*YDIM*ZDIM];
    float *zRaw = new float [XDIM*YDIM*ZDIM];
    
    array_t x = reinterpret_cast<array_t>(*xRaw);
    array_t f = reinterpret_cast<array_t>(*fRaw);
    array_t p = reinterpret_cast<array_t>(*pRaw);
    array_t r = reinterpret_cast<array_t>(*rRaw);
    array_t z = reinterpret_cast<array_t>(*zRaw);
    

    // Initialization
    {
        Timer timer;
        timer.Start();
        InitializeProblem(x, f);
        timer.Stop("Initialization : ");
    }

    // Call Conjugate Gradients algorithm
    timerTotal.Reset();
    timerLaplacian.Reset();
    timerInnerProduct.Reset();
    timerSaxpy.Reset();
    timerNorm.Reset();
    timerCopy.Reset();
    ConjugateGradients(x, f, p, r, z);
    timerLaplacian.Print("Total Laplacian Time : ");
    timerInnerProduct.Print("Total Inner Product Time : ");
    timerSaxpy.Print("Total Saxpy Time : ");
    timerNorm.Print("Total Norm Time : ");
    timerCopy.Print("Total Copy Time : ");
    timerTotal.Print("Total Time : ");
    
    std::cout << "[Total Time (predicted) : " << timerLaplacian.mElapsedTime.count() + timerInnerProduct.mElapsedTime.count() + timerSaxpy.mElapsedTime.count() + timerNorm.mElapsedTime.count() + timerCopy.mElapsedTime.count() << "ms]" << std::endl;
    return 0;
}
