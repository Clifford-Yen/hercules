#include <math.h>
#include "basin.h"

int getBasinMaterialProperties(cvmpayload_t *g_props, double z_m) {
    // TODO: This is a hardcoded part. It would be great if these formulas could 
    // be read from a file. But this is pretty hard to be achieved in C.
    g_props->rho = 2140.0 + 0.125*z_m;
    g_props->Vp = 1000.0 + 1.2*z_m;        
    g_props->Vs = 320.0 + 19.0*sqrt(z_m);
    g_props->Qs = 0.4*g_props->Vs; // Not sure why Qs is 0.4*Vs instead of 0.1*Vs like in the 1D velocity model.
    g_props->Qp = 0.5*g_props->Qs;
    return 0;
}