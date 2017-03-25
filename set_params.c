#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void set_cosmological_parameters(double x1, double x2, double x3, double x4, double x5, double x6)
{
  Omega_m = x1;
  Omega_v = x2;
  Omega_b = x3;
  H0      = x4;
  ns      = x5;             
  sigma_8 = x6;
  delta_c = 1.686;	  

  init_global();  
}    /* end set_cosmological_parameters */