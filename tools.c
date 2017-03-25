#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>

#include "allvars.h"
#include "proto.h"

double radius_from_mass(double M)
{
  double result = M * 3.0 / 4.0 / PI / rho_m;	
		
  return pow(result , 1.0/3.0);	
}      /* end radius_from_mass */	

double vir_radius(double M, double z)
{ 
  double rho_crit_z = rho_crit * (Omega_v + Omega_m * pow(1.0 + z, 3)); /* the critical density at z, used to define halos */
  double result = M * 3.0 / 4.0 / PI / (rho_crit_z * 200);  
  
   result = pow(result , 1.0/3.0);
   result = result * (1.0 + z);            /* the co-moving virial radius_from_mass */	
   
  return result;		
}      /* end vir_radius */	
  
double mass_from_r(double R)
{
  return 4.0 * PI / 3.0 * pow(R, 3.0) * rho_m;	
}    /* end mass */

void free_1d_spline()
{
  gsl_spline_free(sp_mu);
  gsl_spline_free(sp_ncdm);
  
  gsl_interp_accel_free(spacc_mu);
  gsl_interp_accel_free(spacc_ncdm);
		
}   /* end free_spline_objects */	

void mark(char *s)
{
 printf("%s\n", s);	
}	
