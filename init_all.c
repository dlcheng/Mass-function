#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "allvars.h"
#include "proto.h"

int init_global()
{
 rho_crit = 3e10*MPCTOM/(8.0*PI*G0)/SUNTOKG;
 rho_m    = rho_crit*Omega_m;
	
 init_growth();	
 init_sigma_m();   /* Calculate the sigma_8 normalization */	     	   
 
 return 1; 	  	
}  /* end init_step_1 */	


int init_same_time(double z)
{
 delta_z = delta_c / growth_factor(z);
 
 init_spline_mu();
 init_spline_ncdm();	
 
 return 1;
}      /* end init_same_time */
