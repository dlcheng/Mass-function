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

int main()
{
/* Step 1
 * 1) Cosmological parameters.
 * 2) 2D table details
 * 3) Global initialization 
 */ 	
 
  set_cosmological_parameters(0.3, 0.7, (0.024/0.7/0.7), 0.7, 0.96, 0.8);   /* parameter of our CDM simulation */
  init_same_time(0);                                                        /* the parameter is the redshift */
  out_put(1e9, 1e15, 20);
  free_1d_spline();

 // set_cosmological_parameters(0.3175, 0.6825, (0.02207/0.674/0.674), 0.674, 0.9624, 0.8344);   /* best fit of Planck CMB */
 // init_same_time(0);                                                        /* the parameter is the redshift */
 // out_put(1e13, 1e17, 20);
 // free_1d_spline();

 // set_cosmological_parameters(0.2814, 0.7186, 0.0464, 0.697, 1.0, 0.9);   /* parameter of our CDM simulation */
 // init_same_time(0); 
 // out_put(1e11, 1e16, 40);
 // free_1d_spline();

  return 0;
}   /* end main */	

void out_put(double m1, double m2, int num)
{
   double log_mass_dis = log10(m2/m1)/(double) num;
   int i;
   double mass;
   double mass_function;

  printf("Mass     Mass-function  Mass*MF\n");

   for(i=0; i<num; i++)
   {
    mass = log10(m1) + (i+0.5) * log_mass_dis;
    mass = pow(10, mass);
    mass_function = ncdm(mass);
    printf("%.5e  %.5e  %.5e\n",mass, mass_function, mass*mass_function*log(10));
   }
}  /* end out_put */ 
