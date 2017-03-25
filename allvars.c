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

/* Cosmological parameters */
 double Omega_m;
 double Omega_v;
 double Omega_b;
 double H0;
 double ns;         
 double sigma_8;
 double delta_c;

/* Halo model  */
 double Aps;                       /* normalization factor for the PS */
 double Agr;                       /* normalization factor for the growth function */

/* Density */
 double rho_crit;
 double rho_m;

/* others */
 double delta_z;

/* 1D Spline related */
 double SplogM[NI];
 double SpMu[NI];
 double SpNcdm[NI];

 gsl_interp_accel * spacc_mu;
 gsl_interp_accel * spacc_ncdm;

 gsl_spline *sp_mu;
 gsl_spline *sp_ncdm;

