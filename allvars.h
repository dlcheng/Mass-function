#ifndef ALLVAR_H
#define ALLVAR_H

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

#include "define.h"


/* Cosmological parameters */
extern double Omega_m;
extern double Omega_v;
extern double Omega_b;
extern double H0;
extern double ns;         
extern double sigma_8;
extern double delta_c;

/* Halo model  */
extern double Aps;                       /* normalization factor for the PS */
extern double Agr;                       /* normalization factor for the growth function */

/* Density */
extern double rho_crit;
extern double rho_m;

/* others */
extern double delta_z;

/* 1D Spline related */
extern double SplogM[NI];
extern double SpMu[NI];
extern double SpNcdm[NI];

extern gsl_interp_accel * spacc_mu;
extern gsl_interp_accel * spacc_ncdm;

extern gsl_spline *sp_mu;
extern gsl_spline *sp_ncdm;


#endif
