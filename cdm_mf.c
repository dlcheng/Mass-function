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


/*********************************     Spline functions related to v(M) relation *************************/
void init_spline_mu()
{
   int i;	
   double sig_m;

/* initialize LogM */	
   for(i=0; i<NI; i++)
      {	
      SplogM[i] = log10(Mmin) + i * (log10(Mmax) - log10(Mmin)) / (double)(NI - 1);
      sig_m = sqrt(sigma_m_sq(pow(10, SplogM[i])));
      SpMu[i] = delta_z / sig_m;
      }
      
   spacc_mu = gsl_interp_accel_alloc();    
   sp_mu = gsl_spline_alloc(gsl_interp_cspline, NI);
   
   gsl_spline_init(sp_mu, SplogM, SpMu, NI);
		
}               /* end init_spline_mu */	

/* Interface */
double mu(double M)
{
  double temp_logm = log10(M);
  double result = gsl_spline_eval(sp_mu, temp_logm, spacc_mu);
  
  return result;	
}              /* end mu */
 
double dlogsigmasq_dlogm(double M)
{
  double temp_logm = log10(M);
  double result = -2.0 / mu(M) / log(10) * gsl_spline_eval_deriv(sp_mu, temp_logm, spacc_mu);
  
  return result;
}            /* end dlogsigmasq_dlogm */
 
/*********************************     Spline functions related to CDM mass function *************************/
void init_spline_ncdm()
{
  int i;
  double temp_mu;
  double temp_M;
  double temp_dsm_dm;
  
  for(i=0; i<NI; i++)
    {
	  temp_M = pow(10, SplogM[i]);
	  temp_mu = mu(temp_M);
	  temp_dsm_dm = dlogsigmasq_dlogm(temp_M);
#ifdef ST_MF        
    SpNcdm[i] = fabs(-0.5 * rho_m / temp_M / temp_M * fst_v(temp_mu) * temp_dsm_dm);
#endif
#ifdef PLANCK_500
    SpNcdm[i] = fabs(-0.5 * rho_m / temp_M / temp_M * fplanck_v(temp_mu) * temp_dsm_dm);
#endif    
    }
    
  spacc_ncdm = gsl_interp_accel_alloc();    
  sp_ncdm = gsl_spline_alloc(gsl_interp_cspline, NI);
   
   gsl_spline_init(sp_ncdm, SplogM, SpNcdm, NI);    
}             /* end init_spline_ncdm */

double fst_v(double v)
{
  double result = Ast * sqrt(2.0/PI) * sqrt(qst) * v;

  result *= 1.0 + pow(sqrt(qst) * v, -2.0 * pst);
  result *= exp(-qst / 2.0 * v * v);

  return result;
}      /* end fst_v */     


double fplanck_v(double v)                    /* the overdensity is 500 to the critical density at z=0 */
{
 double Aplanck = 0.261067162295;
 double aplanck = 2.32644914136;
 double bplanck = 1.45015783928;
 double cplanck = 1.99566917921;

 double result;

 result = Aplanck * (1+ pow((delta_c / v / bplanck), -aplanck)) * exp(- cplanck * v * v / delta_c / delta_c);

 return result;
} /* end fplanck_v */     

/* Interface */
double ncdm(double M)
{
  double temp_logm = log10(M);	
  double result = gsl_spline_eval(sp_ncdm, temp_logm, spacc_ncdm);
  
  return result;	
}	         /* end ncdm */
