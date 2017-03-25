#ifndef ALLVAR_H
 #include "allvars.h"
#endif

void init_spline_mu();
double mu(double M);
double dlogsigmasq_dlogm(double M);
void init_spline_ncdm();
double fst_v(double v);
double fplanck_v(double v);
double ncdm(double M);

double win_1(double x);

double g_factor(double z);
double unnorm_groth(double z);
void init_growth();
double growth_factor(double z);

int init_global();
int init_same_time(double a);

void set_cosmological_parameters(double x1, double x2, double x3, double x4, double x5, double x6);

double prim_tk_bbks(double k);

double prim_tk(double k);

double radius_from_mass(double M);
double vir_radius(double M, double z);
double mass_from_r(double R);
void free_1d_spline();
void mark(char *s);

double sigma_m_int_kernel(double k, void *param);
double unnorm_sigma_m_sq(double M);
void init_sigma_m();
double sigma_m_sq(double M);

void out_put(double m1, double m2, int num);
