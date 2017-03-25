/* 
Unit:
                         1 length = 1 Mpc/h
                         1 mass = 1 M_sun /h
*/

/* Unit transformation constants */
#define SUNTOKG               1.9891e30
#define MPCTOM                3.08567758e22
#define G0                    6.67384e-11
#define PI                    3.14159265358979323

/* ST parameters */
#define Ast                   0.3222
#define pst                   0.3
#define qst                   0.707

/* 1D Spline related parameters */
#define Mmin                  1e7               /* Mmin < 1e5 is enough for Kmax < 100 h/Mpc */
#define Mmax                  1e19              /* Mmax > 1e18 */
#define NI                    500               /* the interpolation point number */
