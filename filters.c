#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

double win_1(double x)
{
  double result = 1;
  
  if(x < 1e-8)
    result = 1.0 - 0.1 * x * x;
  else  
	result = 3.0 / pow(x, 3) * (sin(x) - x*cos(x));
		
   return result;		
}           /* end win_1 */

