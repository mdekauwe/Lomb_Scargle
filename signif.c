#include "lombscargle.h"

double getchi2(double dof, double alpha)
{
    double  tol = 1.0E-3;
    double  ac, lm, rm, eps, chi2, za, x;
    int     i, itmax = 100;
    
    if (dof > 30.0) {
        /* approximate a value if degrees of freedom are > 30 based on eqn 1. Sachs, 1984 */
        za   = -getz(alpha);        /* NB: Eq. requires change of sign for percentile */
         x    = 2.0 / 9.0 / dof;
         chi2 = pow((dof * (1.0 - x + za * sqrt(x))), 3.0);
      } else{
         
        lm   = 0.0;
         rm   = 1000.0;
         if (alpha > 0.5)
               eps = (1.0 - alpha) * tol;
         else
            eps = alpha * tol;
     
        for (i = 1; i <= itmax; i++) {
               chi2 = 0.5 * (lm + rm);
            ac   = 1.0 - gammp(0.5 * dof, 0.5 * chi2);
            
               if (fabs(ac - alpha) <= eps) 
                return (chi2);    
       
               if (ac > alpha){
                  lm = chi2;
               }else {
                  rm = chi2;
            }    
       }
     }
    return (chi2);
}

void getdof(control *c, timeseries *ts)
{
    /* degrees of freedom for given window based on Harris, 1978 */
    
    double rn = (double)ts->n50, neff, c50, c2, denom;
    
    if (c->WINDOW_TYPE == RECTANGULAR) {
        c50    = 0.5;
        c2     = 2.0 * c50 * c50;
        denom  = 1.0 + c2 - c2 / rn;
        neff   = rn / denom;
        ts->dof = 2.0 * neff;
    } else if (c->WINDOW_TYPE == HANNING) {
        c50    = 0.167;
        c2     = 2.0 * c50 * c50;
        denom  = 1.0 + c2 - c2 / rn;
        neff   = rn / denom;
        ts->dof = 2.0 * neff;
    } else {
        check_errors("Unknown window", __LINE__);
    }
    
    return;
}
