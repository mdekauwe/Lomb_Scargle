 #include "lombscargle.h"

void gettau(control *c, timeseries *ts)
{
    /* Estimate a value for tau based on Mudelsee, 2002, Computers and Geosciences, 28, 69-72 */
    
    double *xwk = NULL, *ywk = NULL, rhosum = 0., xfit = 0.0, intercept = 0., abdev = 0., slope = 0.;
    long    j = 0, i, istart;
    
    xwk = allocate_memory_double(ts->nseg+1, __LINE__); 
    ywk = allocate_memory_double(ts->nseg+1, __LINE__); 
    
    for (i = 1; i <= ts->n50; i++) {
        
        /* copy data of i'th segment into workspace */
        istart = (int)((i - 1) * (double)ts->nseg / 2.);
        for (j = 1; j <= ts->nseg; j++) {
            xwk[j] = ts->x[istart + j];
            ywk[j] = ts->y[istart + j];
        }
    
        /* detrend the data */
        if (c->DETREND == TRUE) {
            rmtrend(c, ts, ywk, xwk);
            for (j = 1; j <= ts->nseg; j++) {
                xwk[j] = ts->x[istart + j];
                ywk[j] = ts->y[istart + j];
            }
        } else if (c->ROBUSTDETREND == TRUE) {
            /* 
                fitting method that is more sensitive to outliers in 
                the time-series, I guess the way to go? 
            */
            medfit(ts->nseg, ywk, xwk, &intercept, &slope, &abdev); 
            
            /* Subtract the fitted line from the time-series to obtain a linear detrend */
            for (j = 1; j <= ts->nseg; j++) {
                xfit    = (slope * xwk[j]) + intercept;
                ywk[j] -= xfit;
                xwk[j] = ts->x[istart + j];
                ywk[j] = ts->y[istart + j];
            }
        }       
    
        /* estimate and sum rho for each segment */
        tauest(c, ts, xwk, ywk, ts->nseg);
    
        /* bias correction for rho (Kendall & Stuart, 1967; Vol. 3)) */
        ts->rho = (ts->rho * ((double)ts->nseg - 1.0) + 1.0) / ((double)ts->nseg - 4.0);
        rhosum += ts->rho;
    }
    
    /* average rho */
    ts->rho = rhosum / (double)ts->n50;
    
    /* average tau */
    ts->tau = (-1. * ts->avgdt) / log(ts->rho);
    
    /* tidy up */
    free_array_double(xwk);   
    free_array_double(ywk);

    return;
}

void tauest(control *c, timeseries *ts, double *xwork, double *ywork, int np)
{
    /* 
        Manfred Mudelsee's code for tau estimation, func for persistence 
        estimation for unevenly spaced time seriesautocorrelation coefficient 
        estimation (equidistant data). Note I have added the option to calc
        the coefficent differently (i.e. the standard way) for when one is
        using "equally" spaced data.
    */
    int     j;
    double  amin = 0., scalt = 0., fac = 0., avg, var, *tscal, *xscal;
    
    tscal = allocate_memory_double(np+1, __LINE__); 
    xscal = allocate_memory_double(np+1, __LINE__); 
    
    /* correct time direction? */
    for (j = 1; j <= np; j++) {
        tscal[j] = -1.0 * xwork[np + 1 - j];
        xscal[j] = ywork[np + 1 - j];
        /*tscal[j] = xwork[j];
        xscal[j] = ywork[j];*/
    }
    
    /* scaling of y */
    avevar(&avg, &var, xscal, np);
    
    fac = sqrt(var);
    for (j = 1; j <= np; j++) {
        xscal[j] /= fac;
    }
    
    /* scaling of x */
    ts->dt = (tscal[np] - tscal[1]) / ((double)(np - 1.0));
    rhoest(xscal, c, ts);
    
    if (ts->rho <= RHOMIN)
        ts->rho = 0.05;
    else if (ts->rho > RHOMAX)
        ts->rho = 0.95;
    scalt = -log(ts->rho) / ts->dt;
    for (j = 1; j <= ts->nseg; j++) {
        tscal[j] *= scalt;
    }
    
    /*for (j = 1; j <= np; j++) {
          
        ts->tscal[j] =  tscal[j];
        ts->xscal[j] =  xscal[j];
    }*/
    
    /* estimation */
    amin = minimisation_wrapper(c, ts, tscal, xscal);
    
    /* determine tau */
    ts->tau = -1.0 / (scalt * log(amin));
            
    /* determine corresponding rho */
    ts->rho = exp(-ts->dt / ts->tau);
    
    /* tidy up */
    free_array_double(tscal);    
    free_array_double(xscal);
    
    return;
}    
        
void rhoest(double *tmp, control *c, timeseries *ts)
{
    /* autocorrelation coefficient estimation (equidistant data) */
    
    int   i;
    double sum1 = 0.0, sum2 = 0.0;
    
    for (i = 2; i <= ts->nseg; i++) {
        sum1 += tmp[i] * tmp[i - 1];
        sum2 += tmp[i] * tmp[i];
    }
    
    ts->rho = sum1 / sum2;
    
    return;
}

double minimisation_wrapper(control *c, timeseries *ts, double *tscal, double *xscal)
{
    /* 
        Note unlike the redfit code if I find an error with the
        minimisation I barf rather than continue
    */
    
    double a_ar1 = 0.367879441;                                 /* 1/e */
    double tol   = 1.0E-6;                                      /* Brent's search, precision */
    double tol2  = 1.0E-6;                                      /* multiple solutions, precision */
    double dum1  = 0.0, dum2 = 0.0, dum3 = 0.0, dum4 = 0.0, a_ar11 = 0.0, a_ar12 = 0.0, a_ar13 = 0.0;
    int    nmu   = 0;                                           /* flag to indicate multiple minima found */
    double amin;                                                /* estimated value of a = exp(-scalt/tau) */
     
    
    dum1 = brent(-2.0, a_ar1, 2.0, tol, c, ts, &a_ar11, tscal, xscal);
    dum2 = brent(a_ar1, 0.5 *(a_ar1 + 1.0), 2.0, tol, c, ts, &a_ar12, tscal, xscal);
    dum3 = brent(-2.0, 0.5 * (a_ar1 - 1.0),  a_ar1, tol, c, ts, &a_ar13, tscal, xscal);

    #ifdef debug
        fprintf(stderr, "%f %f %f ** %f %f %f\n", dum1,dum2,dum3, a_ar11,a_ar12,a_ar13);
        fflush(stderr);
    #endif

    if (ts->min_err == 1) {
        check_errors("Error in minimisation_wrapper (call to brent)", __LINE__);
    }
      
    if (a_ar11 == -1.0)
        a_ar11 = ts->rho;
    if (a_ar12 == -1.0)
        a_ar12 = ts->rho;
    if (a_ar13 == -1.0)
        a_ar13 = ts->rho;
    
    if ((fabs(a_ar12 - a_ar11) > tol2 && fabs(a_ar12 - a_ar1) > tol2) ||\
        (fabs(a_ar13 - a_ar11) > tol2 && fabs(a_ar13 - a_ar1) > tol2)) 
          nmu = 1;
  
    if (nmu == 1) {
        check_errors("Least squares solution has multiple minima", __LINE__);
    }
    
    if (dum1 < dum2 && dum1 < dum3)
        dum4 = dum1;
    else if (dum2 < dum1 && dum2 < dum3)
        dum4 = dum2;
    else if (dum3 < dum1 && dum3 < dum2)
        dum4 = dum3;  
        
    if (dum4 == dum2)
         amin = a_ar12;
    else if (dum4 == dum3)
         amin = a_ar13;
    else
         amin = a_ar11;
        
    if (amin <= 0.0) {
        check_errors("Estimation problem: AMIN =< 0.0", __LINE__);
    }
    
    if (amin >= 1.0) {
        check_errors("Estimation problem: AMIN >= 1.0", __LINE__);
    }        
    
    return (amin);
}

double least_sqs(double a, int end_loop, double *tscal, double *xscal)
{
    double ls = 0.0;
    int    i;
    
	for (i = 2; i <= end_loop; i++) {
		ls += pow(xscal[i] - xscal[i - 1] * SIGN(1.0, a) *\
                pow(fabs(a), (tscal[i] - tscal[i - 1])), 2.0);	

    }

    return ls;
} 

double brent(double ax, double bx, double cx, double tol, control *c, timeseries *ts, double *xmin, double *tscal, double *xscal)
{
    /* Modified form of Numerical Recipes Brent 1D min finder, pg.404, 10.2 */
    
    #define ITMAX 100
    #define CGOLD 0.3819660
    #define ZEPS  1.0e-10
    
    double a = 0.0, b = 0.0, d = 0.0, e = 0.0, etemp = 0.0, fu = 0.0, fv = 0.0, fw = 0.0, fx = 0.0; 
    double p = 0.0, q = 0.0, r = 0.0, tol1 = 0.0, tol2 = 0.0, u = 0.0, v = 0.0, w = 0.0, x = 0.0, xm = 0.0;
    int    iter;
    
    a  = (ax < cx ? ax : cx);
    b  = (ax > cx ? ax : cx);
    v  = bx;
    w  = v;
    x  = v;
    e  = 0.;
    fx = least_sqs(x, ts->nseg, tscal, xscal);
    fv = fx;
    fw = fx;
    
    /*fprintf(stderr, "** %f %f %f %f %f %f %f %f %f\n",a,b,v,w,x,e,fx,fv,fw);
    fflush(stderr); */
    
    for (iter = 1; iter <= ITMAX; iter++) {
        
        xm   = 0.5 * (a + b);
        tol1 = tol * fabs(x) + ZEPS;
        tol2 = 2.0 * tol1;
        /*fprintf(stderr, "** %.15f %.15f %.15f\n", xm, tol1, tol2);
        fflush(stderr);
        exit(1);*/
        if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
            *xmin = x;
            return (fx);
        }
        
        if (fabs(e) > tol1) {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
              
            if (q > 0.0) 
                p = -p;
              
            q     = fabs(q);
            etemp = e;
            e     = d;
              
            if (fabs(p) > fabs(0.5 * q * etemp) || p <= q * (a - x) ||\
                p >= q * (b - x)) {
                if (x >= xm)
                    e = a - x;
                else
                    e = b - x;
                d = CGOLD * e;
            } else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2) 
                    d = SIGN(tol1, xm - x);
            }
        } else {
            if (x >= xm)
                e = a - x;
            else
                e = b - x;
            d = CGOLD * e;
        }
        
        if (fabs(d) >= tol1)
           u = x + d;
        else
           u = x + SIGN(tol1, d);
          
        fu = least_sqs(u, ts->nseg, tscal, xscal);
        
        if (fu <= fx) {
            if (u >= x)
                a = x;
            else
                b = x;
            
            v  = w;
            fv = fw;
            w  = x;
            fw = fx;
            x  = u;
            fx = fu;
        } else {
            if (u < x)
                a = u;
            else
                b = u;
            
            if (fu <= fw || w == x) {
                v  = w;
                fv = fw;
                w  = u;
                fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v  = u;
                fv = fu;
            }
        }
        
    }
    ts->min_err = 1.0;
    *xmin = x;
    
    #undef ITMAX
    #undef CGOLD
    #undef ZEPS
    
    return (fx);         
}

double calculate_persistence(timeseries *ts) 
{
    double xbar = 0., arg1 = 0.0, arg2 = 0.0;
    long   j;
    
    for (j = 1; j <= ts->nseg; j++) {
        xbar +=  ts->y[j];
    }
    xbar /= ts->nseg;
    
    for (j = 1; j <= ts->nseg; j++) {
        arg1 += (ts->y[j] - xbar) * (ts->y[j+1] - xbar);
        arg2 += (ts->y[j] - xbar) * (ts->y[j] - xbar);
    }
    
    /* equal dist autocorrelation coefficient http://en.wikipedia.org/wiki/Autocorrelation 
        or in Mudelsee et al 2009 Nonlin, Processes Geophys */
    
    return (arg1 / arg2);
}
