#include "lombscargle.h"

void rmtrend(control *c, timeseries *ts, double *residual, double *xwork)
{
    /* 
        figure out linear trend by means of least squares and subtract
        said trend from the data set - y = mx + b.
        See Numerical recipies pg 661-666 for explanation, although this is my
        formulation of the text.
    */
    double sumx, sumy, sumxy, sumxx, slope, intercept;
    int    i;
    
    sumx  = 0.;
    sumy  = 0.;
    sumxy = 0.;
    sumxx = 0.;
    
    for (i = 1; i <= ts->nseg; i++) {
        sumx  += xwork[i];
         sumy  += residual[i];
        sumxy += xwork[i] * residual[i];
        sumxx += xwork[i] * xwork[i];
    }
    
    slope     = (sumx * sumy - ts->nseg * sumxy) / (sumx * sumx - ts->nseg * sumxx);
    intercept = (sumy - slope * sumx) / ts->nseg;
    
    /* subtract least squares fit */
    for (i = 1; i <= ts->nseg; i++) {
        residual[i] -= (slope * xwork[i] + intercept);
      }
    
    return;
}

void window(control *c, timeseries *ts, double *tmp, double *xwork)
{
    /* 
        hanning window or hann (moderate res window) - used to reduce aliasing, or 
        false translation of power falling into some frequency range. Reduces 
        *leakage* allowing one to distinguish between spectral peacks, where the 
        leakage from a larger component may hide that of a smaller one. This effect can be
        reduced by either increasing the length of the time sequence or by using a windowing
        algorithm
    */
    
    double jeff, tlen, *ww, sumw2 = 0.0, scal;
    int    i;
    double fac3 = (double)ts->nseg - 1.0;
    
    /* set up some space */
    ww = allocate_memory_double(ts->nseg+1, __LINE__);
    tlen = xwork[ts->nseg] - xwork[1];
    
    for (i = 1; i <= ts->nseg; i++) {
        jeff = (((double)ts->nseg) * (xwork[i] - xwork[1])) / tlen;
        
        if (c->WINDOW_TYPE == RECTANGULAR)
            ww[i] = 1.0;
        else if (c->WINDOW_TYPE == HANNING)
            ww[i] = 0.5 * (1.0 - cos(TWO_M_PI * jeff / fac3));
        else {
            check_errors("Unknown window - check!", __LINE__);
        }
        sumw2 += ww[i] * ww[i]; 
    }
    
    scal = sqrt(((double)ts->nseg) / sumw2);
    
    for (i = 1; i <= ts->nseg; i++) {
        ww[i]  *= scal;
        tmp[i] *= ww[i];
      }
    
    /* tidy up */
    free_array_double(ww);
    
    return;
}

void window_bandwidth(control *c, timeseries *ts)
{
    /*
        Estimate the 6db bandwidth accounting for the OFAC
    */
    double bw = 0.0;
    
    if (c->WINDOW_TYPE == RECTANGULAR)
        bw = 1.21;
    else if (c->WINDOW_TYPE == HANNING)
        bw = 2.00;

    ts->bandwidth = ts->df * ts->ofac * bw;
    
    return;
}

void taper_timeseries(control *c, timeseries *ts, double *y)
{
    /* 
        Apply split-cosine-bell taper to time-series y, adapted
        from Bloomfield (1976) pg. 116. Prevent periodogram leakage
        due to the ending of a time series by tapering first 5% and 
        last 5% 
    */
    double stdiff = 0.0;
    double bit    = 0.0;
    double propor = 0.0;    /* proportion of time series to be tapered */
    double weight = 0.0;
    int   i       = 0;
    
    stdiff = ts->x[ts->nseg] - ts->x[1];
    if (stdiff > 0.0) {
        
        bit = stdiff / (10. * ts->nseg);
        
        for (i = 1; i <= ts->nseg; i++) {
            propor = (ts->x[ts->nseg] - ts->x[i] + bit) / stdiff;
            
            if ((propor < 0.05) || (propor > 0.95)) {
                weight = 0.5 - 0.5 * cos(M_PI * (propor / 0.05));
                y[i] *= weight;
            } 
        }
    } else {
        
        stdiff = ts->x[1] - ts->x[ts->nseg];
        bit = stdiff / (10. * ts->nseg);
        
        for (i = 1; i <= ts->nseg; i++) {
            propor = (ts->x[i] - ts->x[ts->nseg] + bit) / stdiff;
            
            if ((propor < 0.05) || (propor > 0.95)) {
                weight = 0.5 - 0.5 * cos(M_PI * (propor / 0.05));
                y[i] *= weight;
            } 
        }
    
    }        
    
    return;
}

