/*========================================================================================
    Lomb-Scargle algorithm    - Estimate the power spectrum from an unevenly sampled time 
                              series, with red noise bias correction or potential to test
                              against white noise hypothesis.
  ========================================================================================
    
    ++ READ THE MANPAGE, CHANGES.TEXT
    ++ Code depends on being able to link against my C library functions...
        
     Martin De Kauwe: 13th January, 2010                                                         
=========================================================================================*/

#include "lombscargle.h"

int main(int argc, char **argv)
{
    /* setup up space to pass stuff around in */
    control *c;
    if ((c = malloc(sizeof(control))) == NULL) {
        check_errors("Control structure: Not allocated enough memory", __LINE__);
     }
    
    timeseries *ts;
    if ((ts = malloc(sizeof(timeseries))) == NULL) {
        check_errors("Timeseries structure: Not allocated enough memory", __LINE__);
     }
    
    meta *m;
    if ((m = malloc(sizeof(meta))) == NULL) {
        check_errors("Meta structure: Not allocated enough memory", __LINE__);
     }
    
    output_options *o;
    if ((o = malloc(sizeof(output_options))) == NULL) {
        check_errors("Output_options structure: Not allocated enough memory", __LINE__);
     }
    
    /* setup initial junk */
    setup_initial_conditions(c, ts, o);
    
    /* send sweep along the command line */
    clparse(argc, argv, c, ts, o);
    
    /* check the input file is on the stdin? */
    check_stdin(__LINE__);
    
    /* read data in from stdin */
    read_input_file(argv, c, ts, o);
    
    /* setup space for working arrays and other misc stuff */
    allocate_working_environment(c, ts);
    
    /* calculate lomb periodgram for raw time-series */
    spectrum_wrapper(c, ts, ts->x, ts->y, ts->frq, ts->gxx);
    
    /* Calculate red noise bias and correct periodgram 
        - need to generate a seed ONCE */
    srand (time(NULL));
    rednoise_bias_wrapper(c, ts, o);
    
    /* dump out what you need */
    print_outputs(c, ts, o);
    
    /* tidy up */
    free_array_double(ts->gxx);
    free_array_double(ts->frq);
    free_array_double(ts->red);
    free_array_double(ts->grr);
    free_array_double(ts->grrsum);
    free_array_double(ts->gredth);
    free_array_double(ts->corr);
    free_array_double(ts->gxxc);
    free_array_double(ts->grravg);
    free_array_double(ts->x);
    free_array_double(ts->y);  
    free(c);
    free(ts);
    free(m);
    free(o);
    
    return (0);
}

void spectrum_wrapper(control *c, timeseries *ts, double *xdata, double *ydata, double *frequency, double *power)
{
    /* 
        wrapper func - to calculate periodgram for given time-series 
        also has the options to remove trends from the data and apply
        a window.
    */
    
    int     i, j, istart;
    double  intercept = 0., abdev = 0., slope = 0.;
    double *ftrx = NULL, *ftix = NULL, *xwk = NULL, *ywk = NULL, scal;
    
    ftrx = allocate_memory_double(ts->lfreq+1, __LINE__);
    ftix = allocate_memory_double(ts->lfreq+1, __LINE__);
    xwk  = allocate_memory_double(ts->nseg+1,  __LINE__); /* periodgram power */
    ywk  = allocate_memory_double(ts->nseg+1,  __LINE__); /* corresponding frequencies */
    
    /* scale autospectrum and setup frequency axis */
    scal = 2.0 / ((double)ts->n50 * (double)ts->nseg * ts->df * ts->ofac);
    ts->factor = scal;    /* store for coherency code */
             
    for (i = 1; i <= ts->n50; i++) {
        
        /* copy data of i'th segment into workspace */
        istart = (int)((i - 1) * (double)ts->nseg / 2.);
        for (j = 1; j <= ts->nseg; j++) {
            xwk[j] = xdata[istart + j];
              ywk[j] = ydata[istart + j];
            
        }
    
        /* detrend the data */
        if (c->DETREND == TRUE) {
            rmtrend(c, ts, ywk, xwk);
        } else if (c->ROBUSTDETREND == TRUE) {
            /* 
                fitting method that is more sensitive to outliers in 
                the time-series, I guess the way to go? 
            */
            medfit(ts->nseg, ywk, xwk, &intercept, &slope, &abdev);
    
            /* Subtract the fitted line from the time-series to obtain a linear detrend */
            for (j = 1; j <= ts->nseg; j++) {
                ywk[j] -= (slope * xwk[j]) + intercept;
            }
        }       
        
        /* apply window to data */
        if (c->WINDOW_TYPE != NO_WINDOW) {
            if (c->WINDOW_TYPE == COSINE_TAPER) {
                taper_timeseries(c, ts, ywk);    /* leads to some bizarre results if used in coherency */
            } else {
                window(c, ts, ywk, xwk);
            }
        }
    
        /* 
            calculate periodgram 
                - there are two options the original implementation in Scargle paper
                  or the Press et al method. They give identical results.
        */
        if (c->PERIODOGRAM_TYPE == SCARGLE) { 
            ft_uneven_data(xwk, ywk, ftrx, ftix, ts->nfreq, ts->nseg, ts->lfreq, ts->wz);
            
            /* sum raw spectra */
            for (j = 1; j <= ts->nout; j++) {
                power[j] += ftrx[j] * ftrx[j] + ftix[j] * ftix[j];
                
            }
        } else if (c->PERIODOGRAM_TYPE == PRESS) {
            period(xwk, ywk, ftrx, ftix, ts->nseg, ts->df, ts->nout);
            
            /* sum raw spectra */
            for (j = 1; j <= ts->nout; j++) {
                power[j] += ftrx[j] * ftrx[j] + ftix[j] * ftix[j];
            }
        }
        for (j = 1; j <= ts->nseg; j++) {
            xwk[j] = 0.0;
              ywk[j] = 0.0;
        }    
    }    
    
    /* rescale periodogram */
    if (c->PERIODOGRAM_TYPE == SCARGLE) { 
        for (i = 1; i <= ts->nout; i++) {
            power[i] *= scal;
            frequency[i] = (i - 1) * ts->df;
        }
    } else if (c->PERIODOGRAM_TYPE == PRESS) {
        for (i = 1; i <= ts->nout; i++) {
            power[i] *= scal;
            /* press algorithm starts at df, rather than 0.0 */
            frequency[i] = i * ts->df;
        }
    }
    
    /* tidy up */
    free_array_double(xwk);   
    free_array_double(ywk);
    free_array_double(ftrx);
    free_array_double(ftix);
    
    return;
}

void rednoise_bias_wrapper(control *c, timeseries *ts, output_options *o)
{
    /* 
        As mentioned already, spectra of time-series shows a decrease in
        spectral amplitude with increasing frequency - red-noise. By fitting
        a first-order autoregressive (AR1) process to the irregularly spaced time-series
        and transforming this to the frequency domain, the spectra can be 
        tested against the hypothesis that the time-series originates from an AR1
        process.
        
        See Schulz and Mudelsee (2002). 
    */
    int    i, k, rprv, rnow; 
    double sum = 0.0, fac90, fac, rho, rhosq, fnyq;

    /* estimate variance for raw periodgram - needed to correct data */
    for (i = 1; i <= ts->nout; i++) {
          sum += ts->gxx[i];
    }
    ts->var = ts->df * sum;
    
    if (c->EQUAL_DIST_DATA == TRUE) {
        ts->equalrho = calculate_persistence(ts);
    } else {
        /* Estimate tau - unless prescribed? */
        if (ts->rho < 0.0) {
            /* Estimate a value for tau */
            gettau(c, ts);
            /* can't have a negative tau */
            if (ts->tau < 0.0) {
                ts->tau = 0.0;
                fprintf(stderr, "Negative tau returned, tau forced to zero: %f\n", ts->tau);
            }    
        } else {
            ts->tau = -(ts->avgdt) / log(ts->rho);
        }
    }     
    
    /* 
        need to feed random numbers with a different seed - one off runs can simply use the clock.
        However, when running on multiple machines, there is potentially a scenerio where the same
        seed will be used, hence there is the option to feed the program with a seed on the cmd
        line. Ideally some bash func that generates a load of NEGATIVE ints.
    */
    if (c->CMD_LINE_SEED == FALSE) {
        ts->idum = -fabs(rand());
    } else {
        if (ts->seed > 0) {
            check_errors("Error: You must seed me with a negative int", __LINE__);
        } else {
            ts->idum = ts->seed;
        }
    }
    
    /* if debugging...*/
    /*ts->idum = -123456789;*/
    
    /* 
        As the lomb-scargle spectrum cannot be directly corrected (dependent on 
        the sampling intervals), a monte-Carlo ensemble, nsim AR1 time-series are generated.
        The difference between the average ensemble spectrum from the theoretical
        specturm is used as the correction.
    */
    
    /* gneratre AR(1) spectra */
    for (i = 1; i <= ts->nsim; i++) {
        
        /* set up AR1 junk */
        make_AR1(c, ts);
        
        /* calculate lomb periodgram for AR1 spectra */
        spectrum_wrapper(c, ts, ts->x, ts->red, ts->frq, ts->grr);
        
        /* scale and sum red-noise spectra */
        sum = 0.0;
        for (k = 1; k <= ts->nout; k++) {
            sum += ts->grr[k];
        }
        
        ts->varx = ts->df * sum;
        fac      = ts->var / ts->varx;
        for (k = 1; k <= ts->nout; k++) {
            ts->grrsum[k] += ts->grr[k] * fac;    
        }
        
        /* Need to zero grr array */
        for (k = 1; k <= ts->nout; k++) {
            ts->grr[k] = 0.0;
        }
            
    }
    
    /* average red-noise spectrum */
    sum = 0.0;
    for (k = 1; k <= ts->nout; k++) {
        ts->grravg[k] = ts->grrsum[k] / ((double)ts->nsim);
        sum += ts->grravg[k]; 
    }
    
    ts->varx = ts->df * sum;
    fac = ts->var / ts->varx;
    
    for (k = 1; k <= ts->nout; k++) {
        ts->grravg[k] *= fac;
    }
    
    /* determine lag-1 autocorrelation coeff, which decays exponetially as a func of time */
    if (c->EQUAL_DIST_DATA == TRUE) {
        rho   = ts->equalrho; 
        rhosq = rho * rho;
    } else {
        rho   = exp((-1.0 * ts->avgdt) / ts->tau);
        rhosq = rho * rho;
    }
    
    
    /* 
        If we want to test the spectrum against the background white noise
        spectrum instead then according to Mann and Lees (1996) Eq. 1, then setting
        rho to zero yields a white-noise process, I guess technically we also need to
        set rho in the AR1 stuff, but this seems fine as gredth is what we use to form
        the significance level 
    */
    if (c->SIGNIF_TEST == WHITE) {
        rho = 0.0;
    }
    
    /* 
        set theoretical spectrum (e.g., Mann and Lees, 1996, Eq. 4) 
         make area equal to that of the input time series
    */
    fnyq = ts->frq[ts->nout];    
           
      sum  = 0.0;
    for (k = 1; k <= ts->nout; k++) {
        /* power spectrum of an AR(1) process is given by ... */
        ts->gredth[k] = (1.0 - rhosq) / (1.0 - 2.0 * rho * cos(M_PI * ts->frq[k] / fnyq) + rhosq);
        sum += ts->gredth[k];         
      }
    
    ts->varx = ts->df * sum;
      fac = ts->var / ts->varx;
    
    /* 
        work out degrees of freedom and signif. The degrees of freedom affect the size of the 
        confidence levels. Greater dof = a smoother spectrum, the lower the variance of the estimates 
        and the narrower the confidence intervals and hence levels (Weedon 2003). 
    */
    getdof(c, ts);
    fac90 = getchi2(ts->dof, ts->alpha) / ts->dof;
    
    /* determine correction factor */
    for (k = 1; k <= ts->nout; k++) {
        ts->gredth[k] *= fac;
        ts->corr[k] = ts->grravg[k] / ts->gredth[k];
        ts->gxxc[k] = ts->gxx[k] / ts->corr[k];
        ts->gredth[k] *= fac90;    
        
        /* store up area under power spectrum to normalise by */
        ts->rel_pow_rescale += ts->gxxc[k];    
    }

    if (c->DO_RUNS_TEST == TRUE) {
        /* 
            Check it is appropriate to use an AR1 model to describe the time-series 
            Using a non-parametric test according to Bendat and Piersol, 1986, Random
            data, 2nd ed, wiley: New York, pg 95.
        */
        if ((c->WINDOW_TYPE != RECTANGULAR) || (ts->ofac >= 1.5) || (ts->n50 != 1)) {
            check_errors("You can only apply the runs test using a rectangular window, ofac = 1.0, and segs = 1", __LINE__);
        }
    
        ts->runs_count = 1;
        rprv = SIGN(1.0, ts->gxxc[1] - ts->gredth[1]);
 
        for (k = 1; k <= ts->nout; k++) {
            rnow = SIGN(1.0, ts->gxxc[k] - ts->gredth[k]);
            if (rnow != rprv) {
                ts->runs_count++;
            }
            rprv = rnow;
        }
    
        /* note I have added 0.5 at the end to make sure rounding when converting to int is right */
        ts->rcrithlo = (int)((pow((-0.79557086 + 1.0088719 * sqrt((double)(ts->nout/2))), 2.0)) + 0.5);
        ts->rcrithi  = (int)((pow(( 0.75751462 + 0.9955133 * sqrt((double)(ts->nout/2))), 2.0)) + 0.5);
        
        print_outputs(c, ts, o);
    }

    return;
}

void setup_initial_conditions(control *c, timeseries *ts, output_options *o)
{
    strcpy(c->versionstamp, "Lombscargle V2.0, Build Date: " __DATE__ " " __TIME__ ", Martin De Kauwe\n");
    c->DO_RUNS_TEST       =  FALSE;
    ts->factor            =  0.0;
    ts->undef_val         = -500.0;        /* undef values in lst are -999., so 500 seems a nice number! */
    c->AVERAGE_MULTI_YEAR =  FALSE;
    c->EQUAL_DIST_DATA    =  FALSE;
    c->AUTO_SET_SAMP_INT  =  FALSE;
    c->CMDAVGDT           =  FALSE;
    c->CMDSI              =  FALSE;
    c->WHAT_IS_NYQ        =  FALSE;
    c->CMD_LINE_SEED      =  FALSE;
    c->WINDOW_TYPE        =  HANNING;
    c->DETREND            =  TRUE;
    c->ROBUSTDETREND      =  FALSE;
    c->ENOUGH_DATA        =  FALSE;
    c->INPUT_FILE_TYPE    =  ASCII;
    c->SIGNIF_TEST        =  RED;
    c->PERIODOGRAM_TYPE   =  SCARGLE;    /* option to use scargle implementation rather than Press et al method */
    ts->num_param         =  0;
    ts->ofac              =  4.0;  /* oversampling frequency */
    ts->hifac             =  1.0;     /* Specify the scaling factor of the average Nyquist frequency */
    ts->min_interval      =  0.0;
    ts->rel_pow_rescale   =  0.0;   /* used to rescale power spectrum to be relative power */
    ts->min_err           =  0;
    ts->n50               =  1;
    c->ifp                =  stdin;
    ts->nsim              =  200;
    ts->alpha             =  CONF90;
    c->row                = -99;   /* set to be negative so we can check it is set on cmd line if code used in binary mode */ 
    c->start_row          = -99;   /* set to be negative so we can check it is set on cmd line if code used in binary mode */ 
    c->end_row            = -99;   /* set to be negative so we can check it is set on cmd line if code used in binary mode */ 
    c->start_col          = -99;   /* set to be negative so we can check it is set on cmd line if code used in binary mode */ 
    c->end_col            = -99;   /* set to be negative so we can check it is set on cmd line if code used in binary mode */ 
    ts->sampling_int      = -99;   /* set to be negative so we can check it is set on cmd line */ 
    c->num_frames         =  0;
    c->USE_SET_MAX_FREQ   =  FALSE;
    ts->rho               = -99.;
    c->year_of_interest   = -99;
    
    /* Set default input directory */
    strcpy(c->dir_fn,   "/users/eow/mgdk/DATA/REGRIDDED_LST/gridded_03/");
    strcpy(c->file_fn,  "lsta_daily_");
    strcpy(c->file_ext, ".gra");
    
    /* output options */
    o->OUTPUT_TYPE = PRINT_FREQ_SCALE;
    
    return;
}




void make_AR1(control *c, timeseries *ts)
{
    /* 
        set up AR1 time series. An AR spectrum models the data using a linear combination of the previous and 
        subsequent time steps (see Mann and Lees, 1996): 
        
                red[i] = rho[i] * y[i-1] + eps[i]     

           where
                                t[i] - t[i-1]
                rho[i] =  exp(--------------)
                                   tau

    
        and 
        
                 eps ~ NV(0, vareps). To ensure Var[red] = 1, we set

                                    2 * (t[i] - t[i-1])
                vareps = 1 -  exp(- -------------------).
                                      tau

        Stationarity of the generated AR(1) time series is assured by dropping
        the first N generated points.
    */
    
    int    k;
    double dt, sigma;
    
    ts->red[1] = gasdev(&ts->idum);
    
    for (k = 2; k <= ts->num_param; k++) {
        dt = ts->x[k] - ts->x[k - 1];
        sigma = sqrt(1.0 - exp(-2.0 * dt / ts->tau));
        if (c->EQUAL_DIST_DATA == TRUE) {
            ts->red[k] = exp(-dt / ts->equalrho) * ts->red[k - 1] + sigma * gasdev(&ts->idum);
        } else { 
            ts->red[k] = exp(-dt / ts->tau) * ts->red[k - 1] + sigma * gasdev(&ts->idum);
        }
    }
    
    return;
}

void allocate_working_environment(control *c, timeseries *ts)
{
    long   i;
    double tp, sumdt = 0.;
    
     ts->nseg = (int)(2. * ts->num_param / (ts->n50 + 1)); /* pts per segment */ 
    
    /* average sampling interval of entire time series */
    for (i = 2; i <= ts->num_param; i++) {
        sumdt += ts->x[i] - ts->x[i-1];
        
    }
    
    /* set sampling int on cmd line */
    if (c->CMDAVGDT == FALSE) {
        ts->avgdt = sumdt / ((double)(ts->num_param - 1.0));
    } 
    
    tp = ts->avgdt * ((double)ts->nseg); /* average period of a segment  */
    
    if (c->AUTO_SET_SAMP_INT == TRUE) {
        if ((ts->avgdt > 0.999999) && (ts->avgdt <= 1.00001)) {
            c->EQUAL_DIST_DATA = TRUE;
        }
    }
    
    if (c->CMDSI == FALSE) { 
        ts->df = 1.0 / (ts->ofac * tp); /* freq. spacing */
    }
    
    ts->wz = 2.0 * M_PI * ts->df; /* omega = 2*pi*f    */
      
    /* 
        the nyquist freq is not well defined for irregularly spaced data
        so use so use average samp interval 
    */
    if (c->USE_SET_MAX_FREQ == FALSE) {
        ts->fnyq = ts->hifac * 1.0 / (2.0 * ts->avgdt); /* average Nyquist freq. */
    } else {
        ts->hifac = ts->max_freq / ( 1.0 / (2.0 * ts->avgdt) );
        ts->fnyq  = ts->hifac * 1.0 / (2.0 * ts->avgdt); /* NYQ freq explictly set */
    }
      
    if (c->WHAT_IS_NYQ == TRUE) {
        printf("%f %f\n", ts->fnyq, ts->df);
        exit(0);
    }
    
    /* seem that the adding of one is to overcome some fortran rounding issue to an int, so isn't required! */
    /*ts->nfreq = (int)(ts->fnyq / ts->df + 1.); */
    ts->nfreq = (int)(ts->fnyq / ts->df); /* number of freq for transform f(1) = f0; f(nfreq) = fNyq */
    
    /* make sure that number of frequencies is even - maybe put this in???? */
    /*if (ts->nfreq % 2 != 0) {
        ts->nfreq--;
    }*/
    ts->lfreq = (int)(ts->nfreq * 2);
      
    /*ts->nout = 0.5 *ts->ofac *ts->hifac * np;*/
    ts->nout  = ts->nfreq;
    
    /* sort out some space */
    ts->gxx    = allocate_memory_double(ts->nout+1,      __LINE__); /* periodgram power */
    ts->frq    = allocate_memory_double(ts->nout+1,      __LINE__); /* corresponding frequencies */
    ts->red    = allocate_memory_double(ts->num_param+1, __LINE__);
    ts->grr    = allocate_memory_double(ts->nout+1,      __LINE__);
    ts->grrsum = allocate_memory_double(ts->nout+1,      __LINE__);
    ts->gredth = allocate_memory_double(ts->nout+1,      __LINE__);
    ts->corr   = allocate_memory_double(ts->nout+1,      __LINE__);
    ts->gxxc   = allocate_memory_double(ts->nout+1,      __LINE__);
    ts->grravg = allocate_memory_double(ts->nout+1,      __LINE__);

    return;
}

