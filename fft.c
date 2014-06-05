#include "fft.h"

void ft_uneven_data(double *tmpxwk, double *tmpywk, double *ftrx, double *ftix, int nfreq, int np, int lfreq, double wz)
{
    /* 
        Estimation of periodgram for the time-series given that the data is likely
        to irregularly spaced! This implementation is based on
        Scargle, 1989, The Astrophysical Journal (a cracker), 343, 874-887. 
    */
    
    double  csum, ssum, wdel, wrun = 0., wuse,  tim, sumtc, sumts, ttt;
    double  tc, tss, watan, wtnew, arg;
      double  tol1 = 1.0E-04, tol2 = 1.0E-08;
    double *wtau = NULL, *tcos = NULL, *tsin = NULL;
    int     i, count, count_start_offset;
    long    offset;
    int     il, iput, nstop;
    double *cross = NULL, *scos2 = NULL, *ssin2 = NULL, *sumr = NULL, *sumi = NULL;
    double  sqrt2 = 1.41421356237309504880168872420969807856967;
    double  fnn, const2, const1 = 1.0 / sqrt2, ftrd, ftid;
    double  si = 1.0, sumt = 0., sumx = 0.;
    double  tzero = 0.0, phase;
    double  cmplx_real, cmplx_im;
    
    /* setup some space */
    wtau  = allocate_memory_double(nfreq+1,        __LINE__);
    tcos  = allocate_memory_double((np * nfreq)+1, __LINE__);
    tsin  = allocate_memory_double((np * nfreq)+1, __LINE__);
    cross = allocate_memory_double(nfreq+1,        __LINE__);
    scos2 = allocate_memory_double(nfreq+1,        __LINE__);
    ssin2 = allocate_memory_double(nfreq+1,        __LINE__);
    sumr  = allocate_memory_double(nfreq+1,        __LINE__);
    sumi  = allocate_memory_double(nfreq+1,        __LINE__);
    
    /* trig stuff */
    wuse  = wz;
    wdel  = wuse;
    wrun  = wuse;
    count = 2;
    count_start_offset = count;
    
    /* start frequency loop */
    while (count <= nfreq) { 
        
        /* calc. tau */
        csum  = 0.0;
        ssum  = 0.0;
        sumtc = 0.0;
        sumts = 0.0;
        
        for (i = 1; i <= np; i++) { 
            ttt    = tmpxwk[i];
            arg    = 2.0 * wrun * ttt;
            tc     = cos(arg);
            tss    = sin(arg);
            csum  += tc;
            ssum  += tss;
            sumtc += ttt * tc;
            sumts += ttt * tss;
            
        }
        
        if (fabs(ssum) > tol1 || fabs(csum) > tol1)
            watan = atan2(ssum, csum);
        else
            watan = atan2(-sumtc, sumts);
        
        wtau[count] = 0.5 * watan;
        wtnew = wtau[count];
        
        /* summations over the sample */
        for (i = 1; i <= np; i++) { 
            offset       = (count - count_start_offset) * np + i;
            tim          = tmpxwk[i];
            arg          = wrun * tim - wtnew;
            tcos[offset] = cos(arg);
            tsin[offset] = sin(arg);    
        }
        wrun += wdel;
        count++;
    }
    
    fnn    = (double)(np);
    const2 = si * const1;
      
    for (i = 1; i <= np; i++) { 
        sumt += tmpxwk[i];
        sumx += tmpywk[i];
    }

    /* initialise for zero frequency */
    ftrx[1] = sumx / sqrt(fnn);            /* real part of DFT */
    ftix[1] = 0.0;                        /* img part of DFT  */
    wdel    = wuse;
    wrun    = wuse;
    count   = 2;
    count_start_offset = count;
    
    /* start frequency loop */
    while  (count <= nfreq) { 
        wtnew = wtau[count];
        
        /* summation over the sample */
        cross[count] = scos2[count] = ssin2[count] = sumr[count] = sumi[count] = 0.0;
        
        for (i = 1; i <= np; i++) { 
            offset        = (count - count_start_offset) * np + i;
            cross[count] += tmpxwk[i] * tcos[offset] * tsin[offset];
            scos2[count] += tcos[offset] * tcos[offset];
            ssin2[count] += tsin[offset] * tsin[offset];
            sumr[count]  += tmpywk[i] * tcos[offset];
            sumi[count]  += tmpywk[i] * tsin[offset];
        }
        
        ftrd = const1 * sumr[count] / sqrt(scos2[count]);
        
        if (ssin2[count] <= tol1) { 
            ftid = const2 * sumx / sqrt(fnn);
            if (fabs(cross[count]) > tol2) 
                ftid = 0.0;
            } else
                ftid = const2 * sumi[count] / sqrt(ssin2[count]);
        
        phase = wtnew - wrun * tzero;
        
        /*
            first have to deal with complex exponential cexp in fortran 
            cexp(cmplx(0.0, phase))
            c.re = exp( z.re );
            c.im = c.re * sin( z.im );
            c.re *= cos( z.im );
        
        */
        cmplx_real  = exp(0.0);
        cmplx_im    = cmplx_real * sin(phase);
        cmplx_real *= cos(phase);
        
        /* 
            ok the next bit is a hack to get around using the complex arrays used 
            in the redfit code.
            
            multiplication of complex numbers =
            z3 = z1 x z2
            real part = (real_1 * real_2) - (imag_1 * imag_2) 
            imag part = (real_1 * imag_2) + (real_2 * imag_1)
            
            from Lennox and chadwick (1977) Mathematics for engineers and applied
            scientists pg. 248
        */
        
        ftrx[count] = (ftrd * cmplx_real) - (ftid * cmplx_im);
        ftix[count] = (ftrd * cmplx_im)   + (cmplx_real * ftid);
        
        wrun += wdel;
        count++;
    }
    
    /* zero-fill transform (oversample inverse) impose symmetry for real data */
    if (2.0 * nfreq > lfreq) {
        check_errors("Error: 2 * nfreq > lfreq", __LINE__);
    }
      
    il = nfreq + 1;    /* gives good results on even data */
    for (i = il; i <= lfreq; i++) {
        ftrx[i] = 0.0;
        ftix[i] = 0.0;
    }
      
    nstop = (int)(lfreq / 2.);
    for (i = 2; i <= nstop; i++) {
        iput       = lfreq - i + 2;
        ftrx[iput] =  ftrx[i];
        ftix[iput] = -ftix[i];
    }
    
    
    /* tidy up */
    free_array_double(wtau);   
    free_array_double(cross);  
    free_array_double(scos2);  
    free_array_double(ssin2);  
    free_array_double(sumr);   
    free_array_double(sumi);   
    free_array_double(tsin);   
    free_array_double(tcos);   
    
    return;
}

void period(double *tmpxwk, double *tmpywk, double *ftrx, double *ftix, double np, double df, int nout)
{
    /* 
        Estimation of periodgram for the time-series given that the data is likely
        to irregularly spaced! This implementation is based on
        Press et al 1992, chp 13.8, pg 579. 
    */
    
    int    i, j;
    double avg, var;
    double cd, cc, cwtau, pnow, s, ss, sumc, sumcy, sums, sumsh;
    double sumsy, swtau, wtau, xave, xdif, xmax, xmin, yy;
    double arg, wtemp, *wi = NULL, *wpi = NULL, *wpr = NULL, *wr = NULL;
    
    /* allocate space for arrays */
    wi  = allocate_memory_double(np+1, __LINE__);
    wpi = allocate_memory_double(np+1, __LINE__);
    wpr = allocate_memory_double(np+1, __LINE__);
    wr  = allocate_memory_double(np+1, __LINE__);

    /* calculate the mean/variance of the time-series */
    avevar(&avg, &var, tmpywk, np);
    
    /* figure out range of x vector */
    xmax = xmin = tmpxwk[1];
    for (j = 1; j <= np; j++) {
        if (tmpxwk[j] > xmax) 
            xmax = tmpxwk[j];
        if (tmpxwk[j] < xmin) 
            xmin = tmpxwk[j];
    }
    xdif  = xmax - xmin;
    xave  = 0.5 * (xmax + xmin);
    
    /* decide on the first frequency uses avg sampling interval, based on Mudelsee et al*/
    pnow = df;                        
    /*pnow = 1.0 / (xdif * ofac); */ /* original scheme from Press et al */
    
    /* initialise values for the trigonometric recurrences at each data point */
    for (j = 1; j <= np; j++) {
        arg    =  TWO_M_PI * ((tmpxwk[j] - xave) * pnow);
        wpr[j] = -2.0 * SQR(sin(0.5 * arg));
        wpi[j] =  sin(arg);
        wr[j]  =  cos(arg);
        wi[j]  =  wpi[j];
    }
    
    /* loop over and evaluate all the frequencies */
    for (i = 1; i <= nout; i++) {
        sumsh = 0.0;
        sumc  = 0.0;
        
        /* first loop over the data to get tau and related quantites */
        for (j = 1; j <= np; j++) {
            cd     = wr[j];
            s      = wi[j];
            sumsh += s * cd;
            sumc  += (cd - s) * (cd + s);
        
        }
        wtau  = 0.5 * atan2(2.0 * sumsh, sumc);
        swtau = sin(wtau);
        cwtau = cos(wtau);
        sums  = 0.0;
        sumc  = 0.0;
        sumsy = 0.0;
        sumcy = 0.0;
        
        /* loop over the data again to get the periodgram values */
        for (j = 1; j <= np; j++) {
            s      = wi[j];
            cd     = wr[j];
            ss     = s  * cwtau - cd * swtau;
            cc     = cd * cwtau + s  * swtau;
            sums  += ss * ss;
            sumc  += cc * cc;
            yy     = tmpywk[j] - avg;
            sumsy += yy * ss;
            sumcy += yy * cc;
            
            /* update the trignometric recurrences */
            wtemp = wr[j];
            wr[j] = (wtemp * wpr[j] - wi[j] * wpi[j]) + wr[j];
            wi[j] = (wi[j] * wpr[j] + wtemp * wpi[j]) + wi[j];
            
        }
        ftrx[i] = sumcy / sumc;
        ftix[i] = sumsy / sums;
        
        /* move onto next frequency */
        pnow += df;
        /*pnow += 1.0 / (xdif * ofac); */ /* original scheme from Press et al */
    }        
    
    /* tidy up */
    free_array_double(wi);
    free_array_double(wpi);
    free_array_double(wpr);
    free_array_double(wr);
    
    return;
}

