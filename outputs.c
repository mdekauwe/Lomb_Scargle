#include "lombscargle.h"

void print_outputs(control *c, timeseries *ts, output_options *o)
{
    int k, start_offset = 0;
    
    if (c->PERIODOGRAM_TYPE == SCARGLE) { 
        start_offset = 2;
    } else if (c->PERIODOGRAM_TYPE == PRESS) {
        start_offset = 1;
    }
    
    
    if (c->DO_RUNS_TEST == TRUE) {
        printf("# runs = %d, critical btw %d and %d\n\n", ts->runs_count,\
                                                          ts->rcrithlo,  \
                                                          ts->rcrithi);
        /* so it only prints this once */
        c->DO_RUNS_TEST = FALSE;
        return;
    }
    
    if (o->OUTPUT_TYPE == FINDPEAKS) {
        window_bandwidth(c, ts);
        find_peaks(c, ts);
    } else if (o->OUTPUT_TYPE == PRINT_FREQ_SCALE) {
        for (k = start_offset; k <= ts->nout; k++) {
            printf("%.10lf %.10lf %.10lf %.10lf\n", ts->frq[k],   \
                                                    ts->gxxc[k],  \
                                                    ts->gredth[k],\
                                                    ts->rel_pow_rescale); 
        }
    } else if (o->OUTPUT_TYPE == PRINT_UNCORRECTED) {
        for (k = start_offset; k <= ts->nout; k++) {
            printf("%.10lf %.10lf %.10lf %.10lf\n", ts->frq[k],  \
                                                    ts->gxx[k],  \
                                                    ts->gxxc[k], \
                                                    ts->gredth[k]); 
        }
    } else if (o->OUTPUT_TYPE == PRINT_RELATIVE_POW) {
        /* As I am outputing 1/frequency - index k = 1, will be 0 and lead to divide by zero! */
        for (k = start_offset; k <= ts->nout; k++) {
            printf("%.10lf %.10lf %.10lf %.10lf\n", 1. / ts->frq[k],                          \
                                                    ts->gxxc[k] / ts->rel_pow_rescale,      \
                                                    ts->gredth[k] / ts->rel_pow_rescale,    \
                                                    ts->rel_pow_rescale); 
        }
    } else if (o->OUTPUT_TYPE == USE_db_SCALE) {
        
        for (k = start_offset; k <= ts->nout; k++) {
            printf("%.10lf %.10lf %.10lf %.10lf\n", ts->frq[k],   \
                                                    10. * log10(ts->gxxc[k]),  \
                                                    10. * log10(ts->gredth[k]),\
                                                    10. * log10(ts->rel_pow_rescale)); 
        }    
    } else if (o->OUTPUT_TYPE == DUMP_INPUT_DATA) {
        if (c->year_of_interest == 2006) {
            for (k = 1; k <= c->cnt2006; k++) {
                printf("%f %f\n", ts->x[k], ts->y[k]);
            }
        } else if (c->year_of_interest == 2007) {
            for (k = (c->cnt2006 + 1); k <= (c->cnt2006 + c->cnt2007); k++) {
                printf("%f %f\n", (ts->x[k]) - 142.0, ts->y[k]);
            }
        } else if (c->year_of_interest == 2008) {
            for (k = (c->cnt2006 + c->cnt2007 + 1); k <= (c->cnt2006 + c->cnt2007 + c->cnt2008); k++) {
                printf("%f %f\n", (ts->x[k]) - 284.0, ts->y[k]);
            }
        } else if (c->year_of_interest == 2009) {
            for (k = (c->cnt2006 + c->cnt2007 + c->cnt2008+ 1); k <= (c->cnt2006 + c->cnt2007 + c->cnt2008 + c->cnt2009); k++) {
                printf("%f %f\n", (ts->x[k]) - 426.0, ts->y[k]);
            }
        } if (c->year_of_interest == -99) {
            for (k = 1; k <= ts->num_param; k++) {
                printf("%f %f\n", (ts->x[k]), ts->y[k]);
            }
        } 
        
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
        
        exit(0);    
    } 
    
    
        
    return;
}

void find_peaks(control *c, timeseries *ts)
{
    /* 
        find the significant peaks...I am sure this could be formulated better, but it works 
        and frankly I can't be arsed to spend time on it.
    */
    double  currentmax = 0.;
    double  descent    = 0.;
    double *max = NULL, *freq = NULL, *pow = NULL, *sig = NULL;
    double  tmp_max    = 0.;
    int     real_start = 1;
    int     i = 0, j = 0, day = 0;
    int     peak_found = FALSE;
    int     end_offset = 0;
    
    if (c->PERIODOGRAM_TYPE == SCARGLE) { 
        end_offset = ts->nout;
    } else if (c->PERIODOGRAM_TYPE == PRESS) {
        end_offset = ts->nout+1;
    }
    
    /* set up some space */
    max  = allocate_memory_double(ts->nout+1, __LINE__);
    freq = allocate_memory_double(ts->nout+1, __LINE__);
    pow  = allocate_memory_double(ts->nout+1, __LINE__);
    sig  = allocate_memory_double(ts->nout+1, __LINE__);
    
    for (j = 1; j <= ts->nout; j++) {
        freq[j] = ts->frq[ts->nout + 1 - j];
        pow[j]  = ts->gxxc[ts->nout + 1 -j];
        sig[j]  = ts->gredth[ts->nout + 1 -j];
    }
    
    /* 
        for the scenario where we start with the power above the sign level
           ignore this info as there is not enough data to resolve a clear peak
           this may not be right, but appears logical to me 
    */
    if (sig[1] > pow[1]) {
        for (j = 1; j <= ts->nout; j++) {
            /* find out when the power is below the sign line - i.e. the real starting pt */
            if (pow[j] < sig[j]) {
                real_start = j;
                break;
            }
        }
    }
    
    for (i = real_start; i <= ts->nout; i++) {
        /* is the power is above the sign line we have found a peak    */
        if (pow[i] > sig[i]) {
            /* however we need to check this is a "real" peak, i.e. full resolved */
            for (j = i; j <= ts->nout; j++) {
                /* check if power goes below sign line again - indicating a peak */
                if (pow[j] < sig[j]) {
                    peak_found = TRUE;
                    descent    = j; 
                    break;
                }
            }
            
            /* if we have found a "real" peak - determine where the top point of this peak is */
            if (peak_found == TRUE) {
                currentmax = pow[real_start];
                for (j = i-1; j <= descent; j++) {
                    if (pow[j] > currentmax){
                        tmp_max = freq[j];
                        day = j;
                        currentmax = pow[j];    
                    }
                }
                
                max[day]   = tmp_max;
                i          = descent;
                currentmax = 0.0;
                peak_found = FALSE;
            }
        }    
    }    
    
    /* dump out sign. peaks*/
    for (i = 1; i < end_offset; i++) {
        if (max[i] > 0.0) {
            printf("%f %f %f %f %f %f\n", freq[i],               \
                                          pow[i],                \
                                          sig[i],             \
                                          ts->rel_pow_rescale,\
                                          pow[i],             \
                                          ts->bandwidth / 2.);    /* divided by 2, see pg59, Weedon 2002. */
        } else {
            printf("%f %f %f %f\n", freq[i],   \
                                    pow[i],    \
                                    sig[i],    \
                                    ts->rel_pow_rescale);
        }
    }        
    
    /* tidy up */
    free_array_double(max);
    free_array_double(freq);
    free_array_double(pow);
    free_array_double(sig);
        
    return;
}

