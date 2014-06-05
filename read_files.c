#include "lombscargle.h"

void read_input_file(char **argv, control *c, timeseries *ts, output_options *o)
{
    int     doy = 1;
    int     count_2006 = 0, count_2007 = 0, count_2008 = 0, count_2009 = 0;
    int     array_doy = 1;
    int     frame;
    int     distance_counter = 0;
    int     len;
    int     still_within_day = 1;
    int     num_undef_days = 0;
    int     num_valid_days = 1;
    int     DODGY_2006_PIXEL = FALSE;
    int     DODGY_2007_PIXEL = FALSE;
    int     DODGY_2008_PIXEL = FALSE;
    int     DODGY_2009_PIXEL = FALSE;
    char   *str;
    char    line[STRING_LENGTH];
    char    stry[STRING_LENGTH];
    char    strm[STRING_LENGTH];
    char    strd[STRING_LENGTH];
    int     year;
    int     day;
    int     month;
    long    offset = 0;
    long    j = 1, i;    /* j = 1 and not 0 because num rec wants the arrays to go from 1 to n. */
    float  *tmp = NULL;
    double  rainf_rate = 0.;
    double  distance;
    
    
    if (c->INPUT_FILE_TYPE == ASCII) {
        /* work out how big the file is */
        while (fgets(line, STRING_LENGTH, c->ifp) != NULL) {
            /* work out how big the file is */
            j++;
        }
    
        /* allocate space for arrays */
        ts->x = allocate_memory_double(j+1, __LINE__);
        ts->y = allocate_memory_double(j+1, __LINE__);
    
        /* sort the file out */
        rewind(c->ifp);
    
        /* read in file */
        j = 1;
        while (fgets(line, STRING_LENGTH, c->ifp) != NULL) {
            if (sscanf(line, "%lf %lf", &(ts->x[j]), &(ts->y[j])) != 2) {
                fprintf(stderr, "%s: badly formatted input in file on line %d\n", *argv, (int)j + 1);
                exit(1);
            }
            j++;
        }
    
        /* take away 1 as counter will have been incremented too far */
        ts->num_param = (j - 1);
        close_file(c->ifp, __LINE__);

    } else if (c->INPUT_FILE_TYPE == BINARY_VERT_LENGTH_SCALES) {
        if (c->start_row == -99 || c->end_row == -99) {
            check_errors("You need to set start_row and end_row on the cmd line if using binary mode, length scales", __LINE__);
        }
        
        /* allocate space for arrays */
        ts->x = allocate_memory_double(c->num_rows+1, __LINE__);
        ts->y = allocate_memory_double(c->num_rows+1, __LINE__);
        tmp = allocate_memory_float((c->num_rows * c->num_cols)+1, __LINE__);
        
        /* read file into memory */
        read_file_float(c->ifp, tmp, c->num_cols * c->num_rows, argv, __LINE__);
        
        distance_counter = 1;
        distance = 1.0;
            
        for (i = c->start_row; i <= c->end_row; i++) {
            offset = i * c->num_cols + c->col;
                
            if (tmp[offset] > ts->undef_val) {
                ts->x[distance_counter] = distance;
                ts->y[distance_counter] = tmp[offset];
                distance_counter++;
            }
            distance += 3.0;    
        }
        
        close_file(c->ifp, __LINE__);
            
    } else if (c->INPUT_FILE_TYPE == BINARY_HORI_LENGTH_SCALES) {
    
        if (c->start_col == -99 || c->end_col == -99) {
            check_errors("You need to set start_col and end_col on the cmd line if using binary mode, H.length scales", __LINE__);
        }
        
        /* allocate space for arrays */
        ts->x = allocate_memory_double(c->num_cols+1, __LINE__);
        ts->y = allocate_memory_double(c->num_cols+1, __LINE__);
        tmp  = allocate_memory_float((c->num_rows * c->num_cols)+1, __LINE__);
        
        /* read file into memory */
        read_file_float(c->ifp, tmp, c->num_cols * c->num_rows, argv, __LINE__);
        
        distance_counter = 1;
        distance = 1.0;
            
        for (i = c->start_col; i <= c->end_col; i++) {
            offset = c->row * c->num_cols + i;
            
            if (tmp[offset] > ts->undef_val) {
                ts->x[distance_counter] = distance;
                ts->y[distance_counter] = tmp[offset];
                distance_counter++;
            }
            distance += 1.0;    
        }
        close_file(c->ifp, __LINE__);
    
    } else if (c->INPUT_FILE_TYPE == BINARY_STANDARD) {        
    
        /* work out how big the file is */
        while (fgets(line, STRING_LENGTH, c->ifp) != NULL) {
            j++;
        }
    
        /* allocate space for arrays */
        ts->x = allocate_memory_double(j+1, __LINE__);
        ts->y = allocate_memory_double(j+1, __LINE__);
        tmp   = allocate_memory_float(j+1, __LINE__);
    
        if (c->num_frames > 0) {
            /* increase space for arrays */
            ts->x = reallocate_memory_double(ts->x, 365, __LINE__);
            ts->y = reallocate_memory_double(ts->y, 365, __LINE__);
            tmp = reallocate_memory_float(tmp, c->num_frames, __LINE__);
        }
    
        /* sort the file out */
        rewind(c->ifp);
    
        /* read in the relevant pixels info */
        j   = 1;
        
        /* start frame to get a lag sum rf from 12:00 to 12:00, frame 5 see below */
        frame = 5;
        still_within_day = 1;
        while ((len = getline(line, MAXLINE)) > 0) {
            
            if (c->num_frames > 0) {
                
                /* SECRET OPTION FOR RAINFALL SPECTRA */
                c->infile_fp = open_file(line, "rb", __LINE__);
                
                /* 
                    Ignore the first frame as this is the last 3 hr slot from the previous year
                    (i.e. 2100-2400) - Phil checked this on the AMMA website. This value is used
                    for interpolation purposes only. So first frame I am interested in is the 5th 
                    frame which will be 12-15:00 and the last frame for the first day is 12, 09-12:00
                    FRAME  TIME
                      0    21-2400 
                      1    24-0300
                      2    03-0600
                      3    06-0900
                      4    09-1200
                     !--------------!
                     ! 5    12-1500 !
                     ! 6    15-1800 !
                     ! 7    18-2100 !  Period of interest...i.e. 8 x 3hr time slots
                     ! 8    21-2400 !
                     ! 9    24-0300 !
                     ! 10   03-0600 !
                     ! 11    06-0900 !
                     ! 12    09-1200 !
                     !--------------!
                */
                /* 
                    This scheme also means there won't be enough data for the final day...which is fine
                    I think as it will be output as 0. And I think this is a reasonable assumption as it
                    won't be raining?
                */
                while (frame < c->num_frames) {
                    
                    /* 
                        read the relevant pixels info into memory, using fseek
                        to jump across the frames (i.e. days)
                    */
                    offset = ((frame * c->num_cols * c->num_rows) + (c->row * c->num_cols + c->col)) * sizeof(float);
                    seek_around_file(c->infile_fp, offset, SEEK_SET, __LINE__);
                    read_file_float(c->infile_fp, &tmp[j], 1, argv,  __LINE__);
                    
                    /* 
                        units of almip/trimm data is kg/m^2 sec -> converting to mm/3hr 
                        
                        mass of water = 1000m^2 
                        depth (m) = mass / density
                        depth (mm) = mass / density * 1000.
                        d=1000*mass/1000
                        
                        1 sec is (60 * 60) hours * 3 (3 hrs)    
     
                    */
                    rainf_rate += tmp[j] * 1000. * (3. * 60.0 * 60.0) / 1000.;
                    
                    /* sum up hourly data to daily data  */
                    if (still_within_day == 8) {
                    
                        ts->x[array_doy] = (double)array_doy;
                        ts->y[array_doy] =  rainf_rate;
                        /* 
                            However, as we come out of this loop
                            it will be incremented by 1 (see below), so we set 
                            this to 0.
                        */
                        array_doy++;
                        rainf_rate = 0.;
                        still_within_day = 0 ;    /* needs to start at one,  but is incremented on the next line */
                    }
                    still_within_day++;
                    frame++;
                }
                
            
                /* tidy up */
                close_file(c->infile_fp, __LINE__);    
                
            } else {
                
                str = (char *)malloc((STRING_LENGTH) * sizeof(char));
                
                /* break the input file up so we know the year, month and day of the binary image */
                strncpy(stry, line, 4);
                year = atoi(stry);
                strncpy(strm, line+4, 2);
                month = atoi(strm);
                strncpy(strd, line+6, 2);
                day = atoi(strd);
            
                /* work out the doy of year? assumes data starts on May 22nd!!! */
                doy = date_to_doy(year, month, day);
            
                /* join the file info back up so we read in the correct file wherever it lives on the system */
                sprintf(str, "%s%s%s%s%s%s", c->dir_fn, c->file_fn, stry, strm, strd, c->file_ext);
            
                /* 
                    read the relevant pixels info into memory, using fseek
                    to jump across the frames (i.e. days)
                */
                offset = (c->row * c->num_cols + c->col) * sizeof(float);
                c->infile_fp = open_file(str, "rb", __LINE__);
                seek_around_file(c->infile_fp, offset, SEEK_SET, __LINE__);
                read_file_float(c->infile_fp, &tmp[j], 1, argv,  __LINE__);
                
                if (tmp[j] > ts->undef_val) {
                    ts->x[array_doy] = (double)doy;
                    ts->y[array_doy] = tmp[j];
                    /*printf("%f %f\n",ts->x[array_doy],ts->y[array_doy]);*/  
                    array_doy++;
                    
                    if (year == 2006) {
                        count_2006++;
                    } else if (year == 2007) {
                        count_2007++;
                    } else if (year == 2008) {
                        count_2008++;
                    } else if (year == 2009) {
                        count_2009++;
                    }
                    num_undef_days = 0;
                } else {
                    /*  There is an issue of continuous missing data particularly in the north of the 
                        domain for 2007. If we find a 'run' of missing days, then set a flag so we can
                        remove the entirety of this years data from the array. In which case we would 
                        only process 2006 and 2008, 2009 joined together */
                        
                    num_undef_days++;
                    if ((num_undef_days >= 50) && (year == 2006)) {
                        DODGY_2006_PIXEL = TRUE;
                    } else if ((num_undef_days >= 50) && (year == 2007)) {
                        DODGY_2007_PIXEL = TRUE;
                    } else if ((num_undef_days >= 50) && (year == 2008)) {
                        DODGY_2008_PIXEL = TRUE;
                    } else if ((num_undef_days >= 50) && (year == 2009)) {
                        DODGY_2009_PIXEL = TRUE;    
                    }
                }
                j++;
                
                /* tidy up */
                close_file(c->infile_fp, __LINE__);
                free(str);    
            }
        }
        
        /* store counts for dumpin input purposes */
        c->cnt2006 = count_2006;
        c->cnt2007 = count_2007;
        c->cnt2008 = count_2008;
        c->cnt2009 = count_2009;
        
        /* account for any prolonged data gaps, i.e. remove a years data if there is no EO data in the array */
        num_valid_days = (array_doy - 1);
        
        remove_missing_years_from_input_array(c, ts, DODGY_2006_PIXEL, DODGY_2007_PIXEL, DODGY_2008_PIXEL, \
                                              DODGY_2009_PIXEL, count_2006, count_2007, count_2008, count_2009, \
                                              &num_valid_days);
        
        if ((c->INPUT_FILE_TYPE == BINARY_VERT_LENGTH_SCALES) || (c->INPUT_FILE_TYPE == BINARY_HORI_LENGTH_SCALES)) {
            /* -1 as counter will have been incremented too far */
            ts->num_param = (distance_counter - 1);
        } else if ((DODGY_2006_PIXEL == TRUE) || (DODGY_2007_PIXEL == TRUE) || (DODGY_2008_PIXEL == TRUE) ||\
                   (DODGY_2009_PIXEL == TRUE)) {
            ts->num_param = num_valid_days;
        } else {
            /* -1 as counter will have been incremented too far */
            ts->num_param = (array_doy - 1);
        }
    
        /*     instead of tacking multiple years of data together and processing it that way, we may also want to form a single
            years time series by averaging the multiple years (i.e. so we are left with an array > 142 days */ 
        if (c->AVERAGE_MULTI_YEAR == TRUE) {
            average_multiyear_input_into_single_time_series(ts);
        }
    }
    
    /* die if there is no valid data...probably water? */
    if (ts->num_param == 0) {
        check_errors("No valid data, all masked/missing, perhaps we are in a river?", __LINE__);
    } else if (c->ENOUGH_DATA == TRUE) {
        /* don't run code if we deem there not to be a suitable number of data points */
        if (ts->num_param <= 50) 
            exit(1);
    }
    
    /* output input file */
    if (o->OUTPUT_TYPE == DUMP_INPUT_DATA) {
        print_outputs(c, ts, o);
    }
    
    /* tidy up */
    free_array_float(tmp);    
    
    return;
}

int date_to_doy(int year, int month, int day)
{                
    /* Given the input file year, month, day, figure out the corresponding day of the year */
    int yr_offset = 0, day_offset = 0;
    int days_in_may = 10;
    int days_in_jun = 30;
    int days_in_jul = 31;
    int days_in_aug = 31;
    int days_in_sep = 30;
    
    if (year == 2006) {
        yr_offset = 0;
    } else if (year == 2007) {
        yr_offset = 142;
    } else if (year == 2008) {
        yr_offset = 284;
    } else if (year == 2009) {
        yr_offset = 426;    
    }
    
    if (month == 5) {
        day_offset = (day - 22) + 1;
    } else if (month == 6) {
        day_offset = days_in_may + day;
    } else if (month == 7) {
        day_offset = days_in_may + days_in_jun + day;
    } else if (month == 8) {
        day_offset = days_in_may + days_in_jun + days_in_jul + day;
    } else if (month == 9) {
        day_offset = days_in_may + days_in_jun + days_in_jul + days_in_aug + day;
    } else if (month == 10) {
        day_offset = days_in_may + days_in_jun + days_in_jul + days_in_aug + days_in_sep + day;
    }
    
    return (yr_offset + day_offset);
}

void remove_missing_years_from_input_array(control *c, timeseries *ts, int DODGY_2006_PIXEL, int DODGY_2007_PIXEL, int DODGY_2008_PIXEL, \
                                           int DODGY_2009_PIXEL, int count_2006, int count_2007, int count_2008, int count_2009,\
                                           int *num_valid_days)
{
    /*  what it says on the tin. If we have a run of missing data (>50 days) from a year, e.g. 2007 in the north of the domain
        then remove the entire year from the array, and resize the arrays so it contains only 2006 and 2008. */
    double *tempx = NULL, *tempy = NULL;
    int     i;
    
    if (DODGY_2006_PIXEL == TRUE) {
        
        tempx = allocate_memory_double(count_2007 + count_2008 + count_2009 + 1, __LINE__);
        tempy = allocate_memory_double(count_2007 + count_2008 + count_2009 + 1, __LINE__);
        
        /* copy valid data */
        for (i = count_2007; i <= count_2009; i++) {
            tempx[i - count_2006] = ts->x[i]; 
            tempy[i - count_2006] = ts->y[i];
        }
        
        /* restore the valid data and get rid of the missing stuff */
        free_array_double(ts->x);    
        free_array_double(ts->y);

        ts->x = allocate_memory_double(count_2007 + count_2008 + count_2009 + 1, __LINE__);
        ts->y = allocate_memory_double(count_2007 + count_2008 + count_2009 + 1, __LINE__);
        
        for (i = 1; i <= count_2007 + count_2008 + count_2009; i++) {
            ts->x[i] = tempx[i]; 
            ts->y[i] = tempy[i];
        }
        
        /* tidy up */
        free_array_double(tempx);
        free_array_double(tempy);
        
        *num_valid_days -= count_2006;
        c->cnt2006 = 0;    
    } 
    
    if (DODGY_2007_PIXEL == TRUE) {
        
        tempx = allocate_memory_double(count_2006 + count_2008 + count_2009 + 1, __LINE__);
        tempy = allocate_memory_double(count_2006 + count_2008 + count_2009 + 1, __LINE__);
        
        /* copy valid data */
        for (i = 1; i <= count_2006; i++) {
            tempx[i] = ts->x[i]; 
            tempy[i] = ts->y[i];
        }
        
        for (i = (count_2006 + count_2007 + 1); i <= count_2006 + count_2007 + count_2008 + count_2009; i++) {
            tempx[i - count_2007] = ts->x[i]; 
            tempy[i - count_2007] = ts->y[i];
            
        }
        
        /* restore the valid data and get rid of the missing stuff */
        free_array_double(ts->x);    
        free_array_double(ts->y);
        
        ts->x = allocate_memory_double(count_2006 + count_2008 + count_2009 + 1, __LINE__);
        ts->y = allocate_memory_double(count_2006 + count_2008 + count_2009 + 1, __LINE__);
        
        for (i = 1; i <= count_2006 + count_2008 + count_2009; i++) {
            ts->x[i] = tempx[i]; 
            ts->y[i] = tempy[i];
        }
    
        /* tidy up */
        free_array_double(tempx);    
        free_array_double(tempy);
        
        *num_valid_days -= count_2007;    
        c->cnt2007 = 0;
    }     
    
    if (DODGY_2008_PIXEL == TRUE) {
        
        tempx = allocate_memory_double(count_2006 + count_2007 + count_2009 + 1, __LINE__);
        tempy = allocate_memory_double(count_2006 + count_2007 + count_2009 + 1, __LINE__);
        
        /* copy valid data */
        for (i = 1; i <= (count_2006 + count_2007); i++) {
            tempx[i] = ts->x[i]; 
            tempy[i] = ts->y[i];
        }
        
        for (i = (count_2006 + count_2007 + count_2008 + 1); i <= (count_2006 + count_2007 + count_2008 + count_2009); i++) {
            tempx[i - count_2008] = ts->x[i]; 
            tempy[i - count_2008] = ts->y[i];
        }
        
        /* restore the valid data and get rid of the missing stuff */
        free_array_double(ts->x);    
        free_array_double(ts->y);
        
        ts->x = allocate_memory_double(count_2006 + count_2007 + count_2009 + 1, __LINE__);
        ts->y = allocate_memory_double(count_2006 + count_2007 + count_2009 + 1, __LINE__);
        
        for (i = 1; i <= (count_2006 + count_2007 + count_2009); i++) {
            ts->x[i] = tempx[i]; 
            ts->y[i] = tempy[i];
        }
        
        /* tidy up */
        free_array_double(tempx);    
        free_array_double(tempy);
        
        *num_valid_days -= count_2008;    
        c->cnt2008 = 0;
    }
    
    if (DODGY_2009_PIXEL == TRUE) {
        
        tempx = allocate_memory_double(count_2006 + count_2007 + count_2008 + 1, __LINE__);
        tempy = allocate_memory_double(count_2006 + count_2007 + count_2008 + 1, __LINE__);
        
        /* copy valid data */
        for (i = 1; i <= (count_2006 + count_2007 + count_2008); i++) {
            tempx[i] = ts->x[i]; 
            tempy[i] = ts->y[i];
        }
        
        /* restore the valid data and get rid of the missing stuff */
        free_array_double(ts->x);    
        free_array_double(ts->y);
        
        ts->x = allocate_memory_double(count_2006 + count_2007 + count_2008 + 1, __LINE__);
        ts->y = allocate_memory_double(count_2006 + count_2007 + count_2008 + 1, __LINE__);
        
        for (i = 1; i <= (count_2006 + count_2007 + count_2008); i++) {
            ts->x[i] = tempx[i]; 
            ts->y[i] = tempy[i];
        }
        
        /* tidy up */
        free_array_double(tempx);    
        free_array_double(tempy);
        
        *num_valid_days -= count_2009;
        c->cnt2009 = 0;    
    }
    
    return;
}
        
void average_multiyear_input_into_single_time_series(timeseries *ts)
{
    /*  what it says on the tin. Average multiple years of data together to make a single time series (i.e. less than 142 days) */
    int    *count = NULL;
    double *tempx = NULL, *tempy = NULL;
    int     days_in_a_year = 142;
    int     size_of_new_array = 1;
    int     i, j, missing_day_offset_counter = 1;
    
    tempx = allocate_memory_double(days_in_a_year+1, __LINE__);
    tempy = allocate_memory_double(days_in_a_year+1, __LINE__);
    count = allocate_memory_int(days_in_a_year+1,    __LINE__);
    
    /* keep a running total of data */
    for (i = 1; i <= ts->num_param; i++) {
        for (j = 1; j <= days_in_a_year; j++) {
            if (j == (int)ts->x[i]) {
                tempx[j]  = (double)j;
                tempy[j] += ts->y[i];
                count[j]++;
                break;
            } 
            if ((j + days_in_a_year) == (int)ts->x[i]) {
                tempx[j]  = (double)j;
                tempy[j] += ts->y[i];
                count[j]++;
                break;
            } 
            if ((j + days_in_a_year * 2) == (int)ts->x[i]) {
                tempx[j]  = (double)j;
                tempy[j] += ts->y[i];
                count[j]++;
                break;
            } 
            if ((j + days_in_a_year * 3) == (int)ts->x[i]) {
                tempx[j]  = (double)j;
                tempy[j] += ts->y[i];
                count[j]++;
                break;
            } else if (j > ts->x[i]) {
                break;
            }
        }
    }
    
    /* tidy up */
    free_array_double(ts->x);    
    free_array_double(ts->y);
    
    for (j = 1; j <= days_in_a_year; j++) {
        if (count[j] >= 1) {
            size_of_new_array++;
        }
    }
    
    ts->x = allocate_memory_double(size_of_new_array+1, __LINE__);
    ts->y = allocate_memory_double(size_of_new_array+1, __LINE__);
    
    missing_day_offset_counter = 1;
    for (j = 1; j <= days_in_a_year; j++) {
        
        if (count[j] >= 1) {
            ts->x[missing_day_offset_counter] = tempx[j];
            ts->y[missing_day_offset_counter] = tempy[j] / count[j];
            
            missing_day_offset_counter++;
        }
    }
    
    /* change num_param according to array size (minus 1 as it will have been incremented one two many above */
    ts->num_param = size_of_new_array - 1;
    
    /* tidy up */
    free_array_double(tempx);    
    free_array_double(tempy);
    free_array_int(count);

    return;
    
}
