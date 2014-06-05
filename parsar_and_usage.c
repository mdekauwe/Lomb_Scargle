#include "lombscargle.h"

void clparse(int argc, char **argv, control *c, timeseries *ts, output_options *o)
{
    int i;
    
    for (i = 1; i < argc; i++) {
        if (*argv[i] == '-') {
            if (!strncasecmp(argv[i],      "-of",      3)) ts->ofac               = atof(argv[++i]);
            else if (!strncasecmp(argv[i], "-hi",      3)) ts->hifac              = atof(argv[++i]);
            else if (!strncasecmp(argv[i], "-nsim",    5)) ts->nsim               = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-undef",   6)) ts->undef_val          = atof(argv[++i]);
            else if (!strncasecmp(argv[i], "-dir",     4)) strcpy(c->dir_fn, argv[++i]);
            else if (!strncasecmp(argv[i], "-fname",   6)) strcpy(c->file_fn, argv[++i]);
            else if (!strncasecmp(argv[i], "-fext",    5)) strcpy(c->file_ext, argv[++i]);
            else if (!strncasecmp(argv[i], "-si",      3)) {
                c->CMDSI = TRUE;
                ts->df = atof(argv[++i]);
            }
            else if (!strncasecmp(argv[i], "-ai",      3)) {
                c->CMDAVGDT = TRUE;
                ts->avgdt = atof(argv[++i]);
            }
            else if (!strncasecmp(argv[i], "-mf",      3)) {
                ts->max_freq = atof(argv[++i]);
                c->USE_SET_MAX_FREQ = TRUE;
            }
            else if (!strncasecmp(argv[i], "-nrows",   6)) c->num_rows           = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-ncols",   6)) c->num_cols           = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-nframes", 6)) c->num_frames         = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-r",       6)) c->row                = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-c",       3)) c->col                = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-strow",   6)) c->start_row          = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-enrow",   6)) c->end_row            = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-stcol",   6)) c->start_col          = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-encol",   6)) c->end_col            = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-year",    5)) c->year_of_interest   = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-my",      3)) c->AVERAGE_MULTI_YEAR = TRUE;
            else if (!strncasecmp(argv[i], "-checkin", 8)) c->ENOUGH_DATA        = TRUE;
            else if (!strncasecmp(argv[i], "-eq",      3)) c->EQUAL_DIST_DATA    = TRUE;
            else if (!strncasecmp(argv[i], "-auto",    3)) c->AUTO_SET_SAMP_INT  = TRUE;
            else if (!strncasecmp(argv[i], "-segs",    5)) ts->n50 = atoi(argv[++i]);
            else if (!strncasecmp(argv[i], "-seed",    5)) {
                c->CMD_LINE_SEED = TRUE;
                ts->seed          = atoi(argv[++i]);
            }
            
            else if (!strncasecmp(argv[i], "-runs",    5)) c->DO_RUNS_TEST = TRUE;
            else if (!strncasecmp(argv[i], "-nowin",   6)) c->WINDOW_TYPE  = NO_WINDOW;
            else if (!strncasecmp(argv[i], "-costap",  7)) c->WINDOW_TYPE  = COSINE_TAPER;
            else if (!strncasecmp(argv[i], "-rec",     4)) c->WINDOW_TYPE  = RECTANGULAR;
            else if (!strncasecmp(argv[i], "-han",     4)) c->WINDOW_TYPE  = HANNING;
            else if (!strncasecmp(argv[i], "-scargle", 8)) c->PERIODOGRAM_TYPE = SCARGLE;
            else if (!strncasecmp(argv[i], "-press",   6)) c->PERIODOGRAM_TYPE = PRESS;
            else if (!strncasecmp(argv[i], "-fconf99", 8)) ts->alpha           = CONF99;
            else if (!strncasecmp(argv[i], "-fconf95", 8)) ts->alpha           = CONF95;
            else if (!strncasecmp(argv[i], "-fconf90", 8)) ts->alpha           = CONF90;
            else if (!strncasecmp(argv[i], "-fconf80", 8)) ts->alpha           = CONF80;
            else if (!strncasecmp(argv[i], "-white",   6)) c->SIGNIF_TEST      = WHITE;
            else if (!strncasecmp(argv[i], "-notre",   6)) c->DETREND          = FALSE;
            else if (!strncasecmp(argv[i], "-robtre",  7)) {
                c->DETREND       = FALSE;
                c->ROBUSTDETREND = TRUE;
            }
            else if (!strncasecmp(argv[i], "-b",       2)) c->INPUT_FILE_TYPE = BINARY_STANDARD;
            else if (!strncasecmp(argv[i], "-vls",     4)) c->INPUT_FILE_TYPE = BINARY_VERT_LENGTH_SCALES;
            else if (!strncasecmp(argv[i], "-hls",     4)) c->INPUT_FILE_TYPE = BINARY_HORI_LENGTH_SCALES;
            else if (!strncasecmp(argv[i], "-fp",      3)) o->OUTPUT_TYPE = FINDPEAKS;
            else if (!strncasecmp(argv[i], "-printo",  7)) o->OUTPUT_TYPE = PRINT_UNCORRECTED;
            else if (!strncasecmp(argv[i], "-printr",  7)) o->OUTPUT_TYPE = PRINT_RELATIVE_POW;
            else if (!strncasecmp(argv[i], "-db",      3)) o->OUTPUT_TYPE = USE_db_SCALE;
            else if (!strncasecmp(argv[i], "-dumpin",  7)) o->OUTPUT_TYPE = DUMP_INPUT_DATA;
            else if (!strncasecmp(argv[i], "-ver",     4)) {
                fprintf(stderr, "%s", c->versionstamp);
                exit(1);
            }
            else if (!strncasecmp(argv[i], "-u",       2) || !strncasecmp(argv[i], "-h", 2)) { 
                usage(argv);
                exit(1);
            } else {
                fprintf(stderr,"%s: unknown argument on command line: %s\n", argv[0], argv[i]); 
                usage(argv);
                exit(1);
             }
         }
      }
    
    return;
}

void usage(char **argv)
{    
    fprintf(stderr, "\n========\n");
    fprintf(stderr, " USAGE:\n");
    fprintf(stderr, "========\n");
    fprintf(stderr, "\t%s [options] < file.txt > outfile.txt\n", argv[0]);
    fprintf(stderr, "\n\nExpected input file is either two columns of ascii data: time, sample.\n");
    fprintf(stderr, "OR an ascii txt file which has the locations of the relevant binary files\n");
    fprintf(stderr, "\nwhere the options are:\n\n");
    fprintf(stderr, "\n++General options:\n" );
    fprintf(stderr, "[-ver          \t Print the version/build number and date.]\n");
    fprintf(stderr, "[-b            \t Use binary mode - i.e. code expects as input a text file with the locations of binary files]\n");
    fprintf(stderr, "[-scargle      \t Generate periodogram using Scargle paper implementation]\n");
    fprintf(stderr, "[-press        \t Generate periodogram using Press et al implementation]\n");
    fprintf(stderr, "[-my           \t Average multiyears of input data to get a single timeseries and process spectra of this]\n");
    fprintf(stderr, "[-dir          \t Set path to binary files, default is /users/eow/mgdk/DATA/REGRIDDED_LST/gridded_03/]\n");
    fprintf(stderr, "[-fname        \t Set file prefix, default is lsta_daily_]\n");
    fprintf(stderr, "[-fext         \t Set file extension, default is .gra]\n");
    fprintf(stderr, "[-nrows int    \t number of rows in binary images]\n");
    fprintf(stderr, "[-ncols int    \t number of cols in binary images]\n");
    fprintf(stderr, "[-r     int    \t row you want processed]\n");
    fprintf(stderr, "[-c     int    \t col you want processed]\n");
    fprintf(stderr, "[-undef double \t Set undefined value in input file - LST = -999., default is -500.]\n");
    fprintf(stderr, "[-white        \t Test significance of generate power spectrum against hypothesis of random white noise]\n");
    fprintf(stderr, "[-runs         \t Apply runs test to data.]\n");
    fprintf(stderr, "[-si    double \t Explicitly set the interval between output spectra i.e. 0.001]\n");
    fprintf(stderr, "[-ai    double \t Explicitly set the avg sampling int]\n");
    fprintf(stderr, "[-of    double \t Oversampling factor, set to 4.0 (typical values >= 4.0)]\n");
    fprintf(stderr, "[-hi    double \t Specify the scaling factor of the average Nyquist frequency\n");
    fprintf(stderr, "[              \t (Hifac x avg.Nyquist freq = Highest frequency to be examined) set to 1.0 by default]\n");
    fprintf(stderr, "[-mf    double \t Explicitly set the maximum frequency sampled upto, i.e. the NYQ]\n");
    fprintf(stderr, "[-segs  int    \t Number of WOSA segments with a 50 percent overlap]\n");
    fprintf(stderr, "[-fp           \t Find the statistically significant peaks and output the bandwidth resolution, plot in gnuplot u  1:5:6 w xerr]\n");
    fprintf(stderr, "[-eq           \t Data is equal distance use standard method for calculating autocorr coef.]\n");
    fprintf(stderr, "[-auto         \t Automatically determine sampling interval and use -eq option if required.]\n");
    fprintf(stderr, "\n++Data tapering options:\n" );
    fprintf(stderr, "[-nowin        \t Don't apply a window to the raw time-series]\n");
    fprintf(stderr, "[-costap       \t Apply a split cosine bell tapering to the time series]\n");
    fprintf(stderr, "[-rec          \t Use a rectangular window on the raw time-series]\n");
    fprintf(stderr, "[-han          \t Use a Hanning window on the raw time-series]\n");
    fprintf(stderr, "\n++Trend options:\n" );
    fprintf(stderr, "[-notre        \t Don't detrend the raw time-series]\n");
    fprintf(stderr, "[-robtre       \t Detrend using more *robust* fitting a line by minimising absolute deviation method]\n");
    fprintf(stderr, "\n++AR1 options:\n" );
    fprintf(stderr, "[-nsim  double \t Number of AR1 simulations to generate, default is 200]\n");
    fprintf(stderr, "[-seed -ve int \t Seed the AR1 numbers from cmd line - note seed must be a negative int]\n");
    fprintf(stderr, "[-fconf99      \t Output significance at 99 per confidence levels]\n");
    fprintf(stderr, "[-fconf95      \t Output significance at 95 per confidence levels]\n");
    fprintf(stderr, "[-fconf90      \t Output significance at 90 per confidence levels]\n");
    fprintf(stderr, "[-fconf80      \t Output significance at 80 per confidence levels]\n");
    fprintf(stderr, "\n++Cross spectrum analysis options:\n" );
    fprintf(stderr, "[-xspec        \t output the information required for cross spec analysis btw 2 time series]\n");
    fprintf(stderr, "\n++Output options:\n" );
    fprintf(stderr, "[-dumpin       \t Dump input file out to stderr]\n");
    fprintf(stderr, "[-year     int \t Dump input file out to stderr, but determine year of interest and time axis runs 1-142.]\n");
    fprintf(stderr, "[-checkin      \t Check the input time series is at least 50 members long, else die]\n");
    fprintf(stderr, "[-printo       \t Print corrected and uncorrected periodgrams]\n");
    fprintf(stderr, "[-printr       \t Print periodgrams, in relative power on the yaxis and period on xaxis]\n");
    fprintf(stderr, "[-db           \t Print periodgrams, in decibel scale]\n");
    fprintf(stderr, "\n++Print this message:\n" );
    fprintf(stderr, "[-u/-h         \t usage/help]\n");
    
    return;
}



