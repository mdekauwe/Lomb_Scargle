#ifndef LOMBSCARGLE
#define LOMBSCARGLE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>

/* My standard C library */
#include <memory.h>
#include <file_io.h>
#include <error_handling.h>

/*static long lminarg1, lminarg2;
static long lmaxarg1, lmaxarg2;*/
static double sqrarg;
static double tempr;

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define LMAX(a,b) (lmaxarg1 = (a),lmaxarg2 = (b), (lmaxarg1 > lmaxarg2 ? lmaxarg1 : lmaxarg2))
#define LMIN(a,b) (lminarg1 = (a),lminarg2 = (b), (lminarg1 < lminarg2 ? lminarg1 : lminarg2))
#define MOD(a,b)    while (a >= b) a -= b
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)
#define STRING_LENGTH 2000
#define MAXLINE 1000
const double TWO_M_PI = 6.2831853071795865;
/*#define TWO_M_PI 6.2831853071795865*/

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define SCARGLE 0
#define PRESS   1

#define NO_WINDOW    0
#define RECTANGULAR  1 
#define HANNING      2
#define COSINE_TAPER 3

#define CONF99 0.01
#define CONF95 0.05
#define CONF90 0.10
#define CONF80 0.20

#define RHOMIN 0.0
#define RHOMAX 1.0
 
#define FINDPEAKS          0 
#define PRINT_FREQ_SCALE   1
#define PRINT_UNCORRECTED  2
#define PRINT_RELATIVE_POW 3 
#define USE_db_SCALE       4 
#define DUMP_INPUT_DATA    5

#define ASCII                     0
#define BINARY_VERT_LENGTH_SCALES 1
#define BINARY_HORI_LENGTH_SCALES 2
#define BINARY_STANDARD           4

#define WHITE 0
#define RED   1

   
typedef struct  {
    long    num_param;
    double  equalrho;
    double *x;                /* array of times */
    double *y;                /* array of corresponding counts */
    double *ac;
    double *bc;
    double *wtau;
    double *tcos;
    double *tsin;
    double *gxx;
    double *frq;
    double *xwk;
    double *ywk;
    double *red;   
    double *freq;  
    double *grr;  
    double *grrsum;
    double *gredth;
    double *corr;
    double *gxxc;  
    double *grravg;
    double *white;
    double *tmp;
    double *cospectrum;
    double *quadrature;
    double  var;
    double  undef_val;
    double  tau;
    double  rho;
    int     n50;
    double  varx;
    int     nsim;
    int     runs_count;
    int     rcrithlo;
    int     rcrithi;
    double  alpha;
    double  dof;
    int     nseg;
    double  signif_lev;
    int     nfreq;
    int     lfreq;
    int     nout;
    int     min_err;
    double  normalise;
    double  xmin;
    long    idum;
    long    seed;
    double  dt;
    double  avgdt;
    double  df;
    double  wz;
    double  fnyq;
    double  sampling_int;
    double  bandwidth;
    double  avg_interval;
    double  min_interval;
    double  band_width;
    double  factor;
    double  prob;             /* false probability of the largest value of the lomb periodgram */
    double  max_freq;         /* set maximum frequency explicitly one wishes to sample upto */
    double  ofac;                /* oversampling factor, the value comes from NR in C */
    double  hifac;             /* Sets how high in frequency to go.*/
    double  rel_pow_rescale;
    long    jmax;                /* max of power */
} timeseries;

typedef struct  {
    FILE   *ifp;    
    FILE   *infile_fp;
    char    ifp_fn[STRING_LENGTH];
    char    dir_fn[STRING_LENGTH];
    char    file_fn[STRING_LENGTH];
    char    file_ext[STRING_LENGTH];
    char    versionstamp[STRING_LENGTH];
    int     WHAT_IS_NYQ;
    int     start_row;
    int     end_row;
    int     start_col;
    int     end_col;
    int     row;
    int     col;
    int     num_rows;
    int     num_cols;
    int     num_frames;
    int     NEG_COSPEC;
    int        NEG_QUADRA;
    int     STORE;
    int     CMDAVGDT;
    int     CMDSI;
    int     INPUT_FILE_TYPE;
    int     ENOUGH_DATA;
    int     EQUAL_DIST_DATA;
    int     OUTPUT_XSPEC_INFO;
    int     ROBUSTDETREND;
    int     DETREND;
    int     GARY;
    int     USE_SET_MAX_FREQ;
    int     SIGNIF_TEST;
    int     PERIODOGRAM;
    int        CMD_LINE_SEED;
    int     SMOOTH;
    int     WINDOW_TYPE;
    int     AVERAGE_MULTI_YEAR;
    double  mina;
} control;

typedef struct  {
    int     OUTPUT_TYPE;
} output_options;


typedef struct  {
    timeseries *time_series;
    control *controls;
    double min;
} meta;


/* generic misc funcs */
void   read_input_file(char **, control *, timeseries *, output_options *);
void   old_read_input_binaryfile(char **, control *, timeseries *);
void   setup_initial_conditions( control *, timeseries *, output_options *);
void   usage( );
void   clparse( int, char **, control *, timeseries *, output_options *);
void   allocate_working_environment(control *, timeseries *);
void   print_outputs(control *c, timeseries *, output_options *);
int    date_to_doy(int, int, int);
void   remove_missing_years_from_input_array(timeseries *, int, int, int, int, int, int, int *);
void   average_multiyear_input_into_single_time_series(timeseries *);

/* periodgram funcs */
void   ft_uneven_data(control *, timeseries *, double *, double *, double *, double *);
void   gettau(control *, timeseries *);
void   rhoest(double *, control *, timeseries *);
double least_sqs(double, int, double *, double *);
double minimisation_wrapper(control *, timeseries *, double *, double *);
void   rmtrend(control *, timeseries *, double *, double *);
void   window(control *, timeseries *, double *, double *);
void   taper_timeseries(control *, timeseries *, double *);
void   make_AR1(control *, timeseries *);
double getchi2(double, double);
void   spectrum_wrapper(control *, timeseries *, double *, double *, double *, double *);
void   rednoise_bias_wrapper(control *, timeseries *);
void   getdof(control *, timeseries *);
double getz(double alpha);
void   generate_white_noise(control *, timeseries *);
void   smooth_periodogram(double *, control *, timeseries *);
void   avg_and_var(control *, timeseries *, double *, double *, double *);
void   white_sig_test(control *, timeseries *);
void   window_bandwidth(control *, timeseries *);
void   find_peaks(control *c, timeseries *);
void   tauest(control *, timeseries *, double *, double *, int);
double calculate_persistence(timeseries *); 

/* numerical recipes funcs */
void   period(control *, timeseries *, double *, double *, double *, double *, double *);
double brent(double , double, double, double, control *, timeseries *, double *, double *, double *);


/* gsl */
int    gsl_brent(control *, timeseries *, double, double, double );
double mini_gsl_wrapper(double, void *);
double fn1 (double , void *);
#endif /* lombscargle */
