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

/* gsl */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
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

typedef struct  {
    long    num_param;
    double *x;                /* array of times */
    double *y;                /* array of corresponding counts */
    double *ac;
    double *bc;
    double *tscal; 
    double *xscal;
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
    double  aa;
    double  abdevt;
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
    int     ENOUGH_DATA;
    int     OUTPUT_XSPEC_INFO;
    int     ROBUSTDETREND;
    int     DETREND;
    int     OUTPUT_dB;
    int     GARY;
    int     USE_BINARY_MODE;
    int     DUMP_INPUT_DATA;
    int     VERT_LENGTH_SCALES;
    int     HORI_LENGTH_SCALES;
    int     PRINT_UNCORRECTED;
    int     PRINT_RELATIVE_POW;
    int     USE_SET_MAX_FREQ;
    int     TEST_HYP_RED_NOISE;
    int     USE_SCARGLE;
    int        COSINE_TAPER;
    int        CMD_LINE_SEED;
    int     FINDPEAKS;
    int     SMOOTH;
    int     WINDOW;
    int     RECTANGULAR;
    int     HANNING;
    double  mina;
} control;

typedef struct  {
    timeseries *time_series;
    control *controls;
    double min;
} meta;



/* generic misc funcs */
void   read_input_asciifile(char **, control *, timeseries *);
void   read_input_binaryfile(char **, control *, timeseries *);
void   old_read_input_binaryfile(char **, control *, timeseries *);
void   setup_initial_conditions( control *, timeseries * );
void   usage( );
void   clparse( int, char **, control *, timeseries * );
void   allocate_working_environment(control *, timeseries *);
void   print_outputs(control *c, timeseries *);

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
/* numerical recipes funcs */
double avevar(control *, timeseries *, double *);
double brent(double , double, double, double, control *, timeseries *, double *, double *, double *);
double gasdev(control *, timeseries *);
double ran1(control *, timeseries *);
double gammp(double, double );
double gammln(double );
void   gcf(double *, double, double, double *);
void   gser(double *, double, double, double *);
void   period(control *, timeseries *, double *, double *, double *, double *, double *);
double erfcc(double x);
void   medfit(control *, timeseries *, double *, double *, double *, double *, double *);
double rofunc(double, control *, timeseries *, double *, double *);
double select_ele(unsigned long , control *, timeseries *);

/* gsl */
int gsl_brent(control *, timeseries *, double, double, double );
double mini_gsl_wrapper(double, void *);
double fn1 (double , void *);
#endif /* lombscargle */
