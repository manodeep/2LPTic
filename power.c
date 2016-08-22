#include <stdlib.h>//for drand48
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "power.h"
#include "allvars.h"
#include "proto.h"


const size_t WORKSIZE = 100000;
#define MAXFILENAMELEN (1000)

const double EPS_REL_TOL = 1e-4;
const double EPS_ABS_TOL = 1e-6;

static double R8;
static double r_tophat;

static double AA, BB, CC;
static double nu;
static double Norm;


static int NPowerTable;

static struct pow_table
{
  double logk, logD;
}
 *PowerTable;

double fermi_dirac_kernel(double x,void *params);
void fermi_dirac_init(void);
double get_fermi_dirac_vel(void);
double get_gaussian_vel(void);
double PowerSpec_Efstathiou(double k);
double PowerSpec_EH(double k);
double PowerSpec_Tabulated(double k);
double PowerSpec_DM_2ndSpecies(double k);
double tk_eh(double k);
double growth(double a);
double growth_int(double a, void *param);
double sigma2_int(double k, void *param);
double TopHatSigma2(double R);
int    compare_logk(const void *a, const void *b);

double PowerSpec(double k)
{
  double power, alpha, Tf;

  switch (WhichSpectrum)
    {
    case 1:
      power = PowerSpec_EH(k);
      break;

    case 2:
      power = PowerSpec_Tabulated(k);
      break;

    default:
      power = PowerSpec_Efstathiou(k);
      break;
    }


  if(WDM_On == 1)
    {
      /* Eqn. (A9) in Bode, Ostriker & Turok (2001), assuming gX=1.5  */
      alpha =
	0.048 * pow((Omega - OmegaBaryon) / 0.4, 0.15) * pow(HubbleParam / 0.65,
							     1.3) * pow(1.0 / WDM_PartMass_in_kev, 1.15);
      Tf = pow(1 + pow(alpha * k * (3.085678e24 / UnitLength_in_cm), 2 * 1.2), -5.0 / 1.2);
      power *= Tf * Tf;
    }

#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)

  if(Type == 2)
    {
      power = PowerSpec_DM_2ndSpecies(k);
    }

#endif

  power *= pow(k, PrimordialIndex - 1.0);

  return power;
}


double PowerSpec_DM_2ndSpecies(double k)
{
  /* at the moment, we simply call the Eistenstein & Hu spectrum
   * for the second DM species, but this could be replaced with
   * something more physical, say for neutrinos
   */

  double power;
  const double tf_eh = tk_eh(k);
  power = Norm * tf_eh * tf_eh;

  return power;
}



void read_power_table(void)
{
  FILE *fd;
  char buf[MAXFILENAMELEN];
  double k, p;
#ifdef USE_CAMB        // Addition by Greg Poole (and Paul Geil) -- Use Camb transfer function
  double pcdm, pbar;   // Allow for individual transfer functions for CDM and Baryons
#endif
  double kmin,kmax;
  int nbytes=0;

  nbytes=snprintf(buf,MAXFILENAMELEN, "%s", FileWithInputSpectrum);
  if(nbytes > MAXFILENAMELEN){
	printf("Not enough storage allocated to print the name of the power-file. nbytes = %d MAXFILENAMELEN = %d .exiting\n",nbytes,MAXFILENAMELEN);
	FatalError(2283);
  }

  fd = fopen(buf,"r");
  if(fd == NULL)
    {
      printf("can't read input spectrum in file '%s' on task %d\n", buf, ThisTask);
      FatalError(17);
    }

  NPowerTable = 0;
  do
    {
#ifndef USE_CAMB   // Addition by Greg Poole (and Paul Geil) -- Use Camb transfer function
      if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
#else
    /* Format is
            k/h TF (CDM, Baryons, Photon, Massless Neutrino, Massive Neutrino, Total Mass)

       These quantities need to be converted to format 2LPT expects, which is Delta^2.

       So 
       
           Delta^2 = 4 * PI * k^3 * P_k = 4 * PI * k^3 * (8 * PI^2 *k * TF^2)

       and because this has been normalised, we just need to set p= TF^2 * k^4.
    */
		if(fscanf(fd,"%lg %lg %lg %lg %lg %lg %lg", &k, &pcdm, &pbar, &p, &p, &p, &p) == 7) 
#endif
		  {
			NPowerTable++;
		  } else {
		  break;
		}
    }
	  while(1);

  fclose(fd);


  if(ThisTask == 0)
    {
      printf("found %d pairs of values in input spectrum table\n", NPowerTable);
      fflush(stdout);
    }


  PowerTable = malloc(NPowerTable * sizeof(struct pow_table));
  if(PowerTable == NULL) {
	printf("Could not allocate memory for storing P(k) for %d entries on task %d\n", NPowerTable, ThisTask);
	FatalError(20);
  }
  snprintf(buf,MAXFILENAMELEN,"%s", FileWithInputSpectrum);
  fd = fopen(buf,"r");
  if(fd == NULL)
    {
      printf("can't read input spectrum in file '%s' on task %d\n", buf, ThisTask);
      FatalError(18);
    }

  NPowerTable = 0;


  kmin = 2 * M_PI / 0.0001;     /* 0.0001 h/Mpc */
  kmax = 2 * M_PI / 10000.0;    /* 10000 h/Mpc  */
  
  do
    {
#ifndef USE_CAMB   // Addition by CBP -- Use Camb transfer function
      if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
#else
    /* Format is
            k/h TF (CDM, Baryons, Photon, Massless Neutrino, Massive Neutrino, Total Mass)

       These quantities need to be converted to format 2LPT expects, which is Delta^2.

       So 
       
           Delta^2 = 4 * PI * k^3 * P_k = 4 * PI * k^3 * (8 * PI^2 *k * TF^2)

       and because this has been normalised, we just need to set p= TF^2 * k^4.
    */
	if(fscanf(fd,"%lg %lg %lg %lg %lg %lg %lg", &k, &pcdm, &pbar, &p, &p, &p, &p) == 7)
#endif
	{

#ifdef USE_CAMB
	  k=log10(k);
	  p=2.*log10(p)+4.*k;
#endif

	  PowerTable[NPowerTable].logk = k;
	  PowerTable[NPowerTable].logD = p;
	  NPowerTable++;
	  k = pow(10.0,k);
	  k /= (InputSpectrum_UnitLength_in_cm/3.085678e24); /* convert to h/Mpc*/
	  
	  if (k < kmin)
	    kmin = k;

	  if (k > kmax)
	    kmax = k;
	}
      else
	break;
    }
  while(1);
  fclose(fd);

  //check if there is sufficient k-coverage 
  double k_Nyquist = M_PI*Nsample/(Box*UnitLength_in_cm/3.085678e24);
  double k_fundamental = 2.0*M_PI/(Box*UnitLength_in_cm/3.085678e24);

  if(kmin > k_fundamental || kmax < k_Nyquist)
    {
      printf("[kmin, kmax] = [%lg,%lg] h/Mpc are not sufficient to cover [k_fund, k_ny] = [%lg,%lg] h/Mpc\n Aborting ...\n ",
	     kmin,kmax,k_fundamental,k_Nyquist);
      FatalError(19);
    }
  
  
  qsort(PowerTable, (size_t) NPowerTable, sizeof(struct pow_table), compare_logk);
}

int compare_logk(const void *a, const void *b)
{
  if(((const struct pow_table *) a)->logk < (((const struct pow_table *) b)->logk))
    return -1;

  if(((const struct pow_table *) a)->logk > (((const struct pow_table *) b)->logk))
    return +1;

  return 0;
}

void initialize_powerspectrum(void)
{
  double res;

  InitTime = 1 / (1 + Redshift);

  AA = 6.4 / ShapeGamma * (3.085678e24 / UnitLength_in_cm);
  BB = 3.0 / ShapeGamma * (3.085678e24 / UnitLength_in_cm);
  CC = 1.7 / ShapeGamma * (3.085678e24 / UnitLength_in_cm);
  nu = 1.13;

  R8 = 8 * (3.085678e24 / UnitLength_in_cm);	/* 8 Mpc/h */


  if(WhichSpectrum == 2)
    read_power_table();

#ifdef DIFFERENT_TRANSFER_FUNC
  Type = 1;
#endif

  Norm = 1.0;
  res = TopHatSigma2(R8);

  if(ThisTask == 0 && WhichSpectrum == 2)
    printf("\nNormalization of spectrum in file:  Sigma8 = %g\n", sqrt(res));

  Norm = Sigma8 * Sigma8 / res;

  if(ThisTask == 0 && WhichSpectrum == 2)
    printf("Normalization adjusted to  Sigma8=%g   (Normfac=%g)\n\n", Sigma8, Norm);

  Dplus = GrowthFactor(InitTime, 1.0);
}

double PowerSpec_Tabulated(double k)
{
  double logk, logD, kold, u, dlogk, Delta2;
  int binlow, binhigh, binmid;

  kold = k;

  k *= (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);	/* convert to desired units */

  logk = log10(k);

  if(logk < PowerTable[0].logk || logk > PowerTable[NPowerTable - 1].logk)
    return 0;

  binlow = 0;
  binhigh = NPowerTable - 1;

  while(binhigh - binlow > 1)
    {
      binmid = (binhigh + binlow) / 2;
      if(logk < PowerTable[binmid].logk)
	binhigh = binmid;
      else
	binlow = binmid;
    }

  dlogk = PowerTable[binhigh].logk - PowerTable[binlow].logk;

  if(dlogk <= 0.0)
    FatalError(777);

  u = (logk - PowerTable[binlow].logk) / dlogk;

  logD = (1 - u) * PowerTable[binlow].logD + u * PowerTable[binhigh].logD;

  Delta2 = pow(10.0, logD);

  return Norm * Delta2 / (4 * M_PI * kold * kold * kold);
}

double PowerSpec_Efstathiou(double k)
{
  return Norm * k / pow(1 + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);
}



double PowerSpec_EH(double k)	/* Eisenstein & Hu */
{
  const double tf_eh = tk_eh(k);
  return Norm * k * tf_eh * tf_eh;
}




double tk_eh(double k)		/* from Martin White */
{
  //MS 04/21/2011 - changed gamma to gamma1 mpicc was complaining about shadowed definition
  double q, theta, ommh2, a, s, gamma1, L0, C0;
  double tmp;
  double omegam, ombh2, hubble;

  /* other input parameters */
  hubble = HubbleParam;

  omegam = Omega;
  ombh2 = OmegaBaryon * HubbleParam * HubbleParam;

  if(OmegaBaryon <= 0.0)
    ombh2 = 0.04 * HubbleParam * HubbleParam;

  k *= (3.085678e24 / UnitLength_in_cm);	/* convert to h/Mpc */

  theta = 2.728 / 2.7;
  ommh2 = omegam * hubble * hubble;
  s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
  a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
    + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
  gamma1 = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
  gamma1 *= omegam * hubble;
  q = k * theta * theta / gamma1;
  L0 = log(2. * exp(1.) + 1.8 * q);
  C0 = 14.2 + 731. / (1. + 62.5 * q);
  tmp = L0 / (L0 + C0 * q * q);
  return (tmp);
}



double TopHatSigma2(double R)
{

  double result, abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;


  r_tophat = R;
  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &sigma2_int;

  /* note: 500/R is here chosen as (effectively) infinity integration boundary */
  gsl_integration_qag(&F, 0.0, 500.0 * 1 / R,
		      EPS_ABS_TOL, EPS_REL_TOL, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return result;
}


double sigma2_int(double k, void *param)
{
  (void) param;
  double kr, kr3, kr2, w, x;

  kr = r_tophat * k;
  kr2 = kr * kr;
  kr3 = kr2 * kr;

  if(kr < 1e-8)
    return 0.0;

  w = 3.0 * (sin(kr) / kr3 - cos(kr) / kr2);
  x = 4.0 * PI * k * k * w * w * PowerSpec(k);

  return x;
}

double GrowthFactor(double astart, double aend)
{
  return growth(aend) / growth(astart);
}


double growth(double a)
{
  double hubble_a;
  double result, abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;

  hubble_a = sqrt(Omega / (a * a * a) + (1 - Omega - OmegaLambda) / (a * a) + OmegaLambda);
  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  F.function = &growth_int;
  gsl_integration_qag(&F, 0.0, a, EPS_ABS_TOL, EPS_REL_TOL, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
  gsl_integration_workspace_free(workspace);

  return hubble_a * result;
}


double growth_int(double a, void *param)
{
  (void) param;
  return pow(a / (Omega + (1 - Omega - OmegaLambda) * a + OmegaLambda * a * a * a), 1.5);
}


double F_Omega(double a)
{
  double omega_a;

  omega_a = Omega / (Omega + a * (1 - Omega - OmegaLambda) + a * a * a * OmegaLambda);

  return pow(omega_a, 0.6);
}


double F2_Omega(double a)
{
  double omega_a;

  omega_a = Omega / (Omega + a * (1 - Omega - OmegaLambda) + a * a * a * OmegaLambda);

  return 2 * pow(omega_a, 4./7.);
}


/*  Here comes the stuff to compute the thermal WDM velocity distribution */


#define LENGTH_FERMI_DIRAC_TABLE 2000
#define MAX_FERMI_DIRAC          20.0

double fermi_dirac_vel[LENGTH_FERMI_DIRAC_TABLE];
double fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE];

double WDM_V0 = -1.0;


double fermi_dirac_kernel(double x,void *params)
{
  (void) params;
  return x * x / (exp(x) + 1);
}

void fermi_dirac_init(void)
{
  int i;
  double result, abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &fermi_dirac_kernel;
  
  for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
    {
      fermi_dirac_vel[i] = MAX_FERMI_DIRAC * i / (LENGTH_FERMI_DIRAC_TABLE - 1.0);
/*       fermi_dirac_cumprob[i] = qromb(fermi_dirac_kernel, 0.0, fermi_dirac_vel[i]); */
      gsl_integration_qag(&F, 0.0, fermi_dirac_vel[i], EPS_ABS_TOL, EPS_REL_TOL, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      fermi_dirac_cumprob[i] = result;
    }
  gsl_integration_workspace_free(workspace);
  
  for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
    fermi_dirac_cumprob[i] /= fermi_dirac_cumprob[LENGTH_FERMI_DIRAC_TABLE - 1];

  WDM_V0 = 0.012 * (1 + Redshift) * pow((Omega - OmegaBaryon) / 0.3, 1.0 / 3) * pow(HubbleParam / 0.65,
										    2.0 / 3) * pow(1.0 /
												   WDM_PartMass_in_kev,
												   4.0 / 3);

  if(ThisTask == 0)
    printf("\nWarm dark matter rms velocity dispersion at starting redshift = %g km/sec\n\n",
	   3.59714 * WDM_V0);

  WDM_V0 *= 1.0e5 / UnitVelocity_in_cm_per_s;

  /* convert from peculiar velocity to gadget's cosmological velocity */
  WDM_V0 *= sqrt(1 + Redshift);
}




double get_gaussian_vel(void)
{
  double x1,x2,w,y1;

  do {
	x1 = 2.0 *drand48() -1.0;
	x2 = 2.0 *drand48() -1.0;
	w = x1*x1 + x2*x2;
  } while (w >= 1.0);

  w = sqrt((-2.0 * log(w) )/w);
  y1 = x1 *w;
  /*double y2 = x2 *w; */

  return y1;
}

double get_fermi_dirac_vel(void)
{
  int i;
  double p, u;

  p = drand48();
  i = 0;

  while(i < LENGTH_FERMI_DIRAC_TABLE - 2)
    if(p > fermi_dirac_cumprob[i + 1])
      i++;
    else
      break;

  u = (p - fermi_dirac_cumprob[i]) / (fermi_dirac_cumprob[i + 1] - fermi_dirac_cumprob[i]);

  return fermi_dirac_vel[i] * (1 - u) + fermi_dirac_vel[i + 1] * u;
}

void add_WDM_thermal_speeds(float *vel)
{
  double v, phi, theta, vx, vy, vz;

  if(WDM_V0 <= 0.0)
    fermi_dirac_init();

#ifndef WDM_GAUSSIAN_VELOCITIES
  v = WDM_V0 * get_fermi_dirac_vel();
#else
  v = WDM_V0 * get_gaussian_vel();
#endif

  phi = 2 * M_PI * drand48();
  theta = acos(2 * drand48() - 1);

  vx = v * sin(theta) * cos(phi);
  vy = v * sin(theta) * sin(phi);
  vz = v * cos(theta);

  vel[0] += vx;
  vel[1] += vy;
  vel[2] += vz;
}
