#include <drfftw_mpi.h>

#include <stdint.h>
#define  PI          3.14159265358979323846 
#define  GRAVITY     6.67408e-8
#define  HUBBLE      3.2407789e-18   /* in h/sec */


#if 0
double PowerSpec(double kmag);
double GrowthFactor(double astart, double aend);
double F_Omega(double a);
double F2_Omega(double a);
int    read_parameter_file(char *fname);
double PowerSpec_EH(double k);
double PowerSpec_Efstathiou(double k);

#ifdef T3E
typedef short int int4byte;	/* Note: int has 8 Bytes on the T3E ! */
typedef unsigned short int uint4byte;	/* Note: int has 8 Bytes on the T3E ! */
#else
typedef int int4byte;
typedef unsigned int uint4byte;
#endif



extern struct io_header_1
{
  uint4byte npart[6];      /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
  double mass[6];          /*!< mass[1] gives the particle mass */
  double time;             /*!< time (=cosmological scale factor) of snapshot */
  double redshift;         /*!< redshift of snapshot */
  int4byte flag_sfr;       /*!< flags whether star formation is used (not available in L-Gadget2) */
  int4byte flag_feedback;  /*!< flags whether feedback from star formation is included */
  uint4byte npartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
                                the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder. */
  int4byte flag_cooling;   /*!< flags whether radiative cooling is included */
  int4byte num_files;      /*!< determines the number of files that are used for a snapshot */
  double BoxSize;          /*!< Simulation box size (in code units) */
  double Omega0;           /*!< matter density */
  double OmegaLambda;      /*!< vacuum energy density */
  double HubbleParam;      /*!< little 'h' */
  int4byte flag_stellarage;     /*!< flags whether the age of newly formed stars is recorded and saved */
  int4byte flag_metals;         /*!< flags whether metal enrichment is included */
  int4byte hashtabsize;         /*!< gives the size of the hashtable belonging to this snapshot file */
  char fill[84];		/*!< fills to 256 Bytes */
}
header, header1;
#endif //commented out with #if 0

extern struct io_header_1
{ 
  int32_t npart[6];                        /*!< number of particles of each type in this file */
  double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
										 stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                         /*!< time of snapshot file */
  double redshift;                     /*!< redshift of snapshot file */
  int32_t flag_sfr;                        /*!< flags whether the simulation was including star formation */
  int32_t flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
  uint32_t npartTotal[6];          /*!< total number of particles of each type in this snapshot. This can be
									 different from npart if one is dealing with a multi-file snapshot. */
  int32_t flag_cooling;                    /*!< flags whether cooling was included  */
  int32_t num_files;                       /*!< number of files in multi-file snapshot */
  double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                       /*!< matter density in units of critical density */
  double OmegaLambda;                  /*!< cosmological constant parameter */
  double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int32_t flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
  int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
  uint32_t npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
  int32_t  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  char fill[60];                   /*!< fills to 256 Bytes */
} header, header1;



extern int      Nglass;
extern int      *Local_nx_table;
extern int      WhichSpectrum;

extern FILE     *FdTmp, *FdTmpInput;

extern int      Nmesh, Nsample;

extern int      SphereMode;

extern long long IDStart;

extern char     GlassFile[500]; 
extern char     FileWithInputSpectrum[500];

extern int      GlassTileFac,GlassTileFacSample; 

extern double   Box;
extern int Seed;

extern long long TotNumPart;

extern int      NumPart;

extern int      NTaskWithN;


extern int      *Slab_to_task;


extern struct part_data 
{
  float Pos[3];
  float Vel[3];
#ifdef  MULTICOMPONENTGLASSFILE                      
  int   Type;
#endif
  long long ID;
} *P;


extern double InitTime;
extern double Redshift;
extern double MassTable[6];


extern char OutputDir[100], FileBase[100];
extern int  NumFilesWrittenInParallel;


extern int      ThisTask, ThisTaskFileNumber, NTask;

extern int      Local_nx, Local_x_start;

extern int  IdStart;

extern unsigned int TotalSizePlusAdditional;
extern rfftwnd_mpi_plan Inverse_plan;
extern rfftwnd_mpi_plan Forward_plan;
//extern fftw_real        *Disp;
extern fftw_real        *Workspace;
//extern fftw_complex     *Cdata;


extern double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
extern double InputSpectrum_UnitLength_in_cm;
extern double G, Hubble;
extern double RhoCrit;

extern double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
extern double OmegaBaryon, HubbleParam;
extern double PrimordialIndex;
extern double ShapeGamma;

extern double Dplus; /* growth factor */


#ifdef DIFFERENT_TRANSFER_FUNC
extern int Type, MinType, MaxType;
#endif

extern int    WDM_On;
extern int    WDM_Vtherm_On;
extern double WDM_PartMass_in_kev;
