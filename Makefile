EXEC   := 2LPTic
SRCS := main.c power.c allvars.c save.c read_param.c  read_glass.c  
OBJS := $(SRCS:.c=.o)
INCL := allvars.h save.h power.h read_glass.h read_param.h  Makefile


OPT += -DPRODUCE_CONSISTENT_IDS       # Set this to set exact same particle IDs as a higher resolution parent simulation (controlled by 
                                      # parameter GlassTileFacSample). Note, this might require 64bit particle IDs since the 
                                      # parent simulation might have required 64bit IDs (even though a contiguous set of IDs
                                      # for the current simulation might fit in 32bits)

OPT   +=  -DUSE_64BITID               # code outputs 64 bit particles IDs. If particle IDs don't fit inside a 
                                      # a signed int, and this 64bit ID option is not set, code will terminate

OPT += -DUSE_CAMB                     # Allow input tabulated transfer function to be in
                                      # CAMB format                                     

#OPT   +=  -DPRODUCEGAS               # Set this to automatically produce gas particles 
                                      # for a single DM species in the input file by interleaved by a half a grid spacing
#OPT   +=  -DMULTICOMPONENTGLASSFILE  # set this if the initial glass file contains multiple components

#OPT   +=  -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on
                                      # particle type. Switches Eisenstein & Hu for 2nd DM species

#OPT   +=  -DCORRECT_CIC              # only switch this on if particles start from a glass (as opposed to grid)

#OPT += -DONLY_ZA                     # switch this on if you want ZA initial conditions (2LPT otherwise)

#OPT += -DWDM_GAUSSIAN_VELOCITIES


OPTIONS :=  $(OPT)
CC       :=  mpicc
OPTIMIZE :=  -O3 -march=native

# Set for TACC Stampede.
GSL_INCL :=  -isystem$(TACC_GSL_INC) -isystem$(TACC_GSL_INC)/gsl
GSL_LIBS :=  "-Wl,-rpath,$(TACC_GSL_LIB)" -L$(TACC_GSL_LIB) 

FFTW_INCL:= -isystem$(TACC_FFTW2_INC) 
FFTW_LIBS:= "-Wl,-rpath,$(TACC_FFTW2_LIB)" -L$(TACC_FFTW2_LIB) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw 

MPICHLIB := -lmpich

CFLAGS := -std=gnu99 -fno-strict-aliasing
CFLAGS += -Wsign-compare -Wall -Wextra -Wshadow -Wunused
CFLAGS += -Wformat=2  -Wpacked  -Wnested-externs -Wpointer-arith  -Wredundant-decls  -Wfloat-equal -Wcast-qual
CFLAGS += -Wcast-align -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs -Wstrict-prototypes

LIBS   :=   -lm  $(MPICHLIB)  $(FFTW_LIBS) -ldrfftw -ldfftw  $(GSL_LIBS)  -lgsl -lgslcblas

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(OPTIMIZE) $(LIBS) -o  $(EXEC)  

%.o: %.c $(INCL) 
	$(CC) $(OPTIMIZE) $(CFLAGS) $(OPTIONS) $(FFTW_INCL) $(GSL_INCL) -c $< -o $@

.PHONY : clean celna clena celan
clean:
	rm -f $(OBJS) $(EXEC)
celan celna clena:clean


