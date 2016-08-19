EXEC   := 2LPTic
SRCS := main.c power.c allvars.c save.c read_param.c  read_glass.c  
OBJS := $(SRCS:.c=.o)
INCL := allvars.h proto.h  Makefile

#OPT   +=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                         # for a single DM species in the input file by interleaved by a half a grid spacing
#OPT   +=  -DMULTICOMPONENTGLASSFILE  # set this if the initial glass file contains multiple components
#OPT   +=  -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on
                                     # particle type
#OPT   +=  -DNO64BITID    # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC  # only switch this on if particles start from a glass (as opposed to grid)
#OPT += -DONLY_ZA # swith this on if you want ZA initial conditions (2LPT otherwise)

OPTIONS :=  $(OPT)

CC       :=  mpicc
OPTIMIZE :=  -O3 -Wall 
GSL_INCL :=  -I$(TACC_GSL_INC) -I$(TACC_GSL_INC)/gsl
GSL_LIBS :=  "-Wl,-rpath,$(TACC_GSL_LIB)" -L$(TACC_GSL_LIB) 
FFTW_INCL:= -I$(TACC_FFTW2_INC) 
FFTW_LIBS:= "-Wl,-rpath,$(TACC_FFTW2_LIB)" -L$(TACC_FFTW2_LIB) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw 
MPICHLIB := -lmpich

CFLAGS :=   $(OPTIONS)  $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL)
LIBS   :=   -lm  $(MPICHLIB)  $(FFTW_LIBS) -ldrfftw -ldfftw  $(GSL_LIBS)  -lgsl -lgslcblas

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(OPTIMIZE) $(LIBS) -o  $(EXEC)  

%.o: %.c $(INCL) 
	$(CC) $(OPTIMIZE) $(CFLAGS) -c $< -o $@

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



