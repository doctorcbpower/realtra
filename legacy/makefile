#OPT += -DNFW_HALO
OPT += -DHERNQUIST_HALO
#OPT += -DLOGARITHMIC_HALO
#OPT += -DHERNQUIST_BULGE
OPT += -DSTELLAR_DISC_EXPONENTIAL
OPT += -DEXPONENTIAL_DISC
OPT += -DGalIC
OPT += -DHDF5
#OPT += -DDEBUG

OPTS = 

ifeq (HDF5,$(findstring HDF5,$(OPT))) 
FC = /opt/homebrew/bin/gfortran-15
HDF5_INCL = -I/opt/homebrew/Cellar/hdf5/2.1.1/include
HDF5_LIBS  = -L/opt/homebrew/Cellar/hdf5/2.1.1/lib \
	-lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
	-L/opt/homebrew/opt/libaec/lib -lsz -lz -ldl -lm
#	/opt/homebrew/Cellar/hdf5/2.1.1/lib/libhdf5_hl_fortran.a /opt/homebrew/Cellar/hdf5/2.1.1/lib/libhdf5_hl.a /opt/homebrew/Cellar/hdf5/2.1.1/lib/libhdf5_fortran.a /opt/homebrew/Cellar/hdf5/2.1.1/lib/libhdf5.a -L/opt/homebrew/opt/libaec/lib -lsz -lz -ldl -lm
else
FC = gfortran-9 $(OPTS)
HDF5_INCL= 
HDF5_LIBS=
endif

OPTS= $(OPT) 
#OPTS += -ffpe-trap=zero,overflow,underflow
OPTS += -O3 -march=native -funroll-loops -ffast-math

SRCS =  compute_background_potential.f90 mk_gas_disc.f90 merge_comps.f90 mk_gas_halo.f90 mk_nbody_components.f90
SOBJ = $(SRCS:.f90=.o)

FILE = constants.f90 gadget_header.f90 nrutils_modules.f90 structure.f90 read_params.f90 gas_halo.f90 io.f90 compute_parameters.f90 make_input.f90
FOBJ = $(FILE:.f90=.o)

all: mk_gas_halo.exe mk_gas_disc.exe mk_nbody_comps.exe merge_comps.exe

merge_comps.exe: $(FOBJ) merge_comps.o
	$(FC) $^ -o merge_comps.exe $(HDF5_LIBS) 

mk_gas_halo.exe: $(FOBJ) mk_gas_halo.o
	$(FC) $^ -o mk_gas_halo.exe $(HDF5_LIBS) 

mk_gas_disc.exe: $(FOBJ) compute_background_potential.o mk_gas_disc.o
	$(FC) $^ -o mk_gas_disc.exe $(HDF5_LIBS) 

mk_nbody_comps.exe: $(FOBJ) mk_nbody_comps.o
	$(FC) $^ -o mk_nbody_comps.exe $(HDF5_LIBS)

%.o: %.f90
	$(FC) -cpp $(OPTS) ${HDF5_INCL} -c $<
	touch $*.o $*.mod

clean: 
	rm mk_gas_disc.exe mk_gas_halo.exe mk_nbody_comps.exe *.mod *.o merge_comps.exe
	touch makefile
