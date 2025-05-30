.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o

# set your Fortran and C compilers below
FC      = mpif90 -ffree-line-length-none
CC      = mpicc
LDR     = mpif90
OPT     = -O2 -fopenmp

BLASLIB = -llapack -lblas
MPFLIB  = -lmpfr
MPD     = ./mpfun-mpfr-v08

# HDF5 support
HDF5_INC = -I/usr/include/hdf5/serial
HDF5_LIB = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5


BUILD   = build
MODDIR  = $(BUILD)/mod
OBJDIR  = $(BUILD)/obj
BINDIR  = $(BUILD)/bin
LIBDIR  = lib

SRC     = zero nan trans default timef say rnd ptype ort lr mat quad tt ttind ttio dmrgg mvn_pdf cos_approx utils
MPF     = mpfuna mpfunf mpfung1 mpinterface mpmodule mpblas ttmp dmrggmp
OBJ     = $(foreach s,$(SRC),$(OBJDIR)/$(s).o)
MPOBJ   = $(foreach s,$(MPF),$(OBJDIR)/$(s).o)


MAIN    = ising

all: | $(BUILD) $(MODDIR) $(OBJDIR) $(BINDIR)
all: \
	$(BINDIR)/test_crs_ising.exe \
	$(BINDIR)/test_crs_stdnorm.exe  \
	$(BINDIR)/test_crs_mvn.exe  \
	$(BINDIR)/test_crs_mvn_complex.exe \
	$(BINDIR)/test_crs_chf.exe \
	$(BINDIR)/test_crs_pdf.exe \
	$(BINDIR)/test_crs_store.exe

$(BUILD) $(MODDIR) $(OBJDIR) $(BINDIR):
	mkdir -p $@

$(OBJDIR)/%.o: %.f90
	$(FC) $(OPT) $(HDF5_INC) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/%.o: $(LIBDIR)/%.f90
	$(FC) $(OPT) $(HDF5_INC) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/%.o: $(LIBDIR)/%.f
	$(FC) $(OPT) $(HDF5_INC) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/%.o: $(LIBDIR)/%.F90
	$(FC) $(OPT) $(HDF5_INC) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/%.o: $(LIBDIR)/%.c
	$(CC) $(OPT) -c $< -o $@

# MPFUN object rules
$(OBJDIR)/mpfuna.o: $(MPD)/mpfuna.f90
	$(FC) $(OPT) -fno-underscoring -J$(MODDIR) -c $< $(MPFLIB) -o $@

$(OBJDIR)/mpfunf.o: $(MPD)/mpfunf.f90
	$(FC) $(OPT) -fno-underscoring -J$(MODDIR) -c $< $(MPFLIB) -o $@

$(OBJDIR)/mpfung1.o: $(MPD)/mpfung1.f90
	$(FC) $(OPT) -fno-underscoring -J$(MODDIR) -c $< $(MPFLIB) -o $@

$(OBJDIR)/mpmodule.o: $(MPD)/mpmodule.f90
	$(FC) $(OPT) -fno-underscoring -J$(MODDIR) -c $< $(MPFLIB) -o $@

$(OBJDIR)/mpinterface.o: $(MPD)/mpinterface.c
	$(FC) $(OPT) -c $< $(MPFLIB) -o $@

$(OBJDIR)/dmrgg.o: $(LIBDIR)/dmrgg.f90
	$(FC) $(OPT) -fallow-argument-mismatch -J$(MODDIR) -c $< -o $@

$(BINDIR)/%.exe: $(OBJ) $(OBJDIR)/%.o
	$(LDR) $(OPT) $^ -o $@ $(BLASLIB) $(HDF5_LIB)

clean:
	rm -rf $(BUILD)
