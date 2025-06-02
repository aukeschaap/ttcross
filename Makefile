# --------------------------------------------------
# File Suffixes
# --------------------------------------------------
.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o

# --------------------------------------------------
# Compiler Configuration
# --------------------------------------------------
FC      = mpif90 -ffree-line-length-none      # Fortran compiler
CC      = mpicc                               # C compiler
LDR     = mpif90                              # Linker
OPT     = -O2 -fopenmp -fno-strict-aliasing   # Optimization flags

# --------------------------------------------------
# External Libraries
# --------------------------------------------------
BLASLIB  = -llapack -lblas
MPFLIB   = -lmpfr
MPD      = ./mpfun-mpfr-v08              # MPFUN source directory

# HDF5 support
HDF5_INC = -I/usr/include/hdf5/serial
HDF5_LIB = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5

# --------------------------------------------------
# Directory Layout
# --------------------------------------------------
BUILD   = build
MODDIR  = $(BUILD)/mod
OBJDIR  = $(BUILD)/obj
BINDIR  = $(BUILD)/bin
LIBDIR  = lib

# --------------------------------------------------
# Source File Lists
# --------------------------------------------------

# Files with no dependencies between them
NO_DEPS = zero nan trans default timef say rnd ptype funcs constants s_vectors cos_approx mvn_pdf quad

# All source modules (in desired build order)
SRC = zero nan trans default timef say rnd ptype ort lr mat quad \
      tt ttind ttio dmrgg mvn_pdf cos_approx utils s_vectors funcs constants coefficients

# MPFUN sources (from mpfun-mpfr-v08/)
MPF = mpfuna mpfunf mpfung1 mpinterface mpmodule mpblas ttmp dmrggmp

# Test programs (in root dir)
TESTS = test_crs_ising \
		test_crs_stdnorm \
		test_crs_mvn \
		test_crs_mvn_complex \
        test_crs_chf \
		test_crs_pdf \
		test_crs_store \
		test_s_vectors \
		test_crs_coscoeff \
		test_print_cos_coeff \
		test_chf_equal

# Object file expansions
OBJ     = $(foreach s,$(SRC),$(OBJDIR)/$(s).o)
MPOBJ   = $(foreach s,$(MPF),$(OBJDIR)/$(s).o)
TESTOBJ = $(foreach t,$(TESTS),$(OBJDIR)/$(t).o)

# --------------------------------------------------
# Default Target
# --------------------------------------------------
all: | $(BUILD) $(MODDIR) $(OBJDIR) $(BINDIR)
all: ordered_compile $(foreach t,$(TESTS),$(BINDIR)/$(t).exe)

# --------------------------------------------------
# Directory Creation
# --------------------------------------------------
$(BUILD) $(MODDIR) $(OBJDIR) $(BINDIR):
	mkdir -p $@

# --------------------------------------------------
# Rule: No-dependency module compilation (lib/*.f90)
# --------------------------------------------------
$(addprefix $(OBJDIR)/,$(addsuffix .o,$(NO_DEPS))): $(OBJDIR)/%.o: $(LIBDIR)/%.f90
	$(FC) $(OPT) $(HDF5_INC) -J$(MODDIR) -I$(MODDIR) -c $< -o $@

# --------------------------------------------------
# Rule: Remaining module compilation (lib/*.f90)
# --------------------------------------------------
$(OBJDIR)/%.o: $(LIBDIR)/%.f90
	$(FC) $(OPT) $(HDF5_INC) -J$(MODDIR) -I$(MODDIR) -c $< -o $@

# --------------------------------------------------
# Rule: Test program compilation (top-level *.f90)
# --------------------------------------------------
$(OBJDIR)/%.o: %.f90
	$(FC) $(OPT) $(HDF5_INC) -I$(MODDIR) -c $< -o $@

# --------------------------------------------------
# Rule: Link test executables
# Links all test programs with all module objects
# --------------------------------------------------
$(BINDIR)/%.exe: $(OBJ) $(OBJDIR)/%.o
	$(LDR) $(OPT) $^ -o $@ $(BLASLIB) $(HDF5_LIB)

# --------------------------------------------------
# Rule: Compile MPFUN sources
# --------------------------------------------------
$(OBJDIR)/mpfuna.o: $(MPD)/mpfuna.f90
	$(FC) $(OPT) -fno-underscoring -J$(MODDIR) -I$(MODDIR) -c $< $(MPFLIB) -o $@

$(OBJDIR)/mpfunf.o: $(MPD)/mpfunf.f90
	$(FC) $(OPT) -fno-underscoring -J$(MODDIR) -I$(MODDIR) -c $< $(MPFLIB) -o $@

$(OBJDIR)/mpfung1.o: $(MPD)/mpfung1.f90
	$(FC) $(OPT) -fno-underscoring -J$(MODDIR) -I$(MODDIR) -c $< $(MPFLIB) -o $@

$(OBJDIR)/mpmodule.o: $(MPD)/mpmodule.f90
	$(FC) $(OPT) -fno-underscoring -J$(MODDIR) -I$(MODDIR) -c $< $(MPFLIB) -o $@

$(OBJDIR)/mpinterface.o: $(MPD)/mpinterface.c
	$(CC) $(OPT) -c $< -o $@

# --------------------------------------------------
# Explicit Dependencies Between Modules
# Ensures correct build order when modules `use` each other
# --------------------------------------------------
$(OBJDIR)/ort.o:           $(LIBDIR)/ort.f90           $(OBJDIR)/nan.o
$(OBJDIR)/mat.o:           $(LIBDIR)/mat.f90           $(OBJDIR)/nan.o
$(OBJDIR)/lr.o:            $(LIBDIR)/lr.f90            $(OBJDIR)/ort.o $(OBJDIR)/default.o
$(OBJDIR)/tt.o:            $(LIBDIR)/tt.f90            $(OBJDIR)/ptype.o $(OBJDIR)/default.o $(OBJDIR)/mat.o $(OBJDIR)/trans.o $(OBJDIR)/timef.o $(OBJDIR)/say.o
$(OBJDIR)/ttind.o:         $(LIBDIR)/ttind.f90         $(OBJDIR)/tt.o
$(OBJDIR)/ttio.o:          $(LIBDIR)/ttio.f90          $(OBJDIR)/tt.o
$(OBJDIR)/dmrgg.o:         $(LIBDIR)/dmrgg.f90         $(OBJDIR)/tt.o $(OBJDIR)/lr.o $(OBJDIR)/rnd.o $(OBJDIR)/default.o $(OBJDIR)/timef.o $(OBJDIR)/ptype.o
$(OBJDIR)/utils.o:         $(LIBDIR)/utils.f90         $(OBJDIR)/tt.o
$(OBJDIR)/coefficients.o:  $(LIBDIR)/coefficients.f90  $(OBJDIR)/s_vectors.o $(OBJDIR)/funcs.o $(OBJDIR)/constants.o

# --------------------------------------------------
# Ordered Compilation Target (optional)
# Forces all SRC files to be compiled in order
# --------------------------------------------------
ordered_compile: $(addprefix $(OBJDIR)/,$(addsuffix .o,$(SRC)))

# --------------------------------------------------
# Clean
# --------------------------------------------------
clean:
	rm -rf $(BUILD)
