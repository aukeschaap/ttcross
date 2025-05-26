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

BUILD   = build
MODDIR  = $(BUILD)/mod
OBJDIR  = $(BUILD)/obj
BINDIR  = $(BUILD)/bin

SRC     = zero nan trans default timef say rnd ptype ort lr mat quad tt ttind ttio dmrgg qmc mc mvn_pdf
MPF     = mpfuna mpfunf mpfung1 mpinterface mpmodule mpblas ttmp dmrggmp
OBJ     = $(foreach s,$(SRC),$(OBJDIR)/$(s).o)
MPOBJ   = $(foreach s,$(MPF),$(OBJDIR)/$(s).o)


MAIN    = ising

all: | $(BUILD) $(MODDIR) $(OBJDIR) $(BINDIR)
all: $(BINDIR)/test_crs_ising.exe $(BINDIR)/test_qmc_ising.exe $(BINDIR)/test_mc_ising.exe $(BINDIR)/test_crs_stdnorm.exe $(BINDIR)/test_crs_mvn.exe $(BINDIR)/test_crs_mvn_complex.exe

$(BUILD) $(MODDIR) $(OBJDIR) $(BINDIR):
	mkdir -p $@

$(OBJDIR)/%.o: %.f90
	$(FC) $(OPT) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/%.o: %.f
	$(FC) $(OPT) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/%.o: %.F90
	$(FC) $(OPT) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/%.o: %.c
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

$(OBJDIR)/dmrgg.o: dmrgg.f90
	$(FC) $(OPT) -fallow-argument-mismatch -J$(MODDIR) -c $< -o $@

$(OBJDIR)/qmc.o: qmc.f90
	$(FC) $(OPT) -fallow-argument-mismatch -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mc.o: mc.f90
	$(FC) $(OPT) -fallow-argument-mismatch -J$(MODDIR) -c $< -o $@

$(BINDIR)/test_mpf_ising.exe: $(OBJ) $(MPOBJ) $(OBJDIR)/test_mpf_ising.o
	$(LDR) $(OPT) $^ -o $@ $(BLASLIB) $(MPFLIB)

$(BINDIR)/%.exe: $(OBJ) $(OBJDIR)/%.o
	$(LDR) $(OPT) $^ -o $@ $(BLASLIB)

clean:
	rm -rf $(BUILD)
