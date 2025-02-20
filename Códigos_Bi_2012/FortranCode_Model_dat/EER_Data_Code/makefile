
PROG =	model
SRCS =  model.f90 mod_nonlinear_solver.f90 mod_fcn_eulor.f90 mod_parameter.f90 mod_interp_locate.f90 mod_global_constant.f90
OBJS =	model.o mod_nonlinear_solver.o mod_fcn_eulor.o mod_parameter.o mod_interp_locate.o mod_global_constant.o
LIBDIR = /opt/intel/Compiler/11.0/084/mkl/lib/em64t
LIBS =	-limf -lm -lmkl -lmkl_lapack -lguide -lpthread -limsl -limslsuperlu -limslscalar -limslblas -L/usr/local/vni/imsl/fnl600/lnxin100e64/lib -I/usr/local/vni/imsl/fnl600/lnxin100e64/include

FFLAGS = 
F90 = ifort
F90FLAGS = -O2 -shared-intel -threads -CB -mcmodel=medium -c ${LIBS}

all: $(PROG)

model: $(OBJS)
	$(F90) -O2 -shared-intel -threads -CB -mcmodel=medium -o $(PROG) $(OBJS) ${LIBS}

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

mod_interp_locate.o:	mod_global_constant.o
mod_parameter.o:		mod_interp_locate.o mod_global_constant.o
mod_fcn_eulor.o:		mod_parameter.o mod_global_constant.o mod_interp_locate.o
mod_nonlinear_solver.o:	mod_fcn_eulor.o mod_global_constant.o
model.o:				mod_nonlinear_solver.o mod_fcn_eulor.o mod_parameter.o mod_global_constant.o