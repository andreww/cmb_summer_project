FC = gfortran
FFLAGS = -Wall -fcheck=all -ffpe-trap=overflow,invalid,zero
LDFLAGS =
F2PY = f2py-2.7

LIBS = tomo_predict.so
MODS = tomo_predict.o
PROGS = tomo_predict_vdh

default: $(LIBS)
libs: $(LIBS)
progs: $(PROGS)
all: libs progs

$(PROGS): $(MODS)

%: $(MODS) %.o
	$(FC) $(LDFLAGS) -o $@ $(MODS) $*.o

tomo_predict_vdh.o: tomo_predict_vdh.f90

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $*.f90

%.so: %.o %.pyf
	$(F2PY) -c $*.pyf $*.f90

%.pyf: %.f90
	$(F2PY) $*.f90 -m $* -h $*.pyf

.PHONY: clean

clean:
	/bin/rm -f *.mod *.o *.pyf *.so tomo_predict_vdh
	/bin/rm -rf *.dSYM
