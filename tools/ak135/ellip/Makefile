# Override default fortran compiler for F77
FC=gfortran
#
TAUD=
TAUL=./libtau.o ./libsun.o
#
#  path of executables
BIN =  
default: ttimel direct

libtau.o: libtau.f ttlim.inc
	gfortran -c libtau.f
#
# ellipticity correction to travel times
TEF = ttimel.o ellip.o      
ttimel: $(TEF) $(TAUL) 
	gfortran -o $(BIN)ttimel $(TEF) $(TAUL)
#
# direct access conversion of ellipticity corrections
DEF = direct.o     
direct: $(DEF) 
	gfortran -o $(BIN)direct $(DEF)

clean:
	rm $(DEF) $(TEF) $(BIN)direct $(BIN)ttimel
