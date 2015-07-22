#
TAUD=
TAUL=./libtau.o ./libsun.o
#
#  path of executables
BIN =  
#
# ellipticity correction to travel times
TEF = ttimel.o ellip.o      
ttimel: $(TEF) $(TAUL) 
#	f77 -o $(BIN)ttimel $(TEF) $(TAUL)
	ifort -o $(BIN)ttimel $(TEF) $(TAUL)
#
# direct access conversion of ellipticity corrections
DEF = direct.o     
direct: $(DEF) 
#	f77 -o $(BIN)direct $(DEF))
	ifort -o $(BIN)direct $(DEF)
