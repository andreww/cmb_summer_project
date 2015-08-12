FC=gfortran -O0 -gdwarf-2
#
IASQ=emtvelin.o
TAUL=libtau.o libsun.o
#
#  path of executables
BIN = 
default: remodlv setbrn ttimes
#
#  resampling of model 
REM = remodl.o
remodlv: $(REM) $(TAUL) $(IASQ)
	$(FC) -o $(BIN)remodlv $(REM) $(IASQ) $(TAUL)
#
#  set branches  program
SBR = setbrn.o
setbrn: $(SBR) $(TAUL) 
	gfortran -o $(BIN)setbrn $(SBR) $(TAUL)
#
# travel time branches program
TBR =  ttimes.o       
ttimes: $(TBR) $(TAUL)
	gfortran -o $(BIN)ttimes $(TBR) $(TAUL)
#
libtau.o: libtau.f ttlim.inc
	$(FC) -c libtau.f
clean:
	rm $(TAUL) $(TBR) $(SBR) $(REM) $(IASQ)
	rm $(BIN)remodlv $(BIN)setbrn $(BIN)ttimes
