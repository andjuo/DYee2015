CC=g++
RUDir=../RooUnfold/src
CFLAGS=-D__CXX__ -I$(RUDIR) -I./  -I $(shell root-config --incdir)
LFLAGS=-L./ -lRooUnfold $(shell root-config --libs)
OBJ=inputs.o DYbinning.o DYmm13TeV.o DYmm13TeV_eff.o DYee13TeV.o crossSection.o

default: $(OBJ)

clean:
	rm -f *.o

veryclean:
	rm -f *.o *.exe

processDYee_dressed.exe: processDYee_dressed.C $(OBJ)
	$(CC) $(CFLAGS) -o $@ $< $(OBJ) $(LFLAGS)


inputs.o: inputs.cc inputs.h
	$(CC) $(CFLAGS) -c -o $@ $<

DYbinning.o: DYbinning.cxx DYbinning.h
	$(CC) $(CFLAGS) -c -o $@ $<

DYmm13TeV.o: DYmm13TeV.cc DYmm13TeV.h DYbinning.o inputs.o
	$(CC) $(CFLAGS) -c -o $@ $<

DYmm13TeV_eff.o: DYmm13TeV_eff.cc DYmm13TeV_eff.h DYbinning.o inputs.o
	$(CC) $(CFLAGS) -c -o $@ $<

DYee13TeV.o: DYee13TeV.cc DYee13TeV.h DYbinning.o
	$(CC) $(CFLAGS) -c -o $@ $<

crossSection.o: crossSection.cc crossSection.h inputs.o \
	$(RUDir)/RooUnfoldResponse.h $(RUDir)/RooUnfoldBayes.h \
	$(RUDir)/RooUnfoldInvert.h
	$(CC) $(CFLAGS) -c -o $@ $<

