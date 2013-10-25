COMP=g++
ROOTLIBS = `root-config --glibs --cflags` -lMinuit -lRooFit -lRooStats -lRooFitCore -lEG
INC= -I.. -I. -I./include  -I${CLHEP}/include
ROOTINC= -I${ROOTSYS}/include -I${ROOFITSYS}/include
LIBS= -L.  ${ROOTLIBS} -L${CLHEP}/lib -L${CLHEP}/lib -L${ROOFITSYS}/lib
SELECTIONLIB = SusyEventAnalyzer.o
EXE = ../macro/Analyze
EXE3 = UEDPlots

# ********** TEMPLATE *************
# mainProg: mainProg.o $(SELECTIONLIB)
#       $(COMP) $(INC) $(ROOTINC) $(LIBS) $(ROOTLIBS) -o $@  $(SELECTIONLIB) $@.o
# *********************************

all: ${EXE}

printEventId: printEventId.o tree.o
	$(COMP) $(INC) $(ROOTINC) $(LIBS) -o $@ tree.o $@.o


clean:
	rm -f *.o *.lo core core.*
	rm -f *~
	rm -f *.exe
	rm -f $(EXE)
	rm -f ../macro/Analyze_Filelist
	rm -f SusyEventDict.*
.cpp.o:
	$(COMP) -m64 -c $(INC) $(ROOTINC) -o $@ $<

.cc.o:
	$(COMP) -m64 -c $(INC) $(ROOTINC) -o $@ $<

.cxx.o:
	$(COMP) -m64 -c $(INC) $(ROOTINC) -o $@ $<

.C.o:
	$(COMP) -m64 -c $(INC) $(ROOTINC) -o $@ $<


UEDPlots: UEDPlots.o roostats_cl95.o
	$(COMP) $(INC) $(ROOTINC) $(LIBS) -o $@ roostats_cl95.o $@.o

makeClosure: makeClosure.o

Analyze: Analyze.o SusyEvent.o SusyEventAnalyzer.o SusyEventPrinter.o SusyEventDict.o
	$(COMP) $(INC) $(ROOTINC) $(LIBS) -o $@ SusyEvent.o SusyEventAnalyzer.o SusyEventPrinter.o ../jec/tmp/*.o SusyEventDict.o $@.o

Analyze_Filelist: Analyze_Filelist.o SusyEvent.o SusyEventAnalyzer.o SusyEventDict.o
	$(COMP) $(INC) $(ROOTINC) $(LIBS) -o $@ SusyEvent.o SusyEventAnalyzer.o ../jec/tmp/*.o SusyEventDict.o $@.o


SusyEvent.o: SusyEvent.h

SusyEventPrinter.o: SusyEventPrinter.h

SusyEventDict.cc: SusyEvent.h SusyNtuplizer_LinkDef.h
	@echo "Generating dictionary ..."
	@rootcint SusyEventDict.cc -c SusyEvent.h SusyNtuplizer_LinkDef.h
