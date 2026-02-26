CPPFLAGS=-g -I$(ROOTSYS)/include -std=c++17 -I $(RATROOT)/include -I $(RATROOT)/include/stlplus
LDFLAGS=-g -L$(ROOTSYS)/lib -L $(RATROOT)/lib
LDLIBS=-lRATEvent `root-config --libs --cflags` -lcurl -pthread

ROOTCPPFLAGS=-g -I$(ROOTSYS)/include -std=c++17
ROOTLDFLAGS=-g -L$(ROOTSYS)/lib
ROOTLDLIBS=`root-config --libs --cflags`

all: pmt_ntuple pmt_ana pmt_event pmt_bdt ambient_bdt

timing_ntuple: timing_ntuple.cc
	$(CXX) $(LDFLAGS) -o timing_ntuple timing_ntuple.cc $(CPPFLAGS) $(LDLIBS)

timing_root: timing_root.cc
	$(CXX) $(LDFLAGS) -o timing_root timing_root.cc $(CPPFLAGS) $(LDLIBS)

pmt_ntuple: pmt_ntuple.cc
	$(CXX) $(LDFLAGS) -o pmt_ntuple pmt_ntuple.cc $(CPPFLAGS) -I"/ratpac-setup/local/include/Geant4" $(LDLIBS)
pmt_ana: pmt_ana.cc
	$(CXX) $(LDFLAGS) -o pmt_ana pmt_ana.cc $(CPPFLAGS) -I"/ratpac-setup/local/include/Geant4" $(LDLIBS)
pmt_event: pmt_event.cc
	$(CXX) $(LDFLAGS) -o pmt_event pmt_event.cc $(CPPFLAGS) -I"/ratpac-setup/local/include/Geant4" $(LDLIBS)
pmt_bdt: pmt_bdt.cc
	$(CXX) $(LDFLAGS) -o pmt_bdt pmt_bdt.cc $(CPPFLAGS) -I"/ratpac-setup/local/include/Geant4" $(LDLIBS)

ambient_bdt: ambient_bdt.cc
	$(CXX) $(ROOTLDFLAGS) -o ambient_bdt ambient_bdt.cc $(ROOTCPPFLAGS) $(ROOTLDLIBS)

# pmt_bdt2: pmt_bdt2.cc
# 	$(CXX) $(LDFLAGS) -o pmt_bdt2 pmt_bdt2.cc $(CPPFLAGS) -I"/ratpac-setup/local/include/Geant4" $(LDLIBS)

timing_root_digitizer: timing_root_digitizer.cc
	$(CXX) $(LDFLAGS) -o timing_root_digitizer timing_root_digitizer.cc $(CPPFLAGS) $(LDLIBS)

timing_ntuple.o: timing_ntuple.cc

timing_root.o: timing_root.cc

pmt_ntuple.o: pmt_ntuple.cc

timing_root.o: timing_root_digitizer.cc

clean:
	$(RM) pmt_ntuple timing_root timing_ntuple timing_root_digitizer
