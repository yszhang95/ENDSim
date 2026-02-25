#include <stdio.h>
#include <iostream>
#include <vector>

#include <RAT/DSReader.hh>
#include <RAT/DS/MC.hh>
//#include <RAT/DS/RunStore.hh>
//#include <RAT/DS/Run.hh>
//#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCSummary.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/EV.hh>

#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>

#include <glob.h>
#include <algorithm>
#include <numeric>

// x in (-25, 30)
// y in (-22, 22)
// z in (-25, 105)
// transverse is y <-> z

using namespace std;

int ENDpmt2dom(int pmt_id){
    return (int)pmt_id/31;
}

void process(std::string pattern, std::string outfilename){

    glob_t glob_result;
    int result = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if(result != 0){
        std::cerr << "Couldn't find files! Exiting!" << std::endl;
        exit(1);
    }


    RAT::DSReader *dsreader = new RAT::DSReader(pattern.c_str());
    for(int i = 0; i < (int)glob_result.gl_pathc; i++)
        dsreader->Add(glob_result.gl_pathv[i]);
    const unsigned int nevents = dsreader->GetT()->GetEntries();
    std::cout << "NEvents : " << nevents << std::endl;
    TH1D* hmean = new TH1D("hmean", ";Mean PMT Hit Time; Events", 75, 0, 300);
    TH1D* hmax = new TH1D("hmax", ";Max #deltat b/w PMTs; Events", 75, 0, 300);
    TH1D* hstd = new TH1D("hstd", ";RMS spread of PMT Hit Times; Events", 75, 0, 300);

    // Loop over all triggered events
    for(size_t iev = 0; iev < nevents; iev++){

        RAT::DS::Root *rds = dsreader->GetEvent(iev);
        if(!rds->ExistMC()) continue;
        RAT::DS::MC *mc = rds->GetMC();
        for (int subev = 0; subev < rds->GetEVCount(); subev++) {
            double pmt_maxtime = 0.;
            double pmt_meantime = 0.;
            double pmt_stdtime = 0.;

            RAT::DS::EV *ev = rds->GetEV(subev);
            auto evid = ev->GetID();
            std::vector<double> hitPMTTime;
            for (int pmtc : ev->GetAllPMTIDs()) {
                RAT::DS::PMT *pmt = ev->GetOrCreatePMT(pmtc);
                hitPMTTime.push_back(pmt->GetTime());
            }
            if(hitPMTTime.size() == 0) continue;
            std::sort(hitPMTTime.begin(), hitPMTTime.end());
            pmt_meantime = std::accumulate(hitPMTTime.begin(),
                                           hitPMTTime.end(), 0.)/(hitPMTTime.size());
            hmean->Fill(pmt_meantime);
            if(hitPMTTime.size() > 1) {
                for(double &time: hitPMTTime){
                    pmt_stdtime += (time - pmt_meantime)*(time - pmt_meantime);
                }
                pmt_stdtime = std::sqrt(pmt_stdtime)/hitPMTTime.size();
                pmt_maxtime = hitPMTTime.at(hitPMTTime.size()-1) - hitPMTTime.at(0);

                hmax->Fill(pmt_maxtime);
                hstd->Fill(pmt_stdtime);
            }
        }
    }
    std::cout << "Writing to File!" << std::endl;
    TFile* outFile = new TFile(outfilename.c_str(), "recreate");
    hmax->Write();
    hmean->Write();
    hstd->Write();
    outFile->Close();
}

int main(int argc, char *argv[]){

    if(argc == 3){
        std::string input_pattern(argv[1]);
        std::string outfilename(argv[2]);
        process(input_pattern, outfilename);
    }
    else{
        std::cout << "Wrong number of arguments." << std::endl;
    }

    return 0;
}
