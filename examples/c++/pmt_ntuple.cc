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

#include <TTree.h>
#include <TFile.h>

#include <glob.h>

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
    int pmt_i, dom_i, npe_i, pmt_x, pmt_y, pmt_z;
    int ev_i, nch_total, npmt;
    std::vector<Int_t> pdgcodes;
    std::vector<int> doms;
    std::vector<int> pmts;
    std::vector<int> pes;
    std::vector<double> mcKEnergies;
    std::vector<double> mcPosx;
    std::vector<double> mcPosy;
    std::vector<double> mcPosz;
    std::vector<double> mcDirx;
    std::vector<double> mcDiry;
    std::vector<double> mcDirz;

    TTree* pmtTree = new TTree("pmt", "pmt");
    pmtTree->Branch("ev_i", &ev_i);
    pmtTree->Branch("pmt_i", &pmt_i);
    pmtTree->Branch("dom_i", &dom_i);
    pmtTree->Branch("npe_i", &npe_i);
    //    pmtTree->Branch("pmt_x", &pmt_x);
    //    pmtTree->Branch("pmt_y", &pmt_y);
    //    pmtTree->Branch("pmt_z", &pmt_z);


    TTree* evTree = new TTree("ev", "ev");
    evTree->Branch("npmt", &npmt);
    evTree->Branch("nch_total", &nch_total);
    evTree->Branch("ev_i", &ev_i);
    evTree->Branch("mcpdgs", &pdgcodes);
    evTree->Branch("pes", &pes);
    evTree->Branch("doms", &doms);
    evTree->Branch("pmts", &pmts);
    evTree->Branch("mcxs", &mcPosx);
    evTree->Branch("mcys", &mcPosy);
    evTree->Branch("mczs", &mcPosz);
    evTree->Branch("mcus", &mcDirx);
    evTree->Branch("mcvs", &mcDiry);
    evTree->Branch("mcws", &mcDirz);
    evTree->Branch("mckes", &mcKEnergies);

    // Loop over all triggered events
    for(size_t iev = 0; iev < nevents; iev++){

        RAT::DS::Root *rds = dsreader->GetEvent(iev);
        if(!rds->ExistMC()) continue;
        RAT::DS::MC *mc = rds->GetMC();
        //        RAT::DS::Run* runBranch = RAT::DS::RunStore::GetRun(rds);
        //        RAT::DS::PMTInfo *pmtinfo = runBranch->GetPMTInfo();

        // Use the MCSummary to get total # ch photons produced
        RAT::DS::MCSummary *mcs = mc->GetMCSummary();
        nch_total = mcs->GetNumCerenkovPhoton();
        ev_i = iev;
        npmt = mc->GetMCPMTCount();
        int mcpcount = mc->GetMCParticleCount();
        for (int pid = 0; pid < mcpcount; pid++) {
            RAT::DS::MCParticle *particle = mc->GetMCParticle(pid);
            pdgcodes.push_back(particle->GetPDGCode());
            mcKEnergies.push_back(particle->GetKE());
            TVector3 mcpos = particle->GetPosition();
            TVector3 mcdir = particle->GetMomentum();
            mcPosx.push_back(mcpos.X());
            mcPosy.push_back(mcpos.Y());
            mcPosz.push_back(mcpos.Z());
            mcDirx.push_back(mcdir.X() / mcdir.Mag());
            mcDiry.push_back(mcdir.Y() / mcdir.Mag());
            mcDirz.push_back(mcdir.Z() / mcdir.Mag());
        }
        // Loop over all hit PMTs
        for(int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++){

            RAT::DS::MCPMT* mcpmt = mc->GetMCPMT(ipmt);
            // Total npe detected by a PMT
            npe_i = mcpmt->GetMCPhotonCount();
            pmt_i = mcpmt->GetID();
            dom_i = ENDpmt2dom(pmt_i);
            pes.push_back(npe_i);
            pmts.push_back(pmt_i);
            doms.push_back(dom_i);
            //            TVector3 position = pmtinfo->GetPosition(pmt_i);
            //            pmt_x = position.X();
            //            pmt_y = position.Y();
            //            pmt_z = position.Z();

            // For End, type=0 is 8'', type=1 is 12''
            int type = mcpmt->GetType();
            if(type == 0) continue; // Skip the HQE PMTs in this example
            pmtTree->Fill();

        }
        evTree->Fill();

        pes.clear();
        pmts.clear();
        doms.clear();
        pdgcodes.clear();
        mcKEnergies.clear();
        mcPosx.clear();
        mcPosy.clear();
        mcPosz.clear();
        mcDirx.clear();
        mcDiry.clear();
        mcDirz.clear();
    }

    std::cout << "Writing to File!" << std::endl;
    TFile* outFile = new TFile(outfilename.c_str(), "recreate");
    pmtTree->Write();
    evTree->Write();
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
