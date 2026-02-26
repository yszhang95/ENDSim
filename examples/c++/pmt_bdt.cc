#include <stdio.h>
#include <iostream>
#include <vector>

#include <RAT/DSReader.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/PMTInfo.hh>
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

using namespace std;

int ENDpmt2dom(int pmt_id){
    return (int)pmt_id/31;
}

int GetDomX(int dom_id) {
    int dom_plane = dom_id/18;
    if (dom_plane == 0){ return 1; }
    if (dom_plane == 1){ return 0; }
    if (dom_plane == 2){ return 1; }
    if (dom_plane == 3){ return 2; }
    if (dom_plane == 4){ return 2; }
    return 0;
}

int GetDomY(int dom_id) {
    int dom_plane = dom_id/18;
    if (dom_plane == 0){ return 2; }
    if (dom_plane == 1){ return 1; }
    if (dom_plane == 2){ return 0; }
    if (dom_plane == 3){ return 0; }
    if (dom_plane == 4){ return 1; }
    return 0;
}

int GetDomZ(int dom_id) {
    return (int)dom_id % 18;
}

void process(std::string pattern, std::string outfilename){

    TFile* xsec = new TFile("xsec_graphs.root", "read");
    TGraph* gr_Si28 = (TGraph*)xsec->Get("nu_mu_Si28/tot_cc");
    TGraph* gr_O16 = (TGraph*)xsec->Get("nu_mu_O16/tot_cc");
    TGraph* gr_H1 = (TGraph*)xsec->Get("nu_mu_H1/tot_cc");
    xsec->Close();
    double xsec_weight = (gr_Si28->Eval(3.5) + 2.*gr_O16->Eval(3.5))/(gr_O16->Eval(3.5) + 2.*gr_H1->Eval(3.5));

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

    // get pmt positions
    TTree* runT = dsreader->GetRunT();
    auto run = new RAT::DS::Run();
    runT->SetBranchStatus("*", 1);
    runT->SetBranchAddress("run", &run);
    runT->GetEntry(0);
    auto pmtinfo = run->GetPMTInfo();
    int NPMTS = 2790;
    std::vector<double> dom_xs, dom_ys, dom_zs;
    double mean_x = 0, mean_y = 0, mean_z = 0;
    for(int i = 0; i < NPMTS; i++){
        TVector3 pmtpos = pmtinfo->GetPosition(i);
        mean_x += pmtpos[0]/31;
        mean_y += pmtpos[1]/31;
        mean_z += pmtpos[2]/31;
        if((i % 31) == 30) {
            dom_xs.push_back(mean_x);
            dom_ys.push_back(mean_y);
            dom_zs.push_back(mean_z);
            mean_x = 0.;
            mean_y = 0.;
            mean_z = 0.;
        }
    }

    // ttree at the dom level
    TH1D* hevts = new TH1D("nevts", "Total Evts", 1, 0, 1);
    TTree* tout = new TTree("doms", "doms");

    float dom_x, dom_y, dom_z;
    int npmts, npe, pe_min, pe_spread, pe_max;
    int npmts_50, npmts_100, npe_50, npe_100;
    float pe_mean, pe_rms;
    float t_min, t_spread, t_mean, t_rms;
    int ndoms;
    int event_id;
    int du_id, dom_id;
    float vtxX, vtxY, vtxZ;
    float momX, momY, momZ;
    float muE;
    float rock_wgt;

    tout->Branch("event_id", &event_id, "event_id/I");
    tout->Branch("dom_id", &dom_id, "dom_id/I");
    tout->Branch("du_id", &du_id, "du_id/I");
    tout->Branch("vtxX", &vtxX, "vtxX/F");
    tout->Branch("vtxY", &vtxY, "vtxY/F");
    tout->Branch("vtxZ", &vtxZ, "vtxZ/F");
    tout->Branch("momX", &momX, "momX/F");
    tout->Branch("momY", &momY, "momY/F");
    tout->Branch("momZ", &momZ, "momZ/F");
    tout->Branch("muE", &muE, "muE/F");
    tout->Branch("rock_wgt", &rock_wgt, "rock_wgt/F");
    tout->Branch("dom_x", &dom_x, "dom_x/F");
    tout->Branch("dom_y", &dom_y, "dom_y/F");
    tout->Branch("dom_z", &dom_z, "dom_z/F");
    tout->Branch("npmts", &npmts, "npmts/I");
    tout->Branch("npe", &npe, "npe/I");
    tout->Branch("pe_min", &pe_min, "pe_min/I");
    tout->Branch("pe_spread", &pe_spread, "pe_spread/I");
    tout->Branch("pe_rms", &pe_rms, "pe_rms/F");
    tout->Branch("t_min", &t_min, "t_min/F");
    tout->Branch("t_spread", &t_spread, "t_spread/F");
    tout->Branch("t_mean", &t_mean, "t_mean/F");
    tout->Branch("t_rms", &t_rms, "t_rms/F");
    tout->Branch("npmts_50", &npmts_50, "npmts_50/I");
    tout->Branch("npe_50", &npe_50, "npe_50/I");
    tout->Branch("npmts_100", &npmts_100, "npmts_100/I");
    tout->Branch("npe_100", &npe_100, "npe_100/I");

    // Loop over all triggered events
    double nwgt_events = 0;
    for(size_t iev = 0; iev < nevents; iev++){

        RAT::DS::Root *rds = dsreader->GetEvent(iev);
        if(!rds->ExistMC()) continue;
        if(rds->GetEVCount() == 0) continue;

        RAT::DS::MC *mc = rds->GetMC();
        RAT::DS::MCSummary *mcs = mc->GetMCSummary();
        RAT::DS::EV *ev = rds->GetEV(0);

        int mcpcount = mc->GetMCParticleCount();
        vtxX = -100000.;
        vtxY = -100000.;
        vtxZ = -100000.;
        momX = -100000.;
        momY = -100000.;
        momZ = -100000.;
        muE = -5;
        for (int pid = 0; pid < mcpcount; pid++) {
            RAT::DS::MCParticle *particle = mc->GetMCParticle(pid);

            if (abs(particle->GetPDGCode()) == 13) {
                TVector3 mcpos = particle->GetPosition();
                vtxX = mcpos.X();
                vtxY = mcpos.Y();
                vtxZ = mcpos.Z();
                TVector3 mcmom = particle->GetMomentum();
                momX = mcmom.X();
                momY = mcmom.Y();
                momZ = mcmom.Z();
                muE = particle->GetKE()/1000.;
            }
        }
        rock_wgt = 1.;
        if(((pattern.find("rockbed") != string::npos) || (pattern.find("equalx") != string::npos) || (pattern.find("aframe") != string::npos)) && (pattern.find("cosmic") == string::npos)){
            // up-weight rock events manually
            if(vtxY > -23000. && vtxY < -13000. && (pattern.find("short") == string::npos)){
                rock_wgt = xsec_weight;
            }
            // if(vtxY > -13000. && vtxY < -7000. && (pattern.find("short") == string::npos) && ((pattern.find("hex") == string::npos) && (pattern.find("box") == string::npos))){
            //     rock_wgt = 0.;
            // }
            if(vtxY > -22500. && vtxY < -12500. && (pattern.find("short") != string::npos)){
                rock_wgt = xsec_weight;
            }
            if(vtxY < -9000 && (pattern.find("DU") != string::npos)){
                rock_wgt = 0.;
            }
            // if(vtxY > -17000. && vtxY < -7000. && ((pattern.find("hex") != string::npos) || (pattern.find("box") != string::npos))){
            //     rock_wgt = xsec_weight;
            // }
        }
        nwgt_events += rock_wgt;

        event_id = (int)iev;
        std::map<int, int> dom_pmts;
        std::map<int, std::vector<int>> dom_pes;
        std::map<int, std::vector<float>> dom_times;

        for(int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++){

            RAT::DS::MCPMT* mcpmt = mc->GetMCPMT(ipmt);
            // Total npe detected by a PMT
            int npe_i = mcpmt->GetMCPhotonCount();
            int pmt_i = mcpmt->GetID();
            int dom_i = ENDpmt2dom(pmt_i);
            // For End, type=0 is 8'', type=1 is 12''
            int type = mcpmt->GetType();
            if(type == 0) continue; // Skip the HQE PMTs in this example
            RAT::DS::PMT *pmt = ev->GetOrCreatePMT(pmt_i);
            if(npe_i >= 1) {
                dom_pes[dom_i].push_back(npe_i);
                dom_pmts[dom_i] += 1;
                dom_times[dom_i].push_back((float)pmt->GetTime());
            }
        }


        for(auto &dom: dom_pmts){
            // dom_x = GetDomX(dom.first);
            // dom_y = GetDomY(dom.first);
            // dom_z = GetDomZ(dom.first);
            dom_id = dom.first;
            du_id = dom.first/18;
            dom_x = dom_xs[dom.first];
            dom_y = dom_ys[dom.first];
            dom_z = dom_zs[dom.first];
            npmts = dom.second;
            npmts_50 = 0;
            npmts_100 = 0;
            npe_50 = 0;
            npe_100 = 0;
            std::vector<int> dom_pe = dom_pes[dom.first];
            std::vector<float> dom_time = dom_times[dom.first];

            std::sort(dom_pe.begin(), dom_pe.end());
            std::sort(dom_time.begin(), dom_time.end());
            float init_time = dom_time.at(0);
            for(int i = 0; i < (int)dom_time.size(); i++){
                float pmt_time = dom_time.at(i);
                int pmt_pe = dom_pe.at(i);
                if(pmt_time - init_time <= 50){
                    npmts_50 += 1;
                    npe_50 += pmt_pe;
                    npmts_100 += 1;
                    npe_100 += pmt_pe;
                }
                else if(pmt_time - init_time <= 100){
                    npmts_100 += 1;
                    npe_100 += pmt_pe;
                    npmts_50 = 0;
                    npe_50 = pmt_pe;
                }
                else {
                    npmts_50 = 0;
                    npe_50 = pmt_pe;
                    npmts_100 = 0;
                    npe_100 = pmt_pe;
                }
                init_time = pmt_time;
            }

            npe = std::accumulate(dom_pe.begin(), dom_pe.end(), 0);
            t_mean = std::accumulate(dom_time.begin(), dom_time.end(), 0)/dom_time.size();
            pe_mean = npe/dom_pe.size();
            pe_min = dom_pe.at(0);
            pe_max = dom_pe.at(dom_pe.size()-1);
            pe_spread = pe_max - pe_min;
            t_spread = dom_time.at(dom_time.size()-1) - dom_time.at(0);
            t_min = dom_time.at(0);
            pe_rms = 0.;
            t_rms = 0.;
            for(int i = 0; i < (int)dom_pe.size(); i++){
                pe_rms += (TMath::Power(dom_pe[i]-pe_mean, 2));
                t_rms += (TMath::Power(dom_time[i]-t_mean, 2));
            }
            pe_rms = TMath::Sqrt(pe_rms)/dom_pe.size();
            t_rms = TMath::Sqrt(t_rms)/dom_time.size();

            tout->Fill();

        }
    }
    hevts->SetBinContent(1, nwgt_events);
    std::cout << "Writing to File!" << std::endl;
    TFile* outFile = new TFile(outfilename.c_str(), "recreate");
    tout->Write();
    hevts->Write();
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
