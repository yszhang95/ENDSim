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
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>

#include <glob.h>
#include <algorithm>

// x in (-25, 30)
// y in (-22, 22)
// z in (-25, 105)
// transverse is y <-> z

using namespace std;

int ENDpmt2dom(int pmt_id){
    return (int)pmt_id/31;
}

void process(std::string pattern, std::string outfilename){

    TFile* xsec = new TFile("xsec_graphs.root", "read");
    TGraph* gr_Si28 = (TGraph*)xsec->Get("nu_mu_Si28/tot_cc");
    TGraph* gr_O16 = (TGraph*)xsec->Get("nu_mu_O16/tot_cc");
    TGraph* gr_H1 = (TGraph*)xsec->Get("nu_mu_H1/tot_cc");
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

    TH2D* hdom = new TH2D("hdom", "", 5, 0, 5, 18, 0, 18);
    std::vector<double> bound_x = {-25000., 30000.};
    std::vector<double> bound_y = {-23000., 23000.};
    std::vector<double> bound_z = {-25000., 105000.};
    std::vector<double> lenbound_x = {0., 60000.};
    std::vector<double> lenbound_y = {0., 23000.};
    std::vector<double> lenbound_z = {0., 130000.};
    std::vector<double> detbound_x = {0., 6000.};
    std::vector<double> detbound_y = {-500., 500.};
    std::vector<double> detbound_z = {0, 76500.};

    if(pattern.find("cosmic") != std::string::npos){
        bound_x = (std::vector<double>){-35000., 35000.};
        bound_y = (std::vector<double>){0, 5000.};
        bound_z = (std::vector<double>){-35000., 105000.};
    }

    if (pattern.find("trans") != std::string::npos){
        bound_y = (std::vector<double>){-25000., 105000.};
        bound_z = (std::vector<double>){-22000., 22000.};
        detbound_x = (std::vector<double>){0., 6000.};
        detbound_z = (std::vector<double>){-500., 500.};
        detbound_y = (std::vector<double>){0, 76500.};
    }
    // if (pattern.find("large") != std::string::npos){
    //     bound_z = (std::vector<double>){-25000., 200000.};
    //     detbound_x = (std::vector<double>){0., 18000.};
    //     detbound_z = (std::vector<double>){0, 153000.};
    // }
    if (pattern.find("3d") != std::string::npos){
        detbound_x = (std::vector<double>){-3000., 3000.};
        detbound_y = (std::vector<double>){-3000., 3000.};
        detbound_z = (std::vector<double>){0., 40500.};
    }
    if (pattern.find("pen_small") != std::string::npos){
        detbound_x = (std::vector<double>){-2000., 2000.};
        detbound_y = (std::vector<double>){-2000., 2000.};
        detbound_z = (std::vector<double>){0., 76500.};
    }
    if (pattern.find("pen_middle") != std::string::npos){
        detbound_x = (std::vector<double>){-5000., 5000.};
        detbound_y = (std::vector<double>){-5000., 5000.};
        detbound_z = (std::vector<double>){0., 76500.};
    }
    if (pattern.find("pen_large") != std::string::npos){
        detbound_x = (std::vector<double>){-8000., 8000.};
        detbound_y = (std::vector<double>){-8000., 8000.};
        detbound_z = (std::vector<double>){0., 76500.};
    }
    if (pattern.find("pen_tiny") != std::string::npos){
        detbound_x = (std::vector<double>){-1000., 1000.};
        detbound_y = (std::vector<double>){-1000., 1000.};
        detbound_z = (std::vector<double>){0., 76500.};
    }
    double x_bins[26], y_bins[26], z_bins[26];
    x_bins[0] = bound_x[0];
    y_bins[0] = bound_y[0];
    z_bins[0] = bound_z[0];
    for(int i = 1; i <= 25; i++){
        if(i <= 12){
            x_bins[i] = bound_x[0] + i*std::abs(detbound_x[0] - bound_x[0])/12.;
            y_bins[i] = bound_y[0] + i*std::abs(detbound_y[0] - bound_y[0])/12.;
            z_bins[i] = bound_z[0] + i*std::abs(detbound_z[0] - bound_z[0])/12.;
        }
        else{
            x_bins[i] = detbound_x[1] + (i-13)*std::abs(detbound_x[1] - bound_x[1])/12.;
            y_bins[i] = detbound_y[1] + (i-13)*std::abs(detbound_y[1] - bound_y[1])/12.;
            z_bins[i] = detbound_z[1] + (i-13)*std::abs(detbound_z[1] - bound_z[1])/12.;
        }
    }

    TH2D* hxy = new TH2D("hxy", "", 100, bound_x[0], bound_x[1], 100, bound_y[0], bound_y[1]);
    TH2D* hyz = new TH2D("hyz", "", 100, bound_z[0], bound_z[1], 100, bound_y[0], bound_y[1]);
    TH2D* hlenxy = new TH2D("hlenxy", "", 100, lenbound_x[0], lenbound_x[1], 100, lenbound_y[0], lenbound_y[1]);
    TH2D* hxz = new TH2D("hxz", "", 100, bound_z[0], bound_z[1], 100, bound_x[0], bound_x[1]);
    TH2D* hTotlenxy = new TH2D("hTotlenxy", "", 100, lenbound_x[0], lenbound_x[1], 100, lenbound_y[0], lenbound_y[1]);
    TH2D* hTotxz = new TH2D("hTotxz", "", 100, bound_z[0], bound_z[1], 100, bound_x[0], bound_x[1]);
    TH2D* he_xy3 = new TH2D("he_xy_dom3", "", 100, bound_x[0], bound_x[1], 100, bound_y[0], bound_y[1]);
    TH2D* he_xy4 = new TH2D("he_xy_dom4", "", 100, bound_x[0], bound_x[1], 100, bound_y[0], bound_y[1]);
    TH2D* he_xy5 = new TH2D("he_xy_dom5", "", 100, bound_x[0], bound_x[1], 100, bound_y[0], bound_y[1]);
    TH2D* he_xy6 = new TH2D("he_xy_dom6", "", 100, bound_x[0], bound_x[1], 100, bound_y[0], bound_y[1]);
    TH2D* he_xy7 = new TH2D("he_xy_dom7", "", 100, bound_x[0], bound_x[1], 100, bound_y[0], bound_y[1]);
    TH2D* hTotxy = new TH2D("hTotxy", "", 100, bound_x[0], bound_x[1], 100, bound_y[0], bound_y[1]);
    TH2D* hTotyz = new TH2D("hTotyz", "", 100, bound_z[0], bound_z[1], 100, bound_y[0], bound_y[1]);
    TH1D* hnDoms = new TH1D("hnDoms", "", 20, 0, 20);
    std::vector<TH3D*> hdoms_eff;
    // TH3D* hdoms_eff_denom = new TH3D("hdom_tot", "", 25, x_bins, 25, y_bins, 25, z_bins);
    TH3D* hdoms_eff_denom = new TH3D("hdom_tot", "", 100, bound_x[0], bound_x[1], 100, bound_y[0], bound_y[1], 100, bound_z[0], bound_z[1]);
    for(int i = 0; i < 15; i++){
        // auto hnum = new TH3D(TString::Format("hdom_thresh_%d", i+1), "", 25, x_bins, 25, y_bins, 25, z_bins);
        auto hnum = new TH3D(TString::Format("hdom_thresh_%d", i+1), "", 100, bound_x[0], bound_x[1], 100, bound_y[0], bound_y[1], 100, bound_z[0], bound_z[1]);
        hdoms_eff.push_back(hnum);
    }
    TH2D* hmuE = new TH2D("muE", "muE:nDOMs", 100, 0, 20, 15, 1, 16);

    std::cout << "NEvents : " << nevents << std::endl;

    // Loop over all triggered events
    for(size_t iev = 0; iev < nevents; iev++){

        RAT::DS::Root *rds = dsreader->GetEvent(iev);
        if(!rds->ExistMC()) continue;
        RAT::DS::MC *mc = rds->GetMC();
        //        RAT::DS::Run* runBranch = RAT::DS::RunStore::GetRun(rds);
        //        RAT::DS::PMTInfo *pmtinfo = runBranch->GetPMTInfo();

        // Use the MCSummary to get total # ch photons produced
        RAT::DS::MCSummary *mcs = mc->GetMCSummary();
        int mcpcount = mc->GetMCParticleCount();
        double vtxX = -100000.;
        double vtxY = -100000.;
        double vtxZ = -100000.;
        double lenX = -100000.;
        double lenY = -100000.;
        double lenZ = -100000.;
        double muE = -5.;
        for (int pid = 0; pid < mcpcount; pid++) {
            RAT::DS::MCParticle *particle = mc->GetMCParticle(pid);

            if (abs(particle->GetPDGCode()) == 13) {
                TVector3 mcpos = particle->GetPosition();
                vtxX = mcpos.X();
                vtxY = mcpos.Y();
                vtxZ = mcpos.Z();
                TVector3 mcendpos = particle->GetEndPosition();
                lenX = std::abs(mcendpos.X()-vtxX);
                lenY = std::abs(mcendpos.Y()-vtxY);
                lenZ = std::abs(mcendpos.Z()-vtxZ);
                muE = particle->GetKE() + 105.7;
            }
        }
        double rock_wgt = 1.;
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

        hTotxy->Fill(vtxX, vtxY, rock_wgt);
        hTotyz->Fill(vtxZ, vtxY, rock_wgt);
        hTotlenxy->Fill(lenX, lenY, rock_wgt);
        hTotxz->Fill(vtxZ, vtxX, rock_wgt);
        hdoms_eff_denom->Fill(vtxX, vtxY, vtxZ, rock_wgt);
        // Loop over all hit PMTs
        // std::vector<int> dom_ids;
        std::map<int, int> dom_pes;
        std::set<int> unique_doms;
        // int tot_npe = 0;
        // for(int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++){
        //
        //     RAT::DS::MCPMT* mcpmt = mc->GetMCPMT(ipmt);
        //     // Total npe detected by a PMT
        //     int npe_i = mcpmt->GetMCPhotonCount();
        //     int type = mcpmt->GetType();
        //     if(type == 0) continue; // Skip the HQE PMTs in this example
        //     tot_npe += npe_i;
        // }
        // if(tot_npe >= 5) {
            for(int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++){

                RAT::DS::MCPMT* mcpmt = mc->GetMCPMT(ipmt);
                // Total npe detected by a PMT
                int npe_i = mcpmt->GetMCPhotonCount();
                int pmt_i = mcpmt->GetID();
                int dom_i = ENDpmt2dom(pmt_i);
                // For End, type=0 is 8'', type=1 is 12''
                int type = mcpmt->GetType();
                if(type == 0) {
                    std::cout << "Found 8in" << "\n";
                    continue; // Skip the HQE PMTs in this example
                }
                if(npe_i >= 1) {
                    dom_pes[dom_i] += npe_i;
                    if(dom_pes[dom_i] >= 3)
                        unique_doms.insert(dom_i);
                }
                // if(npe_i > 0)
                //     unique_doms.insert(dom_i);
            }
        // }
        // std::set<int> unique_doms(dom_ids.begin(), dom_ids.end());
        int ntot_denoms = unique_doms.size();
        int idx = ntot_denoms >= 15 ? 14 : ntot_denoms-1;
        // if(ntot_denoms > 0)
            hnDoms->Fill(ntot_denoms, rock_wgt);
        if(idx >= 0)
            hdoms_eff[idx]->Fill(vtxX, vtxY, vtxZ, rock_wgt);
        for(int dom: unique_doms){
            int dom_x = (int)dom/18;
            int dom_y = dom % 18;
            hdom->Fill(dom_x, dom_y, 1.*rock_wgt/nevents);
        }
        hxy->Fill(vtxX, vtxY, unique_doms.size()*rock_wgt);
        hyz->Fill(vtxZ, vtxY, unique_doms.size()*rock_wgt);
        hlenxy->Fill(lenX, lenY, unique_doms.size()*rock_wgt);
        hxz->Fill(vtxZ, vtxX, unique_doms.size()*rock_wgt);
        // Let's try 5 DOM cut as a trigger
        if(unique_doms.size() > 3) {
            he_xy3->Fill(vtxX, vtxY, rock_wgt);
        }
        if(unique_doms.size() > 4) {
            he_xy4->Fill(vtxX, vtxY, rock_wgt);
        }
        if(unique_doms.size() > 5) {
            he_xy5->Fill(vtxX, vtxY, rock_wgt);
        }
        if(unique_doms.size() > 6) {
            he_xy6->Fill(vtxX, vtxY, rock_wgt);
        }
        if(unique_doms.size() > 7) {
            he_xy7->Fill(vtxX, vtxY, rock_wgt);
        }
        if(vtxY > 12000.)
            hmuE->Fill(muE*1.E-3, ntot_denoms, rock_wgt);
    }

    hxy->Divide(hTotxy);
    hyz->Divide(hTotyz);
    hlenxy->Divide(hTotlenxy);
    hxz->Divide(hTotxz);

    std::cout << "Writing to File!" << std::endl;
    TFile* outFile = new TFile(outfilename.c_str(), "recreate");
    hnDoms->Write("hnDoms");
    hxy->Write("fid_xy");
    hyz->Write("fid_yz");
    hlenxy->Write("len_xy");
    hxz->Write("fid_xz");
    he_xy3->Write("fid_event3_xy");
    he_xy4->Write("fid_event4_xy");
    he_xy5->Write("fid_event5_xy");
    he_xy6->Write("fid_event6_xy");
    he_xy7->Write("fid_event7_xy");
    hTotlenxy->Write("len_tot_xy");
    hTotxz->Write("fid_tot_xz");
    hTotxy->Write("fid_tot_xy");
    hTotyz->Write("fid_tot_yz");
    hdoms_eff_denom->Write("fid_tot_xyz");
    for(auto h: hdoms_eff)
        h->Write();
    hmuE->Write();
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
