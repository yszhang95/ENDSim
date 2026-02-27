#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <numeric>

#include <RAT/DSReader.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCSummary.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/EV.hh>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TMath.h>

#include <glob.h>

using namespace std;

static const int NPMTS         = 2790;
static const int NPMTS_PER_DOM = 31;
static const int NDOMS         = NPMTS / NPMTS_PER_DOM;  // 90

int ENDpmt2dom(int pmt_id) { return pmt_id / NPMTS_PER_DOM; }

// Compact summary of one ambient frame entry (one PMT's hits in a time window)
struct FrameData {
    int   npe;     // number of photo-electrons in this frame
    float t_first; // earliest hit time relative to frame start (ns)
};

// Load DOM positions from a text file (dom_id dom_x dom_y dom_z per line).
// Lines starting with '#' are treated as comments and skipped.
std::map<int, std::tuple<float,float,float>> LoadDomPositions(const std::string& path)
{
    std::map<int, std::tuple<float,float,float>> m;
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        std::cerr << "Cannot open DOM positions file: " << path << std::endl;
        exit(1);
    }
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int   id;
        float x, y, z;
        if (ss >> id >> x >> y >> z)
            m[id] = std::make_tuple(x, y, z);
    }
    std::cout << "Loaded " << m.size() << " DOM positions from " << path << std::endl;
    return m;
}

void process(const std::string& signal_pattern,
             const std::string& ambient_file,
             const std::string& domposfile,
             const std::string& outfilename,
             unsigned int seed)
{
    // ------------------------------------------------------------------
    // Load ambient frame pool
    // ------------------------------------------------------------------
    TFile* fin_amb = new TFile(ambient_file.c_str(), "read");
    if (!fin_amb || fin_amb->IsZombie()) {
        std::cerr << "Cannot open ambient file: " << ambient_file << std::endl;
        exit(1);
    }
    TTree* tree_amb = (TTree*)fin_amb->Get("frame_light");
    if (!tree_amb) {
        std::cerr << "Cannot find tree 'frame_light' in " << ambient_file << std::endl;
        exit(1);
    }

    tree_amb->SetBranchStatus("*", 0);
    tree_amb->SetBranchStatus("mcpehittime_rel", 1);

    vector<double>* mcpehittime_rel = nullptr;
    TBranch* b_mcpehittime_rel      = nullptr;
    tree_amb->SetBranchAddress("mcpehittime_rel", &mcpehittime_rel, &b_mcpehittime_rel);

    Long64_t nambient = tree_amb->GetEntries();
    std::cout << "NAmbient frames: " << nambient << std::endl;

    vector<FrameData> amb_frames(nambient);
    for (Long64_t e = 0; e < nambient; e++) {
        tree_amb->GetEntry(e);
        amb_frames[e].npe     = (int)mcpehittime_rel->size();
        amb_frames[e].t_first = 0.f;
        if (amb_frames[e].npe > 0)
            amb_frames[e].t_first = (float)*min_element(mcpehittime_rel->begin(),
                                                         mcpehittime_rel->end());
        if (e % 100000 == 0)
            std::cout << "  Loading ambient entry " << e << " / " << nambient << "\r" << std::flush;
    }
    std::cout << "\nAmbient pool loaded." << std::endl;
    fin_amb->Close();

    // ------------------------------------------------------------------
    // Load DOM geometry database
    // ------------------------------------------------------------------
    auto dompos = LoadDomPositions(domposfile);

    // ------------------------------------------------------------------
    // Load cross-section weights (same as pmt_bdt.cc)
    // ------------------------------------------------------------------
    TFile* xsec    = new TFile("xsec_graphs.root", "read");
    TGraph* gr_Si28 = (TGraph*)xsec->Get("nu_mu_Si28/tot_cc");
    TGraph* gr_O16  = (TGraph*)xsec->Get("nu_mu_O16/tot_cc");
    TGraph* gr_H1   = (TGraph*)xsec->Get("nu_mu_H1/tot_cc");
    xsec->Close();
    double xsec_weight = (gr_Si28->Eval(3.5) + 2.*gr_O16->Eval(3.5)) /
                         (gr_O16->Eval(3.5)  + 2.*gr_H1->Eval(3.5));

    // ------------------------------------------------------------------
    // Open signal files with RAT DSReader
    // ------------------------------------------------------------------
    glob_t glob_result;
    int gresult = glob(signal_pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if (gresult != 0) {
        std::cerr << "Couldn't find signal files! Exiting!" << std::endl;
        exit(1);
    }

    RAT::DSReader* dsreader = new RAT::DSReader(signal_pattern.c_str());
    for (int i = 0; i < (int)glob_result.gl_pathc; i++)
        dsreader->Add(glob_result.gl_pathv[i]);
    const unsigned int nevents = dsreader->GetT()->GetEntries();
    std::cout << "NSignal events: " << nevents << std::endl;

    // ------------------------------------------------------------------
    // Set up output tree (same schema as pmt_bdt.cc)
    // ------------------------------------------------------------------
    TH1D* hevts = new TH1D("nevts", "Total Evts", 1, 0, 1);
    TTree* tout  = new TTree("doms", "doms");

    int   event_id, dom_id, du_id;
    float dom_x, dom_y, dom_z;
    float vtxX, vtxY, vtxZ;
    float momX, momY, momZ;
    float muE, rock_wgt;
    int   npmts, npe, pe_min, pe_spread, pe_max;
    int   npmts_50, npmts_100, npe_50, npe_100;
    float pe_mean, pe_rms;
    float t_min, t_spread, t_mean, t_rms;

    tout->Branch("event_id",  &event_id,  "event_id/I");
    tout->Branch("dom_id",    &dom_id,    "dom_id/I");
    tout->Branch("du_id",     &du_id,     "du_id/I");
    tout->Branch("vtxX",      &vtxX,      "vtxX/F");
    tout->Branch("vtxY",      &vtxY,      "vtxY/F");
    tout->Branch("vtxZ",      &vtxZ,      "vtxZ/F");
    tout->Branch("momX",      &momX,      "momX/F");
    tout->Branch("momY",      &momY,      "momY/F");
    tout->Branch("momZ",      &momZ,      "momZ/F");
    tout->Branch("muE",       &muE,       "muE/F");
    tout->Branch("rock_wgt",  &rock_wgt,  "rock_wgt/F");
    tout->Branch("dom_x",     &dom_x,     "dom_x/F");
    tout->Branch("dom_y",     &dom_y,     "dom_y/F");
    tout->Branch("dom_z",     &dom_z,     "dom_z/F");
    tout->Branch("npmts",     &npmts,     "npmts/I");
    tout->Branch("npe",       &npe,       "npe/I");
    tout->Branch("pe_min",    &pe_min,    "pe_min/I");
    tout->Branch("pe_spread", &pe_spread, "pe_spread/I");
    tout->Branch("pe_rms",    &pe_rms,    "pe_rms/F");
    tout->Branch("t_min",     &t_min,     "t_min/F");
    tout->Branch("t_spread",  &t_spread,  "t_spread/F");
    tout->Branch("t_mean",    &t_mean,    "t_mean/F");
    tout->Branch("t_rms",     &t_rms,     "t_rms/F");
    tout->Branch("npmts_50",  &npmts_50,  "npmts_50/I");
    tout->Branch("npe_50",    &npe_50,    "npe_50/I");
    tout->Branch("npmts_100", &npmts_100, "npmts_100/I");
    tout->Branch("npe_100",   &npe_100,   "npe_100/I");

    // ------------------------------------------------------------------
    // RNG for ambient sampling and global time offset
    // ------------------------------------------------------------------
    TRandom3 rng(seed);

    // ------------------------------------------------------------------
    // Event loop
    // ------------------------------------------------------------------
    double nwgt_events = 0;
    for (size_t iev = 0; iev < nevents; iev++) {
        RAT::DS::Root* rds = dsreader->GetEvent(iev);
        if (!rds->ExistMC()) continue;
        if (rds->GetEVCount() == 0) continue;

        RAT::DS::MC* mc = rds->GetMC();
        RAT::DS::EV* ev = rds->GetEV(0);

        // Truth extraction (same as pmt_bdt.cc)
        int mcpcount = mc->GetMCParticleCount();
        vtxX = -100000.f;
        vtxY = -100000.f;
        vtxZ = -100000.f;
        momX = -100000.f;
        momY = -100000.f;
        momZ = -100000.f;
        muE  = -5.f;
        for (int pid = 0; pid < mcpcount; pid++) {
            RAT::DS::MCParticle* particle = mc->GetMCParticle(pid);
            if (abs(particle->GetPDGCode()) == 13) {
                TVector3 mcpos = particle->GetPosition();
                vtxX = mcpos.X();
                vtxY = mcpos.Y();
                vtxZ = mcpos.Z();
                TVector3 mcmom = particle->GetMomentum();
                momX = mcmom.X();
                momY = mcmom.Y();
                momZ = mcmom.Z();
                muE  = particle->GetKE() / 1000.f;
            }
        }

        // Rock upweighting (same as pmt_bdt.cc)
        rock_wgt = 1.f;
        if (((signal_pattern.find("rockbed") != string::npos) ||
             (signal_pattern.find("equalx")  != string::npos) ||
             (signal_pattern.find("aframe")  != string::npos)) &&
            (signal_pattern.find("cosmic") == string::npos)) {
            if (vtxY > -23000.f && vtxY < -13000.f &&
                (signal_pattern.find("short") == string::npos))
                rock_wgt = (float)xsec_weight;
            if (vtxY > -22500.f && vtxY < -12500.f &&
                (signal_pattern.find("short") != string::npos))
                rock_wgt = (float)xsec_weight;
            if (vtxY < -9000.f &&
                (signal_pattern.find("DU") != string::npos))
                rock_wgt = 0.f;
        }
        nwgt_events += rock_wgt;

        event_id = (int)iev;

        // Random global time offset for this event (ns)
        float toffset = (float)rng.Uniform(-10000., 40000.);

        // ------------------------------------------------------------------
        // Collect signal PMT hits (type==1 only, with toffset applied)
        // ------------------------------------------------------------------
        std::map<int, int>   sig_npe;   // pmt_id -> total signal photons
        std::map<int, float> sig_time;  // pmt_id -> trigger time + toffset

        for (int ipmt = 0; ipmt < mc->GetMCPMTCount(); ipmt++) {
            RAT::DS::MCPMT* mcpmt = mc->GetMCPMT(ipmt);
            int npe_i = mcpmt->GetMCPhotonCount();
            int pmt_i = mcpmt->GetID();
            int type  = mcpmt->GetType();
            if (type == 0) continue;  // skip HQE 8'' PMTs
            RAT::DS::PMT* pmt = ev->GetOrCreatePMT(pmt_i);
            if (npe_i >= 1) {
                sig_npe[pmt_i]  = npe_i;
                sig_time[pmt_i] = (float)pmt->GetTime() + toffset;
            }
        }

        // ------------------------------------------------------------------
        // Sample ambient light for all 2790 PMTs
        // ------------------------------------------------------------------
        vector<Long64_t> amb_seq(NPMTS);
        for (int i = 0; i < NPMTS; i++)
            amb_seq[i] = (Long64_t)(rng.Uniform() * nambient);

        // ------------------------------------------------------------------
        // Build per-DOM hit lists, keeping signal and ambient separate
        // ------------------------------------------------------------------
        std::map<int, std::vector<int>>   dom_pes_sig;
        std::map<int, std::vector<float>> dom_times_sig;
        std::map<int, std::vector<int>>   dom_pes_amb;
        std::map<int, std::vector<float>> dom_times_amb;

        for (int pmt_id = 0; pmt_id < NPMTS; pmt_id++) {
            const FrameData& af  = amb_frames[amb_seq[pmt_id]];
            bool has_sig = (sig_npe.count(pmt_id) > 0);
            bool has_amb = (af.npe >= 1);

            if (!has_sig && !has_amb) continue;

            int dom_i = ENDpmt2dom(pmt_id);

            if (has_sig) {
                dom_pes_sig[dom_i].push_back(sig_npe[pmt_id]);
                dom_times_sig[dom_i].push_back(sig_time[pmt_id]);
            }

            if (has_amb) {
                dom_pes_amb[dom_i].push_back(af.npe);
                dom_times_amb[dom_i].push_back(af.t_first);
            }
        }

        // ------------------------------------------------------------------
        // Collect all DOM IDs with either signal or ambient data
        // ------------------------------------------------------------------
        std::set<int> dom_ids_with_data;
        for (auto& entry : dom_pes_sig) dom_ids_with_data.insert(entry.first);
        for (auto& entry : dom_pes_amb) dom_ids_with_data.insert(entry.first);

        // ------------------------------------------------------------------
        // Per-DOM analysis: combine signal and ambient vectors
        // ------------------------------------------------------------------
        for (int dom_id : dom_ids_with_data) {
            du_id = dom_id / 18;

            auto it = dompos.find(dom_id);
            if (it == dompos.end()) {
                std::cerr << "DOM id " << dom_id << " not found in geometry database" << std::endl;
                throw std::runtime_error("Missing DOM position for dom_id " + std::to_string(dom_id));
            }
            dom_x = std::get<0>(it->second);
            dom_y = std::get<1>(it->second);
            dom_z = std::get<2>(it->second);

            // Combine signal and ambient into single vectors for analysis
            std::vector<int>   dom_pe;
            std::vector<float> dom_time;

            if (dom_pes_sig.count(dom_id)) {
                dom_pe.insert(dom_pe.end(), dom_pes_sig[dom_id].begin(), dom_pes_sig[dom_id].end());
                dom_time.insert(dom_time.end(), dom_times_sig[dom_id].begin(), dom_times_sig[dom_id].end());
            }

            if (dom_pes_amb.count(dom_id)) {
                dom_pe.insert(dom_pe.end(), dom_pes_amb[dom_id].begin(), dom_pes_amb[dom_id].end());
                dom_time.insert(dom_time.end(), dom_times_amb[dom_id].begin(), dom_times_amb[dom_id].end());
            }

            if (dom_pe.empty()) continue;

            npmts     = (int)dom_pe.size();
            npmts_50  = 0;
            npmts_100 = 0;
            npe_50    = 0;
            npe_100   = 0;

            std::sort(dom_pe.begin(),   dom_pe.end());
            std::sort(dom_time.begin(), dom_time.end());

            float init_time = dom_time.at(0);
            for (int i = 0; i < (int)dom_time.size(); i++) {
                float pmt_time = dom_time.at(i);
                int   pmt_pe   = dom_pe.at(i);
                if (pmt_time - init_time <= 50) {
                    npmts_50  += 1;
                    npe_50    += pmt_pe;
                    npmts_100 += 1;
                    npe_100   += pmt_pe;
                } else if (pmt_time - init_time <= 100) {
                    npmts_100 += 1;
                    npe_100   += pmt_pe;
                    npmts_50   = 0;
                    npe_50     = pmt_pe;
                } else {
                    npmts_50  = 0;
                    npe_50    = pmt_pe;
                    npmts_100 = 0;
                    npe_100   = pmt_pe;
                }
                init_time = pmt_time;
            }

            npe      = std::accumulate(dom_pe.begin(),   dom_pe.end(),   0);
            t_mean   = std::accumulate(dom_time.begin(), dom_time.end(), 0.f) / dom_time.size();
            pe_mean  = (float)npe / dom_pe.size();
            pe_min   = dom_pe.at(0);
            pe_max   = dom_pe.at(dom_pe.size() - 1);
            pe_spread = pe_max - pe_min;
            t_spread  = dom_time.at(dom_time.size() - 1) - dom_time.at(0);
            t_min     = dom_time.at(0);

            pe_rms = 0.f;
            t_rms  = 0.f;
            for (int i = 0; i < (int)dom_pe.size(); i++) {
                pe_rms += (float)TMath::Power(dom_pe[i]   - pe_mean, 2);
                t_rms  += (float)TMath::Power(dom_time[i] - t_mean,  2);
            }
            pe_rms = TMath::Sqrt(pe_rms) / dom_pe.size();
            t_rms  = TMath::Sqrt(t_rms)  / dom_time.size();

            tout->Fill();
        }

        if ((iev + 1) % 100 == 0 || iev == nevents - 1)
            std::cout << "Event " << iev + 1 << " / " << nevents << "\r" << std::flush;
    }

    hevts->SetBinContent(1, nwgt_events);
    std::cout << "\nWriting to file!" << std::endl;
    TFile* outFile = new TFile(outfilename.c_str(), "recreate");
    tout->Write();
    hevts->Write();
    outFile->Close();
}


int main(int argc, char* argv[]) {
    if (argc < 5 || argc > 6) {
        std::cout << "Usage: " << argv[0]
                  << " <signal_pattern> <ambient_root> <dom_positions> <output_root> [seed]" << std::endl;
        std::cout << "  signal_pattern : glob pattern for RAT signal files" << std::endl;
        std::cout << "  ambient_root   : path to output_ambient_light_10us.root" << std::endl;
        std::cout << "  dom_positions  : text file with DOM positions (dom_id x y z)" << std::endl;
        std::cout << "  output_root    : output file path" << std::endl;
        std::cout << "  seed           : optional TRandom3 seed (default: 0 = time-based)" << std::endl;
        return 1;
    }

    std::string signal_pattern(argv[1]);
    std::string ambient_file(argv[2]);
    std::string domposfile(argv[3]);
    std::string outfilename(argv[4]);
    unsigned int seed = (argc == 6) ? (unsigned int)atoi(argv[5]) : 0;

    process(signal_pattern, ambient_file, domposfile, outfilename, seed);
    return 0;
}
