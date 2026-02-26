#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TMath.h>

using namespace std;

static const int NPMTS         = 2790;
static const int NPMTS_PER_DOM = 31;
static const int NDOMS         = NPMTS / NPMTS_PER_DOM;  // 90

int ENDpmt2dom(int pmt_id) { return pmt_id / NPMTS_PER_DOM; }

// Compact summary of one frame entry (one PMT's ambient hits in a time window)
struct FrameData {
    int   npe;     // number of photo-electrons (hits) in this frame
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

void process(const std::string& infilename,
             const std::string& outfilename,
             const std::string& domposfile,
             int M,
             unsigned int seed)
{
    // ------------------------------------------------------------------
    // Open input and set up the frame_light tree
    // ------------------------------------------------------------------
    TFile* fin = new TFile(infilename.c_str(), "read");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open input file: " << infilename << std::endl;
        exit(1);
    }
    TTree* tree = (TTree*)fin->Get("frame_light");
    if (!tree) {
        std::cerr << "Cannot find tree 'frame_light' in " << infilename << std::endl;
        exit(1);
    }

    // Only read the branch we need
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("mcpehittime_rel", 1);

    vector<double>* mcpehittime_rel = nullptr;
    TBranch* b_mcpehittime_rel = nullptr;
    tree->SetBranchAddress("mcpehittime_rel", &mcpehittime_rel, &b_mcpehittime_rel);

    Long64_t nentries = tree->GetEntries();
    std::cout << "NFrames (entries): " << nentries << std::endl;

    // ------------------------------------------------------------------
    // Pre-load all frame entries
    // Each entry corresponds to one PMT's hits in a single time window.
    // No PMT ID filtering needed — all hits in the entry belong to one PMT.
    // ------------------------------------------------------------------
    vector<FrameData> frames(nentries);
    for (Long64_t e = 0; e < nentries; e++) {
        tree->GetEntry(e);
        frames[e].npe     = (int)mcpehittime_rel->size();
        frames[e].t_first = 0.f;
        if (frames[e].npe > 0)
            frames[e].t_first = (float)*min_element(mcpehittime_rel->begin(),
                                                     mcpehittime_rel->end());
        if (e % 100000 == 0)
            std::cout << "  Loading entry " << e << " / " << nentries << "\r" << std::flush;
    }
    std::cout << "\nPre-loading done." << std::endl;
    fin->Close();

    // ------------------------------------------------------------------
    // Load DOM geometry database
    // ------------------------------------------------------------------
    auto dompos = LoadDomPositions(domposfile);

    // ------------------------------------------------------------------
    // Set up output tree (same branch structure as pmt_bdt.cc doms tree,
    // minus muon-truth branches)
    // ------------------------------------------------------------------
    TFile* fout = new TFile(outfilename.c_str(), "recreate");
    TTree* tout = new TTree("doms", "doms");

    int   event_id, dom_id, du_id;
    float dom_x, dom_y, dom_z;
    int   npmts, npe, pe_min, pe_spread, pe_max;
    int   npmts_50, npmts_100, npe_50, npe_100;
    float pe_mean, pe_rms;
    float t_min, t_spread, t_mean, t_rms;

    tout->Branch("event_id",   &event_id,   "event_id/I");
    tout->Branch("dom_id",     &dom_id,     "dom_id/I");
    tout->Branch("du_id",      &du_id,      "du_id/I");
    tout->Branch("dom_x",      &dom_x,      "dom_x/F");
    tout->Branch("dom_y",      &dom_y,      "dom_y/F");
    tout->Branch("dom_z",      &dom_z,      "dom_z/F");
    tout->Branch("npmts",      &npmts,      "npmts/I");
    tout->Branch("npe",        &npe,        "npe/I");
    tout->Branch("pe_min",     &pe_min,     "pe_min/I");
    tout->Branch("pe_spread",  &pe_spread,  "pe_spread/I");
    tout->Branch("pe_rms",     &pe_rms,     "pe_rms/F");
    tout->Branch("t_min",      &t_min,      "t_min/F");
    // tout->Branch("t_spread",   &t_spread,   "t_spread/F");
    tout->Branch("t_mean",     &t_mean,     "t_mean/F");
    tout->Branch("t_rms",      &t_rms,      "t_rms/F");
    tout->Branch("npmts_50",   &npmts_50,   "npmts_50/I");
    tout->Branch("npe_50",     &npe_50,     "npe_50/I");
    tout->Branch("npmts_100",  &npmts_100,  "npmts_100/I");
    tout->Branch("npe_100",    &npe_100,    "npe_100/I");

    // ------------------------------------------------------------------
    // M-iteration loop
    // ------------------------------------------------------------------
    TRandom3 rng(seed);

    for (int m = 0; m < M; m++) {
        event_id = m;

        // Generate random frame assignment: seq[pmt_id] = frame index
        vector<Long64_t> seq(NPMTS);
        for (int i = 0; i < NPMTS; i++)
            seq[i] = (Long64_t)(rng.Uniform() * nentries);

        // Build per-DOM hit lists
        // dom_pes[d]   : npe per hit PMT in DOM d  (one entry per hit PMT)
        // dom_times[d] : t_first per hit PMT in DOM d
        vector<vector<int>>   dom_pes(NDOMS);
        vector<vector<float>> dom_times(NDOMS);

        for (int pmt_id = 0; pmt_id < NPMTS; pmt_id++) {
            const FrameData& fd = frames[seq[pmt_id]];
            if (fd.npe >= 1) {
                int d = ENDpmt2dom(pmt_id);
                dom_pes[d].push_back(fd.npe);
                dom_times[d].push_back(fd.t_first);
            }
        }

        // Per-DOM analysis — identical to pmt_bdt.cc lines 228-291
        for (int d = 0; d < NDOMS; d++) {
            if (dom_pes[d].empty()) continue;

            dom_id = d;
            du_id  = d / 18;

            // Look up DOM position from geometry database
            auto it = dompos.find(d);
            if (it == dompos.end()) {
                std::cerr << "DOM id " << d << " not found in geometry database" << std::endl;
                throw std::runtime_error("Missing DOM position for dom_id " + std::to_string(d));
            }
            dom_x = std::get<0>(it->second);
            dom_y = std::get<1>(it->second);
            dom_z = std::get<2>(it->second);

            npmts     = (int)dom_pes[d].size();
            npmts_50  = 0;
            npmts_100 = 0;
            npe_50    = 0;
            npe_100   = 0;

            vector<int>&   dom_pe   = dom_pes[d];
            vector<float>& dom_time = dom_times[d];

            sort(dom_pe.begin(),   dom_pe.end());
            sort(dom_time.begin(), dom_time.end());

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

            npe      = accumulate(dom_pe.begin(),   dom_pe.end(),   0);
            t_mean   = accumulate(dom_time.begin(), dom_time.end(), 0.f) / dom_time.size();
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

        if ((m + 1) % 100 == 0 || m == M - 1)
            std::cout << "Event " << m + 1 << " / " << M << "\r" << std::flush;
    }
    std::cout << "\nWriting to file!" << std::endl;

    tout->Write();
    fout->Close();
}


int main(int argc, char* argv[]) {
    if (argc < 5 || argc > 6) {
        std::cout << "Usage: " << argv[0]
                  << " <input_root> <output_root> <dom_positions> <M> [seed]" << std::endl;
        std::cout << "  input_root    : path to output_ambient_light_10us.root" << std::endl;
        std::cout << "  output_root   : output file path" << std::endl;
        std::cout << "  dom_positions : text file with DOM positions (dom_id x y z)" << std::endl;
        std::cout << "  M             : number of sampling iterations" << std::endl;
        std::cout << "  seed          : optional TRandom3 seed (default: 0 = time-based)" << std::endl;
        return 1;
    }

    std::string infilename(argv[1]);
    std::string outfilename(argv[2]);
    std::string domposfile(argv[3]);
    int M = atoi(argv[4]);
    unsigned int seed = (argc == 6) ? (unsigned int)atoi(argv[5]) : 0;

    process(infilename, outfilename, domposfile, M, seed);
    return 0;
}
