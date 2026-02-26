// extract_dom_positions.C
// Usage: root -l -b -q 'extract_dom_positions.C("input.root", "dom_positions.txt")'
//
// Reads the "doms" tree from a pmt_bdt output ROOT file and writes one line
// per unique dom_id:
//   dom_id  dom_x  dom_y  dom_z
// to the specified text file.

#include <TFile.h>
#include <TTree.h>
#include <map>
#include <tuple>
#include <fstream>
#include <iostream>

void extract_dom_positions(const char* infile =
    "/home/yousen/Documents/darpa_daq/ENDSim/examples/c++/"
    "pmt_bdt_signal_numu_100000_rootracker_randomVtx_aframe_spacing5m_hex.root",
    const char* outfile = "dom_positions.txt")
{
    TFile* f = TFile::Open(infile);
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open: " << infile << std::endl;
        return;
    }

    TTree* t = (TTree*)f->Get("doms");
    if (!t) {
        std::cerr << "Cannot find tree 'doms' in " << infile << std::endl;
        f->Close();
        return;
    }

    int   dom_id;
    float dom_x, dom_y, dom_z;

    t->SetBranchStatus("*",      0);
    t->SetBranchStatus("dom_id", 1);
    t->SetBranchStatus("dom_x",  1);
    t->SetBranchStatus("dom_y",  1);
    t->SetBranchStatus("dom_z",  1);

    t->SetBranchAddress("dom_id", &dom_id);
    t->SetBranchAddress("dom_x",  &dom_x);
    t->SetBranchAddress("dom_y",  &dom_y);
    t->SetBranchAddress("dom_z",  &dom_z);

    std::map<int, std::tuple<float,float,float>> dommap;

    Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
        t->GetEntry(i);
        if (dommap.find(dom_id) == dommap.end())
            dommap[dom_id] = std::make_tuple(dom_x, dom_y, dom_z);
    }
    f->Close();

    std::ofstream ofs(outfile);
    if (!ofs.is_open()) {
        std::cerr << "Cannot write: " << outfile << std::endl;
        return;
    }

    ofs << "# dom_id  dom_x  dom_y  dom_z\n";
    for (auto& kv : dommap) {
        ofs << kv.first
            << " " << std::get<0>(kv.second)
            << " " << std::get<1>(kv.second)
            << " " << std::get<2>(kv.second)
            << "\n";
    }

    std::cout << "Wrote " << dommap.size() << " DOM positions to " << outfile << std::endl;
}
