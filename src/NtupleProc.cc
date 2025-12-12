#include <NtupleProc.hh>
#include <RAT/DB.hh>
#include <RAT/DBLink.hh>
#include <RAT/DS/DiscreteSignal.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/Log.hh>
#include <RAT/Processor.hh>
#include <RAT/Rat.hh>
#include <TTree.h>
#include <TFile.h>

namespace END {

NtupleProc::NtupleProc() : RAT::OutNtupleProc() {
  // The registered this->name is used by the macro when calling
  // /rat/proc endntuple
  this->name = "endntuple";
  RAT::DB *db = RAT::DB::Get();
  auto dblink = db->GetLink("ENDNTUPLE", "end_ntuple");
  // std::cerr << "NtupleProc: Getting config from ENDNTUPLE:end_ntuple" << std::endl;
  // this->fSaveDiscreteSignal = dblink->GetZ("save_discrete_signal");


    try {
    RAT::DBLinkPtr dblink = db->GetLink("ENDNTUPLE", "end_ntuple");
    std::cerr << "NtupleProc: Getting config from ENDNTUPLE:end_ntuple" << std::endl;

    try {
      fSaveDiscreteSignal = dblink->GetZ("save_discrete_signal");
      std::cerr << "NtupleProc: save_discrete_signal = "
                << (fSaveDiscreteSignal ? "true" : "false") << std::endl;
    } catch (RAT::DBNotFoundError &e) {
      RAT::warn << "ENDNTUPLE[end_ntuple] has no 'save_discrete_signal', "
                << "defaulting to false.\n";
      fSaveDiscreteSignal = false;
    }
  } catch (RAT::DBNotFoundError &e) {
    RAT::warn << "ENDNTUPLE[end_ntuple] not found in DB, "
              << "discrete signal ntuple disabled.\n";
    fSaveDiscreteSignal = false;
  }
}

NtupleProc::~NtupleProc() { this->FillMeta();  }

void NtupleProc::AssignAdditionalAddresses() {
  if (this->fSaveDiscreteSignal) {
    dwfTree = new TTree("dwf", "Discrete waveform");
    dwfTree->Branch("evid", &evid);
    dwfTree->Branch("pmtid", &dwf_pmtid);
    dwfTree->Branch("inWindowPulseTimes", &dwf_inWindowPulseTimes);
    dwfTree->Branch("inWindowPulseCharges", &dwf_inWindowPulseCharges);
    dwfTree->Branch("waveform", &dwf_waveform);
  }
}

void NtupleProc::AssignAdditionalMetaAddresses() {
  this->metaTree->Branch("geo_index", &geo_index);
  this->metaTree->Branch("geo_file", &geo_file);
  this->metaTree->Branch("experiment", &experiment);
  this->metaTree->Branch("source_pos_x", &source_pos_x);
  this->metaTree->Branch("source_pos_y", &source_pos_y);
  this->metaTree->Branch("source_pos_z", &source_pos_z);
  this->metaTree->Branch("source_rot_x", &source_rot_x);
  this->metaTree->Branch("source_rot_y", &source_rot_y);
  this->metaTree->Branch("source_rot_z", &source_rot_z);
}

void NtupleProc::FillEvent(RAT::DS::Root *ds, RAT::DS::EV *ev) {
  RAT::OutNtupleProc::FillEvent(ds, ev);
  if (this->fSaveDiscreteSignal) {

    RAT::DS::MC *mc = ds->GetMC();
    runBranch = RAT::DS::RunStore::GetRun(ds);
    RAT::DS::PMTInfo *pmtinfo = runBranch->GetPMTInfo();
    const RAT::DS::ChannelStatus *channel_status = runBranch->GetChannelStatus();

    RAT::DS::DiscreteSignal sampler = ev->GetWaveformSampler();
    double readout_window_min = 0;
    double readout_window_max = sampler.GetNSamples() * sampler.GetTimeStepNS();
      for (auto const &pair : sampler.GetAllWaveforms()) {
        dwf_pmtid = pair.first;
        dwf_waveform = pair.second;
        dwf_inWindowPulseTimes.clear();
        dwf_inWindowPulseCharges.clear();
        if (mc->GetMCPMTCount() == 0) {
          // std::cerr << "No MC information found, skipping MCPhoton matching.\n";
        } // if there's no MC information, skip it
        else if (dwf_pmtid < 0) {
        }  // these are nonPMT channels. No MC info
        else {
          auto it = std::find(mcpmtid.begin(), mcpmtid.end(), dwf_pmtid);
          if (it == mcpmtid.end())
            RAT::warn << "No MC information found for PMTID = " << dwf_pmtid
                 << " but waveform exists for some reason..." << newline;
          else {
            double time_offset = channel_status->GetCableOffsetByPMTID(dwf_pmtid);
            RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(it - mcpmtid.begin());
            for (int ipe = 0; ipe < mcpmt->GetMCPhotonCount(); ipe++) {
              RAT::DS::MCPhoton *mcph = mcpmt->GetMCPhoton(ipe);
              Double_t time = mcph->GetFrontEndTime() - ev->GetCalibratedTriggerTime() + time_offset;
              Double_t charge = mcph->GetCharge();
              if (time > readout_window_min && time < readout_window_max) {
                dwf_inWindowPulseTimes.push_back(time);
                dwf_inWindowPulseCharges.push_back(charge);
              }
              // std::cerr << "PMTID " << dwf_pmtid << " MCPhoton " << ipe
              //           << " time " << time << " charge " << charge << "\n";
              // std::cerr << "    front_end_time " << mcph->GetFrontEndTime()
              //           << " cal_trig_time " << ev->GetCalibratedTriggerTime()
              //           << " cable_offset " << time_offset << "\n";
              // std::cerr << "    readout_window_min " << readout_window_min
              //           << " readout_window_max " << readout_window_max << "\n";
            } // END loop over mcphotons
            // std::cerr << "-------------------------- PMTID " << dwf_pmtid
            //           << " total mcphotons " << mcpmt->GetMCPhotonCount() << "\n";
          }
          // std::cerr << "-------------------------- dfw_pmtid " << dwf_pmtid
          //           << " " << dwf_inWindowPulseTimes.size()
          //           << " pulses in window\n";

          // for (size_t i=0; i<dwf_waveform.size(); i++) {
          //   if (i<10 || i>(dwf_waveform.size()-10)) {
          //     std::cerr << dwf_waveform[i] << " ";
          //   }
          //   if (i==10) {
          //     std::cerr << "... ";
          //   }
          // }
          // std::cerr << "\n";
        }  // END IF
        dwfTree->Fill();
      }
  }
}

void NtupleProc::FillNoTriggerEvent(RAT::DS::Root* ds) {}

void NtupleProc::FillMeta() {
  RAT::DB* db = RAT::DB::Get();
  RAT::DBLinkPtr ldetector = db->GetLink("DETECTOR");
  try {
    experiment = ldetector->GetS("experiment");
  } catch (RAT::DBNotFoundError& e) {
    experiment = "";
    RAT::info << "Experiment not found." << newline;
  }

  try {
    geo_file = ldetector->GetS("geo_file");
  } catch (RAT::DBNotFoundError& e) {
    geo_file = "";
    RAT::info << "Geometry file not found." << newline;
  }

  try {
    geo_index = ldetector->GetD("geo_index");
  } catch (RAT::DBNotFoundError& e) {
    geo_index = -9999;
    RAT::info << "Geometry index not found." << newline;
  }

  // This only works if the cal source is named specifically
  try {
    RAT::DBLinkPtr lsource = db->GetLink("GEO", "source_mother");
    std::vector<double> source_pos = lsource->GetDArray("position");
    std::vector<double> source_rot = lsource->GetDArray("rotation");
    source_pos_x = source_pos[0];
    source_pos_y = source_pos[1];
    source_pos_z = source_pos[2];
    source_rot_x = source_rot[0];
    source_rot_y = source_rot[1];
    source_rot_z = source_rot[2];
  } catch (RAT::DBNotFoundError& e) {
    source_pos_x = -9999;
    source_pos_y = -9999;
    source_pos_z = -9999;
    source_rot_x = -9999;
    source_rot_y = -9999;
    source_rot_z = -9999;
    RAT::info << "Source location/rotation not identified." << newline;
  }
}

void NtupleProc::EndOfRun(RAT::DS::Run* run) {
  if (this->outputFile) {
    this->outputFile->cd();
    if (this->fSaveDiscreteSignal) {
      std::cerr << "----------------------------------- write to output\n";
      dwfTree->Write();
    }
  }
  OutNtupleProc::EndOfRun(run);
}
}  // namespace END
