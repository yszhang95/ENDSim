#include <NtupleProc.hh>
#include <RAT/DB.hh>
#include <RAT/DBLink.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Root.hh>
#include <RAT/Log.hh>
#include <RAT/Processor.hh>
#include <RAT/Rat.hh>

namespace END {

NtupleProc::NtupleProc() : RAT::OutNtupleProc() {
  // The registered this->name is used by the macro when calling
  // /rat/proc endntuple
  this->name = "endntuple";
}

NtupleProc::~NtupleProc() { this->FillMeta();  }

void NtupleProc::AssignAdditionalAddresses() {}

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

}  // namespace END
