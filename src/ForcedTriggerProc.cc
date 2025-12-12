/*
Forced trigger causes a trigger regardless of the total number of PMTs hit
 */
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/RunStore.hh>
#include <END/ForcedTriggerProc.hh>
#include <cmath>

namespace END {

ForcedTriggerProc::ForcedTriggerProc() : Processor("endforcedtrigger") {}

void ForcedTriggerProc::BeginOfRun(RAT::DS::Run *run) {
  // Trigger Specifications
  ldaq = RAT::DB::Get()->GetLink("ENDDAQ", "endforcedtrigger");
  fEventCounter = 0;
  fSamplerType = ldaq->GetS("sampler_name");
  std::cerr << "---------------------------------------------sample_waveforms: " << ldaq->GetS("sampler_name") << "\n" << std::endl;
  std::cerr << "---------------------------------------------sample_waveforms: " << ldaq->GetZ("sample_waveforms")  << "\n" << std::endl;
  fSampleWaveforms = ldaq->GetZ("sample_waveforms");

  fSampler = new RAT::WaveformSampler(fSamplerType);

  if (fSampleWaveforms) {
    RAT::DS::PMTInfo *pmtinfo = run->GetPMTInfo();
    const size_t numModels = pmtinfo->GetModelCount();
    for (size_t i = 0; i < numModels; i++) {
      const std::string &modelName = pmtinfo->GetModelName(i);
      fSampler->AddWaveformGenerator(modelName);
    }
  }
}

RAT::Processor::Result ForcedTriggerProc::DSEvent(RAT::DS::Root *ds) {
  RAT::DS::MC *mc = ds->GetMC();
  RAT::DS::Run *run = RAT::DS::RunStore::Get()->GetRun(ds);
  RAT::DS::PMTInfo *pmtinfo = run->GetPMTInfo();
  // Prune the previous EV branchs if one exists
  if (ds->ExistEV()) ds->PruneEV();

  RAT::DS::EV *ev = ds->AddNewEV();
  ev->SetID(fEventCounter++);
  ev->SetCalibratedTriggerTime(0.0);
  ev->SetUTC(mc->GetUTC());
  double totalEVCharge = 0;
  // Loop over the mcpmts and fill the pmt branch assuming we've triggered
  for (int imcpmt = 0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
    RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
    mcpmt->SortMCPhotons();
    int pmtID = mcpmt->GetID();
    double integratedCharge = 0;
    double time = 0;
    for (int pidx = 0; pidx < mcpmt->GetMCPhotonCount(); pidx++) {
      RAT::DS::MCPhoton *photon = mcpmt->GetMCPhoton(pidx);
      // Use the first PE time as the pmt time
      if (pidx == 0) {
        time = photon->GetFrontEndTime();
      }
      integratedCharge += photon->GetCharge();
    }
    RAT::DS::PMT *pmt = ev->GetOrCreatePMT(pmtID);
    pmt->SetTime(time);
    pmt->SetCharge(integratedCharge);
    totalEVCharge += integratedCharge;
    if (fSampleWaveforms) {
      fSampler->SamplePMT(mcpmt, pmtID, 0.0, pmtinfo);
    }
  }
  if (fSampleWaveforms) {
    fSampler->WriteToEvent(ev);
  }
  ev->SetTotalCharge(totalEVCharge);

  return RAT::Processor::OK;
}

}  // namespace END
