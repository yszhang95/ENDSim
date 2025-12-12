#include <CLHEP/Random/RandGauss.h>

#include <END/DS/DiscreteSignal.hh>
#include <END/WaveformSampler.hh>

namespace END {

WaveformSampler::WaveformSampler(std::string samplerName) { SetSamplerType(samplerName); }

void WaveformSampler::SetSamplerType(std::string samplerName) {
  fSamplerName = samplerName;
  fLdaq = DB::Get()->GetLink("SAMPLER", fSamplerName);

  fSamplingRate = fLdaq->GetD("sampling_rate");  // In GHz
  fOffset = fLdaq->GetD("offset");               // vertical offset in mV
  fVhigh = fLdaq->GetD("volt_high");             // in mV
  fVlow = fLdaq->GetD("volt_low");               // in mV
  fNoiseAmpl = fLdaq->GetD("noise_amplitude");  // discrete signal noise, in mV
  fNSamples = fLdaq->GetD("nsamples");
  fTerminationOhms = fLdaq->GetD("termination_ohms");

  detail << dformat("  WaveformSampler: Channel Noise: ............ %6.2f volts\n", fNoiseAmpl);
  detail << dformat("  WaveformSampler: Sampling Rate: ............ %6.2f ns\n", fSamplingRate);
  detail << dformat("  WaveformSampler: Total Number of Samples: .. %d \n", fNSamples);
  detail << dformat("  WaveformSampler: Voltage offset: ........... %6.2f mV\n", fOffset);
}

void WaveformSampler::AddWaveformGenerator(std::string modelName) {
  fPMTWaveformGenerators[modelName] = new PMTWaveformGenerator(modelName);
}

void WaveformSampler::SamplePMT(RAT::DS::MCPMT* mcpmt, int pmtID, double triggerTime, RAT::DS::PMTInfo* pmtinfo) {
  PMTWaveform pmtwfm = fPMTWaveformGenerators[pmtinfo->GetModelNameByID(pmtID)]->GenerateWaveforms(mcpmt, triggerTime);
  AddChannel(pmtID, pmtwfm);
}

void WaveformSampler::ClearWaveforms() { fDiscreteSignal.clear(); }

void WaveformSampler::WriteToEvent(DS::EV* ev) {
  DS::DiscreteSignal discreteSignal;

  auto waveforms = fDiscreteSignal;
  for (auto it = waveforms.begin(); it != waveforms.end(); it++) {
    discreteSignal.SetWaveform(it->first, waveforms[it->first]);
  }

  discreteSignal.SetSamplerName(fSamplerName);
  discreteSignal.SetNSamples(uint32_t(fNSamples));
  discreteSignal.SetDynamicRange((fVhigh - fVlow));
  discreteSignal.SetSamplingRate(fSamplingRate);
  discreteSignal.SetTerminationOhms(fTerminationOhms);

  ev->SetWaveformSampler(discreteSignal);
  ClearWaveforms();
}

// Add channel to sampler and immdediatly sample analogue waveform
void WaveformSampler::AddChannel(int ichannel, PMTWaveform pmtwf) {
  // Reset
  fDiscreteSignal[ichannel].clear();

  double timeres = 1.0 / fSamplingRate;  // in ns

  // Second, sample on analogue waveform
  double currenttime = 0;
  for (int isample = 0; isample < fNSamples; isample++) {
    double voltage = pmtwf.GetHeight(currenttime);
    voltage += fNoiseAmpl * CLHEP::RandGauss::shoot();           // add electronic noise

    // Manage voltage saturation
    if (voltage < fVlow) {
      voltage = fVlow;
    }
    else if (voltage > fVhigh) {
      voltage = fVhigh;
    }

    // Save sample
    fDiscreteSignal[ichannel].push_back(voltage);

    // Step on time
    currenttime += timeres;
  }
}
}  // namespace END
