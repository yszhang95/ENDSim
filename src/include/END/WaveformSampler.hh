////////////////////////////////////////////////////////////////////
/// \class END::WaveformSampler
///
/// \brief   Waveform sampler class for discrete signal
///
/// \author Yousen Zhang <yzhang11@bnl.gov>
///
/// REVISION HISTORY:\n
///     23 Nov 2025: Initial commit
///
/// \details
/// This class provides full support for generating a discrete signal.
/// It samples PMTWaveforms at fixed intervals with a given sampling rate,
/// number of samples. It also adds electronic noise.
////////////////////////////////////////////////////////////////////
#ifndef __END_WaveformSampler__
#define __END_WaveformSampler__

#include <RAT/DB.hh>
#include <END/DS/DiscreteSignal.hh>
#include <END/DS/EV.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/PMTWaveformGenerator.hh>
#include <map>

namespace END {

class WaveformSampler {
 public:
  WaveformSampler(){};
  virtual ~WaveformSampler(){};
  WaveformSampler(std::string);

  virtual void SetSamplerType(std::string);
  virtual void SamplePMT(RAT::DS::MCPMT *mcpmt, int pmtID, double triggerTime, RAT::DS::PMTInfo *pmtinfo);
  virtual void ClearWaveforms();
  virtual void WriteToEvent(DS::EV *ev);
  virtual void AddChannel(int ichannel, PMTWaveform pmtwf);

  void AddWaveformGenerator(std::string modelName);

  std::string fSamplerName;   // Sampler type
  double fSamplingRate;     // Sampling rate in GHz
  int fNSamples;            // Total number of samples per sampling
  double fTerminationOhms;  // Input impedence of the sampler
  // Channel:Discrete waveform for each channel
  std::map<int, std::vector<Float_t>> fDiscreteSignal;

  std::map<std::string, PMTWaveformGenerator *> fPMTWaveformGenerators;

 protected:
  DBLinkPtr fLdaq;
  DBLinkPtr fSampledSiganl;

  double fOffset; // Sampler offset
  double fVhigh;  // Maximum voltage allowed in the sampler
  double fVlow; // Minimum voltage allowed in the sampler
  double fNoiseAmpl;  // Electronic noise width
};

}  // namespace END

#endif
