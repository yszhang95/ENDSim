/**
 * @class DS::EV
 * Data Structure: Triggered event
 *
 * @author Yousen Zhang <yzhang11@bnl.gov>
 *
 * This class represents a detected event as defined by the trigger
 * simulation. It extends the RAT::DS::EV class by adding support for
 * waveform samplers.
 */

#ifndef __END_DS_EV__
#define __END_DS_EV__
#ifdef __MAKECINT__
#define __signed signed
#endif

#include <RAT/DS/EV.hh>
#include <END/DS/DiscreteSignal.hh>

namespace END {
namespace DS {

class EV : public RAT::DS::EV {
 public:
  EV() : RAT::DS::EV() {}
  virtual ~EV() {}

  // Setters
  void SetWaveformSampler(const DiscreteSignal &ds) { fWaveformSampler.push_back(ds); }

  // Getters
  const DiscreteSignal &GetWaveformSampler() const { return fWaveformSampler.at(0); }

  // Check if the sampler exists
  virtual bool SamplerExists() const { return !fWaveformSampler.empty(); }

  // Prune sampler information
  virtual void PruneSampler() { fWaveformSampler.clear(); fWaveformSampler.shrink_to_fit(); }

  ClassDef(EV, 6);

 protected:
   std::vector<DiscreteSignal> fWaveformSampler;
};
} // namespace DS
} // namespace END

#endif
