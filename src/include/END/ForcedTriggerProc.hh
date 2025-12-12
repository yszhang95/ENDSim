#ifndef __RAT_ForcedTriggerProc__
#define __RAT_ForcedTriggerProc__


#include <RAT/DB.hh>
#include <RAT/WaveformSampler.hh>
#include <RAT/Processor.hh>
#include <string>

namespace END {

class ForcedTriggerProc : public RAT::Processor {
 public:
  ForcedTriggerProc();
  virtual ~ForcedTriggerProc(){};
  virtual RAT::Processor::Result DSEvent(RAT::DS::Root *ds);

  void BeginOfRun(RAT::DS::Run *run);

 protected:
  int fEventCounter;
  bool fSampleWaveforms;

  RAT::DBLinkPtr ldaq;

  RAT::WaveformSampler *fSampler;
  std::string fSamplerType;
};

}  // namespace END

#endif
