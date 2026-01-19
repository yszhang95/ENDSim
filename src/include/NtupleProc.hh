#ifndef __END_NtupleProc__
#define __END_NtupleProc__

#include <Config.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Root.hh>
#include <RAT/OutNtupleProc.hh>
#include <RAT/Processor.hh>
#include <string>

namespace END {

class NtupleProc : public RAT::OutNtupleProc {
 public:
  NtupleProc();
  ~NtupleProc();
  void AssignAdditionalAddresses() override;
  void AssignAdditionalMetaAddresses() override;
  void FillNoTriggerEvent(RAT::DS::Root*) override;
  void FillMeta() override;

protected:
  double geo_index;
  std::string geo_file;
  std::string experiment;
  double source_pos_x;
  double source_pos_y;
  double source_pos_z;
  double source_rot_x;
  double source_rot_y;
  double source_rot_z;
};

}  // namespace END
#endif
