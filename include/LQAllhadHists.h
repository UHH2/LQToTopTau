#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/LQAnalysis/include/TTbarFullhadRecoHypothesisDiscriminators.h"
#include "UHH2/LQAnalysis/include/TTbarFullhadRecoHypothesis.h"


/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class LQAllhadHists: public uhh2::Hists {
 public:
  // use the same constructor arguments as Hists for forwarding:
  LQAllhadHists(uhh2::Context & ctx, const std::string & dirname);
  bool is_data;
  virtual void fill(const uhh2::Event & ev) override;

 protected:
  uhh2::Event::Handle<std::vector<TTbarFullhadRecoHypothesis>> h_hyps;
  uhh2::Event::Handle<std::vector<TTbarFullhadRecoHypothesis>> h_hadr_hyps;


  virtual ~LQAllhadHists();
};
