#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/LQAnalysis/include/TTbarFullhadRecoHypothesis.h"


/**  @short This file defines analysis modules which add discriminator values to reconstruction hypotheses
 *
 * As input, all modules need a vector of ReconstructionHypothesis which must be in the event
 * when the module is called. The name of this vector can be specified at time of construction.
 * Further parameter can be defined by passing a 'cfg' object which typically controls the name of the
 * discriminator in th ReconstructionHypothesis or the name of additional event input.
 * 
 * The convention is that smaller values of the discriminator flag are better (as e.g. naturally the case
 * for chi-square). The exact meaning of the discriminators depends on the method, though.
 */


/** \brief Get the best hypothesis, i.e. the one with the smallest discriminator value
 * 
 * If no hypothesis exists with that name or if the smallest discriminator is infinite, returns
 * nullptr.
 * 
 * label is the discriminator label, e.g. "Chi2".
 */
const TTbarFullhadRecoHypothesis * get_best_hypothesis(const std::vector<TTbarFullhadRecoHypothesis> & hyps, const std::string & label);


/** \brief Calculate the chi-square reconstruction discriminator
 * 
 * The Chi-square value is calculated from leptonic and hadronic reconstructed top-quark masses. This
 * is the default reconstruction-level method used in the 8TeV semi-leptonic high-mass CMS analyses.
 * 
 * Per default, fills discriminators "Chi2", "Chi2_tlep" and "Chi2_thad" which are the over chi-square,
 * the chi-square only for the leptonic lep and the chi-square only for the hadronic leg, resp.
 * The name / prefix "Chi2" can be overridden via cfg::discriminator_label.
 * 
 * For numeric values of the means and widths for the masses used see the implementation in the .cxx file;
 * they are the 8TeV values.
 */



class TTbarFullhadRecoChi2Discriminator: public uhh2::AnalysisModule {
 public:
  struct cfg {
    std::string discriminator_label;
  cfg(): discriminator_label("Chi2Hadronic"){}
  };
    
  TTbarFullhadRecoChi2Discriminator(uhh2::Context & ctx, const std::string & rechyps_name, const cfg & config = cfg());
  virtual bool process(uhh2::Event & event) override;
    
 private:
  uhh2::Event::Handle<std::vector<TTbarFullhadRecoHypothesis>> h_hyps;
  cfg config;
};
