#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/LQAnalysis/include/TTbarFullhadRecoHypothesis.h"
#include "TMinuit.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/TopJetIds.h"

typedef std::function< std::vector<LorentzVector>  (const LorentzVector & lepton, const LorentzVector & met)> NeutrinoReconstructionMethod;



/** \brief Make a list of ttbar reconstruction hypotheses as used in high-mass semileptonic ttbar 8 TeV analysis
 *
 * Make a list of ttbar reconstruction hypothesis using all (~ 3^Njet) combinations
 * of assigning jets to either the leptonic top, the hadronic top, or none of them;
 * hypotheses not assigning any jet to either the hadronic or leptonic top
 * are discarded.
 * For the leptonic side, the primary lepton and the neutrino reconstruction
 * according to the neutrinofunction parameter is done, which typically doubles
 * the number of hypotheses.
 *
 * Make sure to run an appropriate cleaner to keep only jets which should be used
 * in the hypothesis. Only works for events with Njets >= 2.
 *
 * neutrinofunction can be e.g. NeutrinoReconstruction or NeutrinoFitPolar
 *
 * label = name of the hypotheses list in the event / output tree
 */


class HighMassHadronicTTbarReco: public uhh2::AnalysisModule {
 public:

  explicit HighMassHadronicTTbarReco(uhh2::Context & ctx, const std::string & label="HighMassHadronicTTbarFullhadReco");

  virtual bool process(uhh2::Event & event) override;

  virtual ~HighMassHadronicTTbarReco();

 private:
  uhh2::Event::Handle<std::vector<TTbarFullhadRecoHypothesis>> h_hadr_recohyps;
};



std::vector<LorentzVector> LQNeutrinoReconstruction(const LorentzVector & lepton, const LorentzVector & met);

