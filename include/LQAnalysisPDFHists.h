#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/LQAnalysis/include/TTbarFullhadRecoHypothesisDiscriminators.h"
#include "UHH2/LQAnalysis/include/TTbarFullhadRecoHypothesis.h"
#include "UHH2/common/include/TTbarGen.h"



/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class LQAnalysisPDFHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    LQAnalysisPDFHists(uhh2::Context & ctx, const std::string & dirname, bool use_pdf_weights_ = false);

    virtual void fill(const uhh2::Event & ev) override;
    std::string histo_names[100];
    std::string histo_names2[100];
    std::string histo_names3[100];
    std::string histo_names4[100];

  protected:
    bool use_pdf_weights;
    bool is_mc;
    uhh2::Event::Handle<std::vector<TTbarFullhadRecoHypothesis>> h_hyps;
    uhh2::Event::Handle<std::vector<TTbarFullhadRecoHypothesis>> h_hadr_hyps;

    virtual ~LQAnalysisPDFHists();
};


