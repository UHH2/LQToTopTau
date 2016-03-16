#include "UHH2/LQAnalysis/include/TTbarFullhadRecoHypothesisDiscriminators.h"
#include "UHH2/core/include/Utils.h"

#include <set>

using namespace uhh2;
using namespace std;

namespace {
    
// invariant mass of a lorentzVector, but save for timelike / spacelike vectors
float inv_mass(const LorentzVector & p4){
    if(p4.isTimelike()){
            return p4.mass();
    }
    else{
        return -sqrt(-p4.mass2());
    }
}

}


const TTbarFullhadRecoHypothesis * get_best_hypothesis(const std::vector<TTbarFullhadRecoHypothesis> & hyps, const std::string & label){
    const TTbarFullhadRecoHypothesis * best = nullptr;
    float current_best_disc = numeric_limits<float>::infinity();
    for(const auto & hyp : hyps){
        if(!hyp.has_discriminator(label)) continue;
        auto disc = hyp.discriminator(label);
        if(disc < current_best_disc){
            best = &hyp;
            current_best_disc = disc;
        }
    }
    if(std::isfinite(current_best_disc)){
      //cout << "minimal Chi2: " << current_best_disc << endl;
        return best;
    }
    else{
        return nullptr;
    }
}





TTbarFullhadRecoChi2Discriminator::TTbarFullhadRecoChi2Discriminator(Context & ctx, const std::string & rechyps_name, const cfg & config_): config(config_){
    h_hyps = ctx.get_handle<vector<TTbarFullhadRecoHypothesis>>(rechyps_name);
}


bool TTbarFullhadRecoChi2Discriminator::process(uhh2::Event & event){
    auto & hyps = event.get(h_hyps);
    const double mass_thad1 = 181;
    const double mass_thad1_sigma = 15;

    for(auto & hyp: hyps){
        double mass_thad1_rec = inv_mass(hyp.tophad1_v4());
        //double mass_thad2_rec = inv_mass(hyp.tophad2_v4());

	double chi2 = pow((mass_thad1_rec-mass_thad1) / mass_thad1_sigma,2)/* + pow((mass_thad2_rec-mass_thad1) / mass_thad1_sigma,2)*/;

        hyp.set_discriminator(config.discriminator_label, chi2); // modified

        //hyp.set_discriminator(config.discriminator_label + "_whad1", chi2_whad1);// added
        //hyp.set_discriminator(config.discriminator_label + "_whad2", chi2_whad2);// added
  }
  return true;
}
