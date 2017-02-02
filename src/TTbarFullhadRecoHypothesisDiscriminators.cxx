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

	//double chi2 = pow((mass_thad1_rec-mass_thad1) / mass_thad1_sigma,2);// + pow((mass_thad2_rec-mass_thad1) / mass_thad1_sigma,2);
	double chi2 = hyp.chi2();

        hyp.set_discriminator(config.discriminator_label, chi2); // modified

        //hyp.set_discriminator(config.discriminator_label + "_whad1", chi2_whad1);// added
        //hyp.set_discriminator(config.discriminator_label + "_whad2", chi2_whad2);// added
  }
  return true;
}



LQCorrectMatchDiscriminator::LQCorrectMatchDiscriminator(Context & ctx, const std::string & rechyps_name, const cfg & config_): config(config_){
  h_hyps = ctx.get_handle<vector<TTbarFullhadRecoHypothesis>>(rechyps_name);
  h_ttbargen = ctx.get_handle<TTbarGen>(config.ttbargen_name);
  //h_LQLQbargen = ctx.get_handle<LQGen>(config.LQLQbargen_name);
}

namespace {
  // match particle p to one of the jets (Delta R < 0.3); return the deltaR
  // of the match.
  template<typename T> // T should inherit from Particle
  float match_dr(const Particle & p, const std::vector<T> & jets, int& index){
    float mindr = infinity;
    index = -1;
    for(unsigned int i=0; i<jets.size(); ++i){
      float dR = deltaR(p, jets.at(i));
      if( dR <0.4 && dR<mindr) {
        mindr=dR;
        index=i;
	//cout << "dR: " << dR << endl;
      }
    }
    return mindr;
  }
}

bool LQCorrectMatchDiscriminator::process(uhh2::Event & event){  // replaced 'infinity' with '999999'
  auto & hyps = event.get(h_hyps);
  const auto & ttbargen = event.get(h_ttbargen);
  //const auto & LQLQbargen = event.get(h_LQLQbargen);
  auto dec = ttbargen.DecayChannel();
  for(auto & hyp: hyps){
    hyp.set_discriminator(config.discriminator_label, 999999);
  }
  if(!(dec==TTbarGen::e_muhad || dec==TTbarGen::e_ehad || dec==TTbarGen::e_tauhad || dec==TTbarGen::e_had)){ // if not semileptonic channel
    for(auto & hyp: hyps){
      hyp.set_discriminator(config.discriminator_label, 999999);
      LorentzVector dummy = {0,0,0,0};
      hyp.set_WJet1_v4(dummy);
      hyp.set_WJet2_v4(dummy);
    }
    return true;
  }

  // note that it is allowed that two partons from the hadronic ttbar decay match the same jet.
  for(auto & hyp: hyps){
    LorentzVector dummy = {0,0,0,0};
    hyp.set_WJet1_v4(dummy);
    hyp.set_WJet2_v4(dummy);

    auto hadr1_jets = hyp.tophad1_jets();
    //auto hadr2_jets = hyp.tophad2_jets();

    if(hadr1_jets.size() > 3){ // < 3 is allowed ...
      hyp.set_discriminator(config.discriminator_label, 999999);
      continue;
    }

    //index lists of jets that can be matched to partons
    std::set<int> matched_hadr_jets, matched_hadr_jets1, matched_hadr_jets2;
    float temp_dr1, temp_dr2;

    // match b jets
    int index_h, index_h1, index_h2;
    float correct_dr=0;
    if(dec==TTbarGen::e_muhad || dec==TTbarGen::e_ehad || dec==TTbarGen::e_tauhad){
      correct_dr = match_dr(ttbargen.BHad(), hadr1_jets, index_h);
      if(index_h >= 0) matched_hadr_jets.insert(index_h);
      //match quarks from W decays
      correct_dr += match_dr(ttbargen.Q1(), hadr1_jets, index_h);
      if(index_h >= 0){
	matched_hadr_jets.insert(index_h);
	hyp.set_WJet1_v4(hadr1_jets[index_h].v4());
      }
      correct_dr += match_dr(ttbargen.Q2(), hadr1_jets, index_h);
      if(index_h >= 0) {
	matched_hadr_jets.insert(index_h);
	hyp.set_WJet2_v4(hadr1_jets[index_h].v4());
      }

      if(matched_hadr_jets.size() != hadr1_jets.size()){
	hyp.set_discriminator(config.discriminator_label, 999999);
	continue;
      }
    }
    if(dec==TTbarGen::e_had){
      //cout << endl << endl << "Antitop test" << endl;
      //cout << "antib quark" << endl;
      temp_dr1 = match_dr(ttbargen.bAntitop(), hadr1_jets, index_h1);
      if(index_h1 >= 0) matched_hadr_jets1.insert(index_h1);
      //match quarks from W decays
      //cout << "antiq1" << endl;
      temp_dr1 += match_dr(ttbargen.WMinusdecay1(), hadr1_jets, index_h1);
      if(index_h1 >= 0){
	matched_hadr_jets1.insert(index_h1);
	hyp.set_WJet1_v4(hadr1_jets[index_h1].v4());
      }
      //cout << "antiq2" << endl;
      temp_dr1 += match_dr(ttbargen.WMinusdecay2(), hadr1_jets, index_h1);
      if(index_h1 >= 0){
	matched_hadr_jets1.insert(index_h1);
	hyp.set_WJet2_v4(hadr1_jets[index_h1].v4());
      }
      //cout << "Top test" << endl;
      //cout << "b quark" << endl;
      temp_dr2 = match_dr(ttbargen.bTop(), hadr1_jets, index_h2);
      if(index_h2 >= 0) matched_hadr_jets2.insert(index_h2);
      //match quarks from W decays
      //cout << "q1" << endl;
      temp_dr2 += match_dr(ttbargen.Wdecay1(), hadr1_jets, index_h2);
      if(index_h2 >= 0){
	matched_hadr_jets2.insert(index_h2);
	hyp.set_WJet1_v4(hadr1_jets[index_h2].v4());
      }
      //cout << "q2" << endl;
      temp_dr2 += match_dr(ttbargen.Wdecay2(), hadr1_jets, index_h2);
      if(index_h2 >= 0){
	matched_hadr_jets2.insert(index_h2);
	hyp.set_WJet2_v4(hadr1_jets[index_h2].v4());
      }
      if(temp_dr1<temp_dr2){
	correct_dr = temp_dr1;
	if(matched_hadr_jets1.size() != hadr1_jets.size()){
	  hyp.set_discriminator(config.discriminator_label, 999999);
	  continue;
	}
      }
      else{
	correct_dr = temp_dr2;
	if(matched_hadr_jets2.size() != hadr1_jets.size()){
	  hyp.set_discriminator(config.discriminator_label, 999999);
	  continue;
	}
      }

    }



    // if not all jets of the hadronic side of the reconstruction could be matched: infinite
    // value:
    /*
    if(matched_hadr_jets.size() != hadr1_jets.size()){
      hyp.set_discriminator(config.discriminator_label, 999999);
      continue;
    }
    */



    //set final dr as discriminator value
    hyp.set_discriminator(config.discriminator_label, correct_dr);
    


  }
  return true;
}


/*
LQCorrectMatchDiscriminator::LQCorrectMatchDiscriminator(Context & ctx, const std::string & rechyps_name, const cfg & config_): config(config_){
  h_hyps = ctx.get_handle<vector<TTbarFullhadRecoHypothesis>>(rechyps_name);
  h_ttbargen = ctx.get_handle<TTbarGen>(config.ttbargen_name);
  //h_LQLQbargen = ctx.get_handle<LQGen>(config.LQLQbargen_name);
}

namespace {
  // match particle p to one of the jets (Delta R < 0.3); return the deltaR
  // of the match.
  template<typename T> // T should inherit from Particle
  float match_dr(const Particle & p, const std::vector<T> & jets, int& index){
    float mindr = infinity;
    index = -1;
    for(unsigned int i=0; i<jets.size(); ++i){
      float dR = deltaR(p, jets.at(i));
      if( dR <0.3 && dR<mindr) {
        mindr=dR;
        index=i;
      }
    }
    return mindr;
  }
}


bool LQCorrectMatchDiscriminator::process(uhh2::Event & event){  // replaced 'infinity' with '999999'
  auto & hyps = event.get(h_hyps);

  const auto & ttbargen = event.get(h_ttbargen);

  //const auto & LQLQbargen = event.get(h_LQLQbargen);
  auto dec = ttbargen.DecayChannel();
  if(dec != TTbarGen::e_muhad){ // if not semilep mu channel
    for(auto & hyp: hyps){
      hyp.set_discriminator(config.discriminator_label, 999999);
    }
    return true;
  }

  // note that it is allowed that two partons from the hadronic ttbar decay match the same jet.
  for(auto & hyp: hyps){
    auto hadr1_jets = hyp.tophad1_jets();
    //auto hadr2_jets = hyp.tophad2_jets();

    if(hadr1_jets.size() > 3){ // < 3 is allowed ...
      hyp.set_discriminator(config.discriminator_label, 999999);
      continue;
    }

    //index lists of jets that can be matched to partons
    std::set<int> matched_hadr_jets;

    // match b jets
    int index_h;
    
    float correct_dr = match_dr(ttbargen.BHad(), hadr1_jets, index_h);
    if(index_h >= 0) matched_hadr_jets.insert(index_h);
    //match quarks from W decays
    correct_dr += match_dr(ttbargen.Q1(), hadr1_jets, index_h);
    if(index_h >= 0) matched_hadr_jets.insert(index_h);
    correct_dr += match_dr(ttbargen.Q2(), hadr1_jets, index_h);
    if(index_h >= 0) matched_hadr_jets.insert(index_h);

    // if not all jets of the hadronic side of the reconstruction could be matched: infinite
    // value:
    if(matched_hadr_jets.size() != hadr1_jets.size()){
      hyp.set_discriminator(config.discriminator_label, 999999);
      continue;
    }
    //set final dr as discriminator value
    hyp.set_discriminator(config.discriminator_label, correct_dr);
    


  }

  return true;
}
*/
