#include "UHH2/LQAnalysis/include/TTbarFullhadReco.h"
#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"

#include <cassert>

using namespace uhh2;
using namespace std;


HighMassHadronicTTbarReco::HighMassHadronicTTbarReco(Context & ctx, const string & label) {
  h_hadr_recohyps = ctx.declare_event_output<vector<TTbarFullhadRecoHypothesis>>(label);
}

HighMassHadronicTTbarReco::~HighMassHadronicTTbarReco() {}

bool HighMassHadronicTTbarReco::process(uhh2::Event & event) {
  assert(event.jets);

  const double mass_thad1 = 181;
  const double mass_thad1_sigma = 15;

  std::vector<TTbarFullhadRecoHypothesis> recoHyps;

  unsigned int n_jets = event.jets->size();
  if(n_jets>7) n_jets=7; //avoid crashes in events with many jets
  // idea: loop over 3^Njet possibilities and write the current loop
  // index j in the 3-base system. The Njets digits represent whether
  // to assign each jet to the hadronic side (0), leptonic side (1),
  // or none of them (2).
  const unsigned int max_j = pow(3, n_jets);

  double chi2 = 99999999999;
  double chi2min = chi2;
  LorentzVector LQ11_v4;
  LorentzVector LQ12_v4;
  LorentzVector LQ21_v4;
  LorentzVector LQ22_v4;

  //loop over neutrino solutions and jet assignments to fill hyotheses
  for (unsigned int j=0; j < max_j; j++) {
    LorentzVector tophad1_v4;
    LorentzVector tophad2_v4;
    int hadjets1=0;
    int hadjets2=0;
    int num = j;
    TTbarFullhadRecoHypothesis hyp;
    for (unsigned int k=0; k<n_jets; k++) {
      if(num%3==0) {
	tophad1_v4 = tophad1_v4 + event.jets->at(k).v4();
	hyp.add_tophad1_jet(event.jets->at(k));
	hadjets1++;
      }

      if(num%3==1) {
	tophad2_v4 = tophad2_v4 + event.jets->at(k).v4();
	hyp.add_tophad2_jet(event.jets->at(k));
	hadjets2++;
      }
      //in case num%3==2 do not take this jet at all
      //shift the trigits of num to the right:
      num /= 3;
    }
   
    double mass_thad1_rec = 0;
    double mass_thad2_rec = 0;

    if(tophad1_v4.isTimelike())      mass_thad1_rec = tophad1_v4.M();
    else      mass_thad1_rec = -sqrt(-tophad1_v4.mass2());
   
    if(tophad2_v4.isTimelike())      mass_thad2_rec = tophad2_v4.M();
    else      mass_thad2_rec = -sqrt(-tophad2_v4.mass2());

    const auto & tau1 = (*event.taus)[0];
    const auto & tau2 = (*event.taus)[1];

    LQ11_v4 = tau1.v4()+tophad1_v4;
    LQ12_v4 = tau1.v4()+tophad2_v4;
    LQ21_v4 = tau2.v4()+tophad1_v4;
    LQ22_v4 = tau2.v4()+tophad2_v4;


    double Mdiff1 = fabs(LQ11_v4.M() - LQ22_v4.M());
    double Mdiff2 = fabs(LQ12_v4.M() - LQ21_v4.M());
    if(Mdiff1<Mdiff2){
      hyp.set_LQ1_v4(LQ11_v4);
      hyp.set_LQ2_v4(LQ22_v4);
    }
    else{
      hyp.set_LQ1_v4(LQ12_v4);
      hyp.set_LQ2_v4(LQ21_v4);
    }

    if(hadjets1>0 && hadjets2>0) {
      chi2 = pow((mass_thad1_rec-mass_thad1) / mass_thad1_sigma,2) + pow((mass_thad2_rec-mass_thad1) / mass_thad1_sigma,2);
      if(chi2 < chi2min){
	chi2min = chi2;
	hyp.set_tophad1_v4(tophad1_v4);
	hyp.set_tophad2_v4(tophad2_v4);

	//cout << mass_thad1_rec << endl;

	if(recoHyps.size() == 1){
	  //cout << "recoHyps-size = 1, popping back the element" << endl;
	  recoHyps.pop_back();
	  //cout << "emplacing back new best hypothesis" << endl;
	  recoHyps.push_back(hyp);
	}
	else if(recoHyps.size() == 0){
	  //cout << "recoHyps empty, emplacing back best (=first) hypothesis" << endl;
	  recoHyps.push_back(hyp);
	}
	else throw runtime_error("size of recoHyps neither 0 nor 1");
    
      }
    }
 
  }
  //cout << endl << endl;

  //cout << recoHyps.at(0).tophad1_v4().M() << endl;

  event.set(h_hadr_recohyps, move(recoHyps));

  return true;

}

