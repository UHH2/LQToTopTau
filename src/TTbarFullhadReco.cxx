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

  const double mass_thad1 = 181;//181
  const double mass_thad1_sigma = 15;//15


  std::vector<TTbarFullhadRecoHypothesis> recoHyps;

  unsigned int n_jets = event.jets->size();
  if(n_jets>6) n_jets=6; //avoid crashes in events with many jets
  // idea: loop over 3^Njet possibilities and write the current loop
  // index j in the 3-base system. The Njets digits represent whether
  // to assign each jet to the hadronic side (0), leptonic side (1),
  // or none of them (2).
  const unsigned int max_j = pow(3, n_jets);

  double chi2 = 99999999999;
  double tchi2 = 99999999999;
  double wchi2 = 99999999999;

  double chi2min = chi2;
  LorentzVector LQ11_v4;
  LorentzVector LQ12_v4;
  LorentzVector LQ21_v4;
  LorentzVector LQ22_v4;

  //loop over neutrino solutions and jet assignments to fill hyotheses
  for (unsigned int j=0; j < max_j; j++) {
    LorentzVector tophad1_v4;
    int hadjets1=0;
    int num = j;
    TTbarFullhadRecoHypothesis hyp;
    for (unsigned int k=0; k<n_jets; k++) {
      if(num%3==0) {
	tophad1_v4 = tophad1_v4 + event.jets->at(k).v4();
	hyp.add_tophad1_jet(event.jets->at(k));
	hadjets1++;
      }
      //in case num%3==2 do not take this jet at all
      //shift the trigits of num to the right:
      num /= 3;
    }

    
    std::vector<LorentzVector> whyps_v4;
    for (unsigned int i=0; i<pow(2,hyp.tophad1_jets().size()); i++) {
      LorentzVector whad1_v4;
      unsigned int whadjets1=0;
      int wnum = i;
      for (unsigned int k=0; k<hyp.tophad1_jets().size(); k++) {
	if(wnum%2==0) {
	  whad1_v4 = whad1_v4 + event.jets->at(k).v4();
	  whadjets1++;
	}
	wnum /= 2;
      }
      if(whadjets1>0 && whadjets1<hyp.tophad1_jets().size() && hadjets1>1){
	whyps_v4.push_back(whad1_v4);
      }
    }
    
    double mass_whad1_rec = 84;
    double mass_whad1 = 84;
    double mass_whad1_sigma = 10.6;
    double wchi2min = 99999999;
    int best_ind = -1;
    
    if(whyps_v4.size()>0){
      for(unsigned int i=0; i<whyps_v4.size(); i++){
	if(whyps_v4[i].isTimelike())      mass_whad1_rec = whyps_v4[i].M();
	else      mass_whad1_rec = -sqrt(-whyps_v4[i].mass2());
	double chi2w = pow((mass_whad1_rec - mass_whad1) / mass_whad1_sigma,2);
	if(chi2w<wchi2min){
	  wchi2min = chi2w;
	  best_ind = i;
	}
      }
    }



    
    double mass_thad1_rec = 0;
 

    if(tophad1_v4.isTimelike())      mass_thad1_rec = tophad1_v4.M();
    else      mass_thad1_rec = -sqrt(-tophad1_v4.mass2());
    

    double pt_thad1_rec = tophad1_v4.Pt();
    double maxDRfit = (-4*pow(10,-9))*pow(pt_thad1_rec,3) + (1.16*pow(10,-5))*pow(pt_thad1_rec,2) + (-0.01)*pt_thad1_rec+3.66;

    int nsubjets = hyp.tophad1_jets().size();
    double maxDR=0.;
    for(int i=0; i<nsubjets; i++){
      for(int j=0; j<nsubjets; j++){
	if(i!=j){
	  double dr_dummy = deltaR(hyp.tophad1_jets().at(i),hyp.tophad1_jets().at(j));
	  if(dr_dummy>maxDR){
	    maxDR = dr_dummy;
	  }
	}
      }
    }
    

    if(hadjets1>0 /*&& hadjets2>0*/) {
      /*    
      if(whyps_v4[best_ind].isTimelike())   mass_whad1_rec = whyps_v4[best_ind].M();
      else      mass_whad1_rec = -sqrt(-whyps_v4[best_ind].mass2());
      tchi2 = pow((mass_thad1_rec-mass_thad1) / mass_thad1_sigma,2);
      wchi2 = pow((mass_whad1_rec-mass_whad1) / mass_whad1_sigma,2);
      */
      //cout << "wmass before: " << mass_whad1_rec << endl;

      chi2 = pow((mass_thad1_rec-mass_thad1) / mass_thad1_sigma,2);// + pow((mass_whad1_rec-mass_whad1) / mass_whad1_sigma,2) + pow((maxDR-maxDRfit+0.29)/0.26,2);
      
      if(chi2 < chi2min){
	chi2min = chi2;
	hyp.set_tophad1_v4(tophad1_v4);
	//hyp.set_whad1_v4(whyps_v4[best_ind]);
	hyp.set_chi2(chi2);
	//cout << mass_thad1_rec << endl;

	hyp.set_tchi2(tchi2);
	hyp.set_wchi2(wchi2);
      

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
	
	
	//recoHyps.emplace_back(move(hyp)); // asdf
	//cout << "tmass: " << hyp.tophad1_v4().M() << endl;
	//cout << "wmass v4:  " << hyp.whad1_v4().M() << endl;
	//cout << "wmass after:  " << mass_whad1_rec << endl;
      }
      

      //get all hypotheses
      /*
	hyp.set_tophad1_v4(tophad1_v4);
	hyp.set_whad1_v4(whyps_v4[best_ind]);
	hyp.set_chi2(chi2);
	//cout << mass_thad1_rec << endl;
      
	hyp.set_tchi2(tchi2);
	hyp.set_wchi2(wchi2);
      
	recoHyps.push_back(hyp);
      */
    }

 
  }

  //cout << "+++++++++++++++++++++++++++++" << endl;

  //cout << endl << endl;

  //cout << recoHyps.at(0).tophad1_v4().M() << endl;

  event.set(h_hadr_recohyps, move(recoHyps));

  return true;

}







/*
AllHadronicTTbarReco::AllHadronicTTbarReco(Context & ctx, const string & label) {
  h_hadr_recohyps = ctx.declare_event_output<vector<TTbarFullhadRecoHypothesis>>(label);
}

AllHadronicTTbarReco::~AllHadronicTTbarReco() {}

bool AllHadronicTTbarReco::process(uhh2::Event & event) {
  assert(event.jets);

  const double mass_thad1 = 181;
  const double mass_thad1_sigma = 15;

  std::vector<TTbarFullhadRecoHypothesis> recoHyps;

  unsigned int n_jets = event.jets->size();
  if(n_jets>6) n_jets=6; 
  const unsigned int max_j = pow(3, n_jets);

  double chi2 = 99999999999;
  double chi2min = chi2;
  LorentzVector LQ11_v4;
  LorentzVector LQ12_v4;
  LorentzVector LQ21_v4;
  LorentzVector LQ22_v4;

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
      num /= 3;
    }
   
    double mass_thad1_rec = 0;
    double mass_thad2_rec = 0;

    if(tophad1_v4.isTimelike())      mass_thad1_rec = tophad1_v4.M();
    else      mass_thad1_rec = -sqrt(-tophad1_v4.mass2());
    
    if(tophad2_v4.isTimelike())      mass_thad2_rec = tophad2_v4.M();
    else      mass_thad2_rec = -sqrt(-tophad2_v4.mass2());
    

    if(hadjets1>0 && hadjets2>0) {
      chi2 = pow((mass_thad1_rec-mass_thad1) / mass_thad1_sigma,2) + pow((mass_thad2_rec-mass_thad1) / mass_thad1_sigma,2);
      if(chi2 < chi2min){
	chi2min = chi2;
	hyp.set_tophad1_v4(tophad1_v4);
	hyp.set_tophad2_v4(tophad2_v4);



	if(recoHyps.size() == 1){

	  recoHyps.pop_back();

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
*/
