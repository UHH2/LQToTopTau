#include "UHH2/LQAnalysis/include/LQAnalysisHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

LQAnalysisHists::LQAnalysisHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  Double_t bins[7] = {1200, 1500, 1800, 2000, 2500, 3000, 4000};
  // ht
  book<TH1F>("MET", "missing E_{T}", 20,0,1000);
  book<TH1F>("HT", "H_{T} Jets", 50, 0, 3500);
  book<TH1F>("HTLep", "H_{T} Lep", 50, 0, 1000);
  book<TH1F>("ST", "S_{T}", 38, 1200, 5000);
  book<TH1F>("ST_binned", "S_{T}", 6, bins);
}


void LQAnalysisHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
 
  hist("MET")->Fill(event.met->pt(), weight);

  auto met = event.met->pt();
  double ht = 0.0;
  for(const auto & jet : *event.jets){
    ht += jet.pt();
  }
  hist("HT")->Fill(ht, weight);

  double ht_lep = 0.0;

  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }
  for(const auto & tau : *event.taus){
    ht_lep += tau.pt();
  }

  hist("HTLep")->Fill(ht_lep, weight);

  hist("ST")->Fill(ht+ht_lep+met, weight);
  hist("ST_binned")->Fill(ht+ht_lep+met, weight);


}

LQAnalysisHists::~LQAnalysisHists(){}
