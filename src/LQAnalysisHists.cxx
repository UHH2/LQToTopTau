#include "UHH2/LQAnalysis/include/LQAnalysisHists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/common/include/TTbarGen.h"

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

LQAnalysisHists::LQAnalysisHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  Double_t bins[9] = {0, 400, 700, 1000, 1300, 1600, 1900, 2200, 2500};
  //double newbin[3] = {700,1500,4000};
  // ht
  book<TH1F>("MET", "missing E_{T}", 20,0,1000);
  book<TH1F>("MET_binned", "missing E_{T}", 50,0,1000);
  book<TH1F>("HT", "H_{T} Jets", 50, 0, 3500);
  book<TH1F>("HTLep", "H_{T} Lep", 50, 0, 1000);
  book<TH1F>("ST", "H_{T}", 50, 0, 5000);
  book<TH1F>("ST_binned", "H_{T}", 8, bins);
  book<TH1F>("ST_testbinned", "H_{T}", 32, 800,4000);
  //book<TH1F>("pt_gentau_hist", "p_{T}^{gen}", 100,0,800);
  //book<TH1F>("HT_weighttest", "H_{T} Jets", 140, 0, 3500);
  //book<TH1F>("ptj1", "p_{T} first jet", 15, 0,1200);



}


void LQAnalysisHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
 
  hist("MET")->Fill(event.met->pt(), weight);
  hist("MET_binned")->Fill(event.met->pt(), weight);

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
  hist("ST_testbinned")->Fill(ht+ht_lep+met, weight);

  


  /*
  Event::Handle<TTbarGen> h_ttbargen;
  std::unique_ptr<AnalysisModule> ttgenprod;
  ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
  ttgenprod->process(event);
  const auto & ttbargen = event.get(h_ttbargen);
  */
  

  /*
  Event::Handle<TTbarGen> h_ttbargen;
  std::unique_ptr<AnalysisModule> ttgenprod;

  ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  //h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
  */

  /*
  const vector<GenParticle> * genparticles = event.genparticles;

  int invalid_daughter = (unsigned short)(-1);
  double result = 0.0;
  for(const auto & gp : *genparticles){
    if(gp.daughter1() != invalid_daughter || gp.daughter2() != invalid_daughter) continue;
    int id = abs(gp.pdgId());
    if((id >= 1 && id <= 5) || (id == 21)){
      result += gp.pt();
    }
  }

  hist("HT_weighttest")->Fill(result, weight);
  */
  

  /*
  const vector<GenParticle> * genparticles = event.genparticles;
  const vector<Tau> * taus = event.taus;
 
  for(auto & genp : *genparticles){ 
    if (abs(genp.pdgId())!=15) continue;
    for(auto & tau : *taus){
      if(deltaR(genp,tau)<0.4){
	hist("pt_gentau_hist")->Fill(genp.pt(),weight);
      }
    }
  }
  */
 
  

}

LQAnalysisHists::~LQAnalysisHists(){}
