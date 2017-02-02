#include "UHH2/LQAnalysis/include/LQAnalysisHists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"


#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
//#include "TMath.h"

using namespace std;
using namespace uhh2;

LQAnalysisHists::LQAnalysisHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  Double_t bins[9] = {0, 400, 700, 1000, 1300, 1600, 1900, 2200, 2500};
  Double_t taubins[5] = {20, 60, 120, 200, 800};
  Double_t bbins[5] = {20, 100, 200, 400, 1000};
  //double newbin[3] = {700,1500,4000};
  // ht
  book<TH1F>("MET", "missing E_{T}", 20,0,1000);
  book<TH1F>("MET_binned", "missing E_{T}", 50,0,1000);
  book<TH1F>("HT", "H_{T} Jets", 50, 0, 3500);
  book<TH1F>("HTLep", "H_{T} Lep", 50, 0, 1000);
  book<TH1F>("ST", "S_{T}", 50, 0, 5000);
  book<TH1F>("ST_binned", "S_{T}", 8, bins);
  book<TH1F>("ST_testbinned", "S_{T}", 20, 0,3500);

  book<TH1F>("pt_firstjet", "p_{T} first jet [GeV]", 20, 0, 1500);

  book<TH1F>("isolation","#tau iso [GeV]",30,0,3);
  book<TH1F>("isolation_sideband","#tau iso [GeV]",15,0,300);
  /*
  double isobins[101];
  isobins[0]=0.5;
  isobins[1]=1.5;
  for(int i=2; i<101; i++){
    isobins[i]=i*4;
  }
  */
  double isobins[15] = {0.5,1.5,3,6,10,15,24,40,65,100,150,200,300,400,600};
  book<TH1F>("isolation_fullrange","#tau iso [GeV]",14,isobins);

  book<TH1F>("MVA_isolation","MVA based #tau iso [GeV]",20,0,2);

  book<TH1F>("isolated_taus","#tau isolation [GeV]",30,0,3);
  book<TH1F>("nonisolated_taus","#tau nonisolated [GeV]",14,isobins);
  book<TH1F>("N_isolated","number of isolated #tau",6,-0.5,5.5);
  book<TH1F>("N_nonisolated","number of nonisolated #tau",6,-0.5,5.5);

  /*
  double pttoprebins[8]={0,80,130,180,230,300,400,1000};
  double pttopbins[5]={0,70,140,200,1000};
  book<TH1F>("M_tophad_own", "M^{top,had} [GeV/c^{2}]", 50, 0, 500 ) ;
  book<TH1F>("Pt_tophad_own", "P_{T}^{top,had} [GeV/c]", 60, 0, 1200 ) ;
  book<TH1F>("Pt_tophad_own_rebin", "P_{T}^{top,had} [GeV/c]", 7, pttoprebins ) ;
  book<TH1F>("Pt_tophad_own_binned", "P_{T}^{top,had} [GeV/c]", 4, pttopbins ) ;

  book<TH1F>("M_tophad_tau1", "M^{top,tau1} [GeV/c^{2}]", 50, 0, 1000 ) ; 
  book<TH1F>("DeltaR_tophad_tau1", "DeltaR_tophad_tau1", 50, 0, 5 ) ; 
  */

  //book<TH1F>("pt_gentau_hist", "p_{T}^{gen}", 100,0,800);
  //book<TH1F>("HT_weighttest", "H_{T} Jets", 140, 0, 3500);
  //book<TH1F>("ptj1", "p_{T} first jet", 15, 0,1200);
  book<TH1F>("M_mumu", "M_{#mu#mu}", 100, 0,1000);
  book<TH1F>("M_mutau", "M_{#mu#tau}", 100, 0,1000);
  book<TH1F>("NbJetsL", "Number of bjets loose", 7, -0.5,6.5);
  book<TH1F>("NbJetsM", "Number of bjets medium", 7, -0.5,6.5);
  book<TH1F>("NbJetsT", "Number of bjets tight", 7, -0.5,6.5);

  book<TH1F>("h_minMT", "M_{T}{b, missing ET}", 16, 0 ,800);
  book<TH1F>("h_MT_btau", "M_{T}{b, tau}", 16, 0 ,800);
  book<TH1F>("h_MT_METtau", "M_{T}{tau, missing ET}", 16, 0 ,800);
  book<TH1F>("M_btau", "M_{b#tau}", 20, 0,1000);
  book<TH1F>("M_tautau", "M_{#tau#tau}", 25, 0,1000);

  book<TH1F>("DeltaR_btau", "DeltaR_{b#tau}", 50, 0,5);
  book<TH1F>("pt_b", "pt bjetM", 50,0,1000);
  book<TH1F>("pt_b_binned", "pt bjetM", 4,bbins);

  book<TH1F>("DeltaR_mutau", "deltaR(#mu1,#tau1)", 50, 0,5);

  book<TH1F>("MT", "M_{T}(mu,missing ET)", 16,0,800);
  book<TH1F>("deltaPhi_mu_met", "#Delta#varphi (#mu, MET)", 25,0,5);
  //book<TH1F>("Weights", "weights", 2,0.5,2.5);
  book<TH1F>("eta_tilde", "|#tilde{#eta}|", 8,0,2.4);
  book<TH1F>("eta_tilde2", "|#tilde{#eta}|", 24,0,2.4);
  double eta_bins[3]={0,0.9,2.5};
  book<TH1F>("eta_tilde_bin1", "|#tilde{#eta}|", 2,eta_bins);
  book<TH1F>("eta_tilde_bin2", "|#tilde{#eta}|", 2,eta_bins);
  book<TH1F>("muon_type", "0 real muon, 1 fake muon", 2,-0.5,1.5);
  book<TH1F>("electron_type", "0 real electron, 1 fake electrron", 2,-0.5,1.5);
  book<TH1F>("tau_type", "0 real tau, 1 fake tau", 2,-0.5,1.5);
  book<TH1F>("tau_categories", "0: fakes, 1: reals, 2: both, 3: others", 4,-0.5,3.5);
  book<TH1F>("N_faketaus", "number of fake taus", 6,-0.5,5.5);
  book<TH1F>("N_realtaus", "number of real taus", 6,-0.5,5.5);
  book<TH1F>("tausign", "DiTau, 0: SS, 1: OS", 2,-0.5,1.5);
  book<TH1F>("N_OStau", "number of OS tau pairs", 6,-0.5,5.5);
  book<TH1F>("N_SStau", "number of SS tau pairs", 6,-0.5,5.5);

  book<TH1F>("pt_real_tau1_binned","p_{T} real tau 1",4,taubins);
  book<TH1F>("pt_fake_tau1_binned","p_{T} fake tau 1",4,taubins);
  book<TH1F>("pt_tau1_binned","p_{T} tau 1",4,taubins);
  book<TH1F>("pt_tau_binned","p_{T} taus",4,taubins);

  book<TH1F>("pt_tau1", "p_{T} first tau [GeV]", 20, 20, 620);

  book<TH1F>("M_jet", "M_{Jet}", 100, 0, 2000);
  book<TH1F>("N_subjets", "N_{Subjets} in a Topjet", 11, -0.5, 10.5);
  book<TH1F>("min_mDisubjet", "Min(m_{ij})", 50, 0, 1000);
  book<TH1F>("N_TopTags", "Number of CMSTopTags",6 ,-0.5, 5.5 );

  book<TH1F>("pt_thad", "P_{T}^{top,had} [GeV]", 20,0,1200);
  book<TH1F>("pt_tlep", "P_{T}^{top,lep} [GeV]", 20,0,1200);
  book<TH1F>("pt_tcom", "P_{T}^{top,combined} [GeV]", 20,0,1200);

  book<TH1F>("faketaupt_realtaupt", "0: p_{T}^{f #tau 1} > p_{T}^{r #tau 1}, 1: p_{T}^{r #tau 1} > p_{T}^{f #tau 1}", 2,-0.5,1.5);
  book<TH1F>("sign", "#mu#tau sign", 2,-1,1);

  book<TH1F>("N", "counting exp.", 1,0.5,1.5);


  double pttoprebins[5]={0,100,200,300,1200};
  double pttopbins[4]={0,100,200,1200};
  double testbin[5] = {0,120,200,380,1200};

  book<TH1F>("M_tophad_own", "M^{top,had} [GeV]", 40, 0, 400 ) ;
  book<TH1F>("Pt_tophad_own", "P_{T}^{top,had} [GeV]", 60, 0, 1200 ) ;
  book<TH1F>("Pt_tophad_own_rebin", "P_{T}^{top,had} [GeV]", 12, 0, 1200 ) ;
  book<TH1F>("Pt_tophad_own_binned", "P_{T}^{top,had} [GeV]", 4, pttoprebins ) ;
  book<TH1F>("Pt_tophad_own_binned2", "P_{T}^{top,had} [GeV]", 3, pttopbins ) ;
  book<TH1F>("Pt_tophad_own_test", "P_{T}^{top,had} [GeV]", 4, testbin ) ;
  book <TH1F>("Chi2", "#chi^{2}", 100, 0, 200);
  book<TH1F>("deltaR_top_gentop", "DeltaR(TopHyp, GenTop)", 50, 0, 5);

  book <TH1F>("top_recjets", "number of jets used for top hypothesis", 6, 0.5, 6.5);
  book<TH1F>("deltaR_recjets", "largest dR of jets from best hypothesis", 50, 0, 5);
  book<TH2F>("Pt_tophad_vs_largestDR", "p_{T}^{top,had} vs largest DR of jets in hypothesis", 60,0,1200, 50,0,5);

  book<TH1F>("M_tophad_tau", "M_{top#tau}", 20, 0,1000);
  book<TH1F>("deltaPhi_topjet_jet", "dPhi(topjet, jet)", 100, 0, 3.2);
  book<TH1F>("M_tophad_jet", "M_{top,jet}", 20, 0,1000);


  book<TH1F>("deltaR_genjets", "largest dR of top partons", 50, 0, 5);
  book<TH2F>("Pt_tophad_vs_largestDR_gen", "p_{T}^{top,had} vs largest DR of partons", 60,0,1200, 50,0,5);
  
  book<TH1F>("pt_top_resolution", "(p_{T}^{top,rec} - p_{T}^{top,gen}) / p_{T}^{top,gen}", 50, -1,1);
  //book<TH1F>("Pt_top_gen", "p_{T}^{top,gen}", 24, 0,1200);

  book<TH1F>("M_whad_own", "M^{W,had} [GeV]", 100, 0, 200 ) ;
  book<TH1F>("WMass_matched", "M^{W,had} matched [GeV]", 200, 0, 200 ) ;
  book<TH1F>("maxDR_minus_maxDRfit", "maxDR - maxDRfit", 200, -5, 5);

  book <TH1F>("WChi2", "W #chi^{2}", 100, 0, 200);
  book <TH1F>("tChi2", "t #chi^{2}", 100, 0, 200);

  
  book<TH2F>("pt_tau1_vs_ST","pt tau 1 vs ST", 8, 0, 800 , 15, 0, 3000);
  book<TH2F>("pt_top_vs_ST","pt top vs ST", 12, 0, 1200 ,15, 0, 3000);
  book<TH2F>("pt_top_vs_pt_tau1","pt top vs pt tau 1", 12, 0, 1200 ,80, 0, 800);
  

  //h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>("HighMassTTbarReconstruction");
  //m_discriminator_name = "Chi2";


  h_hadr_hyps = ctx.get_handle<std::vector<TTbarFullhadRecoHypothesis>>("HighMassHadronicTTbarFullhadReco");

  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  auto dataset_type = ctx.get("dataset_type");
  is_data = dataset_type == "DATA";

}


void LQAnalysisHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  
  std::vector<TTbarFullhadRecoHypothesis> hadr_hyps = event.get(h_hadr_hyps);
  const TTbarFullhadRecoHypothesis* hadr_hyp = get_best_hypothesis( hadr_hyps, "Chi2Hadronic" );

 
  //const auto & ttbargen = event.get(h_ttbargen);



  double mTopHad = hadr_hyp->tophad1_v4().M();
  double ptTopHad = hadr_hyp->tophad1_v4().pt();
  int nsubjets = hadr_hyp->tophad1_jets().size();

  hist("Chi2")->Fill(hadr_hyp->discriminator("Chi2Hadronic"), weight);
  hist("M_tophad_own")->Fill(mTopHad,weight);
  hist("Pt_tophad_own")->Fill(ptTopHad, weight);
  hist("Pt_tophad_own_rebin")->Fill(ptTopHad, weight);
  hist("Pt_tophad_own_binned")->Fill(ptTopHad, weight);
  hist("Pt_tophad_own_binned2")->Fill(ptTopHad, weight);
  hist("Pt_tophad_own_test")->Fill(ptTopHad, weight);
  hist("top_recjets")->Fill(nsubjets,weight);


  double deltaRtop=10;
  if(!is_data){
    for(auto genp : *event.genparticles){
      if(fabs(genp.pdgId())==6){
	double tmp = deltaR(genp,hadr_hyp->tophad1_v4());
	if(tmp<deltaRtop){
	  deltaRtop = tmp;
	}
      }
    }
    hist("deltaR_top_gentop")->Fill(deltaRtop, weight);
  }

  
  double wchi2 = hadr_hyp->wchi2();
  double tchi2 = hadr_hyp->tchi2();
  hist("WChi2")->Fill(wchi2,weight);
  hist("tChi2")->Fill(tchi2,weight);


  
  if(event.taus->size()>0){
    hist("M_tophad_tau")->Fill((hadr_hyp->tophad1_v4()+(*event.taus)[0].v4()).M(), weight);
  }
  
  for(const auto & jet : *event.jets){
    hist("deltaPhi_topjet_jet")->Fill(fabs(hadr_hyp->tophad1_v4().phi()-jet.phi()));
  }
  
  for(const auto & jet : *event.jets){
    if(  fabs(hadr_hyp->tophad1_v4().phi()-jet.phi()) >2.4  ){
      for(int i=0; i<nsubjets; i++){
	auto subjet = hadr_hyp->tophad1_jets().at(i);
	if(deltaR(jet,subjet)<0.4){ continue;}
	hist("M_tophad_jet")->Fill((hadr_hyp->tophad1_v4() + jet.v4()).M(), weight);
      }
    }
  }
  

  if(!is_data){
    int Ntau = event.taus->size();
    int Nreal = 0;
    int Nfake = 0;
    double realtaupt1 = -1;
    double real_dummy = -1;
    double faketaupt1 = -1;
    double fake_dummy = -1;
    for(unsigned int i =0; i<event.taus->size(); ++i){
      Tau tau = event.taus->at(i);
      for(auto genp : *event.genparticles){
	double dR = deltaR(tau,genp);
	if(dR<0.4 && abs(genp.pdgId())==15){
	  Nreal += 1;
	  real_dummy = tau.pt();
	  if(real_dummy > realtaupt1) realtaupt1 = real_dummy;
	}
	else if(dR>=0.4 && abs(genp.pdgId())==15){
	  fake_dummy = tau.pt();
	  if(fake_dummy > faketaupt1) faketaupt1 = fake_dummy;
	}
      }
    }
    Nfake = Ntau - Nreal;
  
    if(Nreal==0){ hist("tau_categories")->Fill(0.,weight);}
    else if(Nfake==0){ hist("tau_categories")->Fill(1.,weight);}
    else if(Nreal>0 && Nfake>0){ hist("tau_categories")->Fill(2.,weight);}
    else{hist("tau_categories")->Fill(3.,weight);}

    hist("N_faketaus")->Fill(Nfake, weight);
    hist("N_realtaus")->Fill(Nreal, weight);

    if(Nfake>0 && Nreal>0){
      if(faketaupt1>realtaupt1) hist("faketaupt_realtaupt")->Fill(0., weight);
      else{ hist("faketaupt_realtaupt")->Fill(1., weight);}
    }

  }

  int N_OSpairs = 0;
  int N_SSpairs = 0;
  for(unsigned int i =0; i<event.taus->size(); ++i){
    Tau taui = event.taus->at(i);
    for(unsigned int j =0; j<event.taus->size(); ++j){
      Tau tauj = event.taus->at(j);
      if(j!=i && j>i){
	if(taui.charge()==tauj.charge()) N_SSpairs++;
	if(taui.charge()!=tauj.charge()) N_OSpairs++;
      }
    }
  }
  hist("N_OStau")->Fill(N_OSpairs, weight);
  hist("N_SStau")->Fill(N_SSpairs, weight);

  if(event.taus->size() > 1){
    const auto & tau1 = (*event.taus)[0];
    const auto & tau2 = (*event.taus)[1];
    hist("M_tautau")->Fill( (tau1.v4() + tau2.v4()).M(), weight);
    if(tau1.charge()==tau2.charge()){
      hist("tausign")->Fill(0., weight);
    }
    else{ hist("tausign")->Fill(1., weight);}
  }

  /*
  if(!is_data){
    const TTbarFullhadRecoHypothesis* hadr_hyp_match = get_best_hypothesis( hadr_hyps, "CorrectMatch" );
    double wmass_matched;
    if(hadr_hyp_match->discriminator("CorrectMatch") < 20){
      if((hadr_hyp_match->WJet1_v4()+hadr_hyp_match->WJet2_v4()).isTimelike())  wmass_matched = (hadr_hyp_match->WJet1_v4()+hadr_hyp_match->WJet2_v4()).M();
      else{ wmass_matched = -sqrt(-(hadr_hyp_match->WJet1_v4()+hadr_hyp_match->WJet2_v4()).mass2());}
      hist("WMass_matched")->Fill(wmass_matched,weight);
    }
  }
  */

  //double mWHad = hadr_hyp->whad1_v4().M();
  double mWHad;
  //if(hadr_hyps.at(0).whad1_v4().isTimelike())      mWHad = hadr_hyps.at(0).whad1_v4().M();
  //else     mWHad  = -sqrt(-hadr_hyps.at(0).whad1_v4().mass2());
  //double mWHad = hadr_hyps.at(0).whad1_v4().M();
  if(hadr_hyp->whad1_v4().isTimelike())      mWHad = hadr_hyp->whad1_v4().M();
  else     mWHad  = -sqrt(-hadr_hyp->whad1_v4().mass2());
  hist("M_whad_own")->Fill(mWHad,weight);

  
  double maxDR=0.;
  for(int i=0; i<nsubjets; i++){
    for(int j=0; j<nsubjets; j++){
      if(i!=j){
        double dr_dummy = deltaR(hadr_hyp->tophad1_jets().at(i),hadr_hyp->tophad1_jets().at(j));
        if(dr_dummy>maxDR){
          maxDR = dr_dummy;
        }
      }
    }
  }
  double maxDRfit = (-4*pow(10,-9))*pow(ptTopHad,3) + (1.16*pow(10,-5))*pow(ptTopHad,2) + (-0.01)*ptTopHad+3.66;

  hist("maxDR_minus_maxDRfit")->Fill(maxDR-maxDRfit,weight);

  hist("deltaR_recjets")->Fill(maxDR,weight);
  ((TH2F*)hist("Pt_tophad_vs_largestDR"))->Fill(ptTopHad, maxDR, weight);

  
  if(!is_data){
    const auto & ttbargen = event.get(h_ttbargen);
    auto dec = ttbargen.DecayChannel();
    if(dec==TTbarGen::e_muhad || dec==TTbarGen::e_ehad || dec==TTbarGen::e_tauhad){
     if(ttbargen.IsTopHadronicDecay()){
       double pt_top_gen = (ttbargen.bTop().v4() + ttbargen.Wdecay1().v4() + ttbargen.Wdecay2().v4()).pt();
       hist("pt_top_resolution")->Fill((ptTopHad-pt_top_gen) / pt_top_gen,weight);
       //hist("Pt_top_gen")->Fill(pt_top_gen, weight);
     }
     if(ttbargen.IsAntiTopHadronicDecay()){
       double pt_top_gen = (ttbargen.bAntitop().v4() + ttbargen.WMinusdecay1().v4() + ttbargen.WMinusdecay2().v4()).pt();
       hist("pt_top_resolution")->Fill((ptTopHad-pt_top_gen) / pt_top_gen,weight);
       //hist("Pt_top_gen")->Fill(pt_top_gen, weight);
     }
    }
    else if(dec==TTbarGen::e_had){
      double pt_top_gen1 = 0;
      double pt_top_gen2 = 0;
      if(ttbargen.IsTopHadronicDecay()){
	pt_top_gen1 = (ttbargen.bAntitop().v4() + ttbargen.WMinusdecay1().v4() + ttbargen.WMinusdecay2().v4()).pt();
      }
      if(ttbargen.IsAntiTopHadronicDecay()){
	pt_top_gen2 = (ttbargen.bAntitop().v4() + ttbargen.WMinusdecay1().v4() + ttbargen.WMinusdecay2().v4()).pt();
      }
      if(abs(ptTopHad-pt_top_gen1)<abs(ptTopHad-pt_top_gen2)){
	hist("pt_top_resolution")->Fill((ptTopHad-pt_top_gen1) / pt_top_gen1,weight);
	//hist("Pt_top_gen")->Fill(pt_top_gen1, weight);
      }
      else{
	hist("pt_top_resolution")->Fill((ptTopHad-pt_top_gen2) / pt_top_gen2,weight);
	//hist("Pt_top_gen")->Fill(pt_top_gen2, weight);
      }
    }
  }
  

  if(!is_data){
    const auto & ttbargen = event.get(h_ttbargen);
    auto dec = ttbargen.DecayChannel();
    double maxDR_gen=0.;
    double pt_top =0.;
    if(dec==TTbarGen::e_muhad || dec==TTbarGen::e_ehad || dec==TTbarGen::e_tauhad){
      if(ttbargen.IsTopHadronicDecay()){
	pt_top = (ttbargen.bTop().v4() + ttbargen.Wdecay1().v4() + ttbargen.Wdecay2().v4()).pt();
	double dr_1 = deltaR(ttbargen.bTop().v4(),ttbargen.Wdecay1().v4());
	double dr_2 = deltaR(ttbargen.bTop().v4(),ttbargen.Wdecay2().v4());
	double dr_3 = deltaR(ttbargen.Wdecay1().v4(),ttbargen.Wdecay2().v4());
	if(dr_1 > dr_2 && dr_1>dr_3){
	  maxDR_gen = dr_1;
	}
	if(dr_2 > dr_1 && dr_2>dr_3){
	  maxDR_gen = dr_2;
	}
	if(dr_3 > dr_1 && dr_3>dr_2){
	  maxDR_gen = dr_3;
	}
      }
      if(ttbargen.IsAntiTopHadronicDecay()){
	pt_top = (ttbargen.bAntitop().v4() + ttbargen.WMinusdecay1().v4() + ttbargen.WMinusdecay2().v4()).pt();
	double dr_1 = deltaR(ttbargen.bAntitop().v4(),ttbargen.WMinusdecay1().v4());
	double dr_2 = deltaR(ttbargen.bAntitop().v4(),ttbargen.WMinusdecay2().v4());
	double dr_3 = deltaR(ttbargen.WMinusdecay1().v4(),ttbargen.WMinusdecay2().v4());
	if(dr_1 > dr_2 && dr_1>dr_3){
	  maxDR_gen = dr_1;
	}
	if(dr_2 > dr_1 && dr_2>dr_3){
	  maxDR_gen = dr_2;
	}
	if(dr_3 > dr_1 && dr_3>dr_2){
	  maxDR_gen = dr_3;
	}
      }
      hist("deltaR_genjets")->Fill(maxDR_gen,weight);
      ((TH2F*)hist("Pt_tophad_vs_largestDR_gen"))->Fill(pt_top, maxDR_gen, weight);
    }
    else if(dec==TTbarGen::e_had){
      double pt_top1 = 0;
      double pt_top2 = 0;
      double maxDR_gen1 = 0;
      double maxDR_gen2 = 0;
      if(ttbargen.IsTopHadronicDecay()){
	pt_top1 = (ttbargen.bTop().v4() + ttbargen.Wdecay1().v4() + ttbargen.Wdecay2().v4()).pt();
	double dr_1 = deltaR(ttbargen.bTop().v4(),ttbargen.Wdecay1().v4());
	double dr_2 = deltaR(ttbargen.bTop().v4(),ttbargen.Wdecay2().v4());
	double dr_3 = deltaR(ttbargen.Wdecay1().v4(),ttbargen.Wdecay2().v4());
	if(dr_1 > dr_2 && dr_1>dr_3){
	  maxDR_gen1 = dr_1;
	}
	if(dr_2 > dr_1 && dr_2>dr_3){
	  maxDR_gen1 = dr_2;
	}
	if(dr_3 > dr_1 && dr_3>dr_2){
	  maxDR_gen1 = dr_3;
	}
      }
      if(ttbargen.IsAntiTopHadronicDecay()){
	pt_top2 = (ttbargen.bAntitop().v4() + ttbargen.WMinusdecay1().v4() + ttbargen.WMinusdecay2().v4()).pt();
	double dr_1 = deltaR(ttbargen.bAntitop().v4(),ttbargen.WMinusdecay1().v4());
	double dr_2 = deltaR(ttbargen.bAntitop().v4(),ttbargen.WMinusdecay2().v4());
	double dr_3 = deltaR(ttbargen.WMinusdecay1().v4(),ttbargen.WMinusdecay2().v4());
	if(dr_1 > dr_2 && dr_1>dr_3){
	  maxDR_gen2 = dr_1;
	}
	if(dr_2 > dr_1 && dr_2>dr_3){
	  maxDR_gen2 = dr_2;
	}
	if(dr_3 > dr_1 && dr_3>dr_2){
	  maxDR_gen2 = dr_3;
	}
      }
      if(maxDR_gen1>maxDR_gen2){
	maxDR_gen = maxDR_gen1;
	pt_top = pt_top1;
      }
      else{ 
	maxDR_gen = maxDR_gen2;
	pt_top = pt_top2;
      }
      hist("deltaR_genjets")->Fill(maxDR_gen,weight);
      ((TH2F*)hist("Pt_tophad_vs_largestDR_gen"))->Fill(pt_top, maxDR_gen, weight);
    }
  }

  /*
  std::vector<ReconstructionHypothesis> hyps = event.get(h_ttbar_hyps); 
  const ReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );

  hist("pt_thad")->Fill( hyp->tophad_v4().Pt(),weight );
  hist("pt_tlep")->Fill( hyp->toplep_v4().Pt(),weight );

  hist("pt_tcom")->Fill( hyp->tophad_v4().Pt(),weight );
  hist("pt_tcom")->Fill( hyp->toplep_v4().Pt(),weight );
  */

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


  int n_isolated=0;
  int n_nonisolated=0;
  for(const auto & tau : *event.taus){
    hist("isolation")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
    hist("isolation_sideband")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
    hist("isolation_fullrange")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
    hist("pt_tau_binned")->Fill(tau.pt(),weight);
    hist("MVA_isolation")->Fill(tau.byIsolationMVArun2v1DBnewDMwLTraw(),weight);
    if(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()<1.5){
      n_isolated++;
      hist("isolated_taus")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
    }
    if(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()>=1.5){
      n_nonisolated++;
      hist("nonisolated_taus")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
    }
  }
  hist("N_isolated")->Fill(n_isolated, weight);
  hist("N_nonisolated")->Fill(n_nonisolated, weight);


  
  const auto muons = event.muons;
  if(muons->size()>1){
    for(unsigned int i=0; i<muons->size(); ++i) 
      {
	Muon muon1 = muons->at(i);
	TLorentzVector Mu1;
	Mu1.SetPtEtaPhiE(muon1.pt() ,muon1.eta() ,muon1.phi() ,muon1.energy() );
	for(unsigned int j=0; j<muons->size(); ++j) 
	  {
	    Muon muon2 = muons->at(j);
	    TLorentzVector Mu2;
	    Mu2.SetPtEtaPhiE(muon2.pt() ,muon2.eta() ,muon2.phi() ,muon2.energy() );
	    TLorentzVector Vec =  Mu1+Mu2;
	    double InvMass = Vec.M(); 
	    hist("M_mumu")->Fill(InvMass, weight);	
	  }
      }
  }

  const auto taus = event.taus;
  if(taus->size()>0){
    hist("pt_tau1")->Fill((*event.taus)[0].pt(), weight);
  }




  if(muons->size()>0){
    if(taus->size()>0){
        const auto & mu1 = (*event.muons)[0];
        const auto & ta1 = (*event.taus)[0];
        hist("DeltaR_mutau")->Fill(deltaR(mu1,ta1) ,weight);
	for(unsigned int i=0; i<muons->size(); ++i) 
	  {
	    Muon muon = muons->at(i);
	    TLorentzVector Mu;
	    Mu.SetPtEtaPhiE(muon.pt() ,muon.eta() ,muon.phi() ,muon.energy() );
	    for(unsigned int j=0; j<taus->size(); ++j) 
	      {
		Tau tau = taus->at(j);
		TLorentzVector Tau;
		Tau.SetPtEtaPhiE(tau.pt() ,tau.eta() ,tau.phi() ,tau.energy() );
		TLorentzVector Vec =  Mu+Tau;
		double InvMass = Vec.M(); 
		hist("M_mutau")->Fill(InvMass, weight);	
	      }
	  }
      }
  }
  


  for(const auto & muon : *event.muons){
    hist("MT")->Fill(sqrt(2*muon.pt()*event.met->pt()* (1-cos(deltaPhi(*event.met,muon)) )), weight);
    hist("deltaPhi_mu_met")->Fill(deltaPhi(*event.met,muon) ,weight);
    break;
  }

  const auto jets = event.jets;

  if(jets->size()>0){
    hist("pt_firstjet")->Fill((*event.jets)[0].pt(), weight);
  }

  vector<Jet> bjetsL, bjetsM, bjetsT;
  for (unsigned int i =0; i<jets->size(); ++i) {
    if(jets->at(i).btag_combinedSecondaryVertex()>0.460) {
      bjetsL.push_back(jets->at(i));
    }
  }
  for (unsigned int i =0; i<jets->size(); ++i) {
    if(jets->at(i).btag_combinedSecondaryVertex()>0.800) {
      bjetsM.push_back(jets->at(i));
    }
  }
  for (unsigned int i =0; i<jets->size(); ++i) {
    if(jets->at(i).btag_combinedSecondaryVertex()>0.935) {
      bjetsT.push_back(jets->at(i));
    }
  }
  int NbJetsL = bjetsL.size();
  int NbJetsM = bjetsM.size();
  int NbJetsT = bjetsT.size();
  hist("NbJetsL")-> Fill(NbJetsL,weight);
  hist("NbJetsM")-> Fill(NbJetsM,weight);
  hist("NbJetsT")-> Fill(NbJetsT,weight);

  ////////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  double minMT = 1000000;
  for (unsigned int i =0; i<=bjetsM.size(); ++i) {
    if (bjetsM.size()> i) {
      Jet bjet = bjetsM[i];
      double tmp =  sqrt(2*bjet.pt()*event.met->pt()* (1-cos(event.met->phi()-bjet.phi())) );
      if(tmp<minMT){
	minMT = tmp;
      }
    }
  }


  double minMT_METtau = 1000000;
  for(const auto & tau : *event.taus){
    double tmp =  sqrt(2*tau.pt()*event.met->pt()* (1-cos(event.met->phi()-tau.phi())) );
    if(tmp<minMT_METtau){
      minMT_METtau = tmp;
    }
  }
  hist("h_MT_METtau")->Fill(minMT_METtau,weight);

  double minMT_btau = 10000000;
  for (unsigned int i=0; i<=(*taus).size(); ++i){
    if((*taus).size()>i){
      Tau tau = (*taus)[i];
      TLorentzVector Tau;
      Tau.SetPtEtaPhiE(tau.pt() ,tau.eta() ,tau.phi() ,tau.energy() );
      for (unsigned int i =0; i<=bjetsM.size(); ++i) {
	if (bjetsM.size()> i) {
	  Jet bjet = bjetsM[i];
	  TLorentzVector BJet;
	  BJet.SetPtEtaPhiE(bjet.pt() ,bjet.eta() ,bjet.phi() ,bjet.energy() );
	  hist("M_btau")->Fill((Tau+BJet).M(),weight);
	  double tmp =  sqrt(2*bjet.pt()*tau.pt()* (1-cos(tau.phi()-bjet.phi())) );
	  if(tmp<minMT_btau){
	    minMT_btau = tmp; 
	  }
	}
      }
    }
  }

  if(NbJetsM>0){
    hist("h_minMT")->Fill(minMT,weight);
    hist("h_MT_btau")->Fill(minMT_btau,weight);
  }

  if(event.taus->size() > 0){
    const auto & tau = (*event.taus)[0];
    for (unsigned int i =0; i<=bjetsM.size(); ++i) {
      if (bjetsM.size()> i) {
	Jet bjet = bjetsM[i];
	double deltaR_btau = deltaR(bjet,tau);
	hist("DeltaR_btau")->Fill(deltaR_btau,weight);
      }
    }
  }

  for (unsigned int i =0; i<=bjetsM.size(); ++i) {
    if (bjetsM.size()> i) {
      Jet bjet = bjetsM[i];
      hist("pt_b")->Fill(bjet.pt(),weight);
      hist("pt_b_binned")->Fill(bjet.pt(),weight);
    }
  }
  

  
  double sum_ele=0;
  double sum_mu=0;
  double sum_tau=0;
  for(const auto & ele : *event.electrons){
    sum_ele+=TMath::ATan(exp(-fabs(ele.eta())));
  }
  for(const auto & muon : *event.muons){
    sum_mu+=TMath::ATan(exp(-fabs(muon.eta())));
  }
  for(const auto & tau : *event.taus){
    sum_tau+=TMath::ATan(exp(-fabs(tau.eta())));
  }
  

  double sum_leptons=(event.electrons->size()+event.muons->size()+event.taus->size());

  //cout << 1/sum_leptons << endl;
  double eta_halil=-TMath::Log(TMath::Tan((1/sum_leptons)*(sum_ele+sum_mu+sum_tau)));
  //cout << eta_halil << endl;

  hist("eta_tilde")->Fill( eta_halil ,weight);
  hist("eta_tilde2")->Fill( eta_halil ,weight);
  if(eta_halil<0.9)  hist("eta_tilde_bin1")->Fill( eta_halil ,weight);
  if(eta_halil>=0.9)  hist("eta_tilde_bin2")->Fill( eta_halil ,weight);


  if(!is_data){
    double dR = 1000;
    for(const auto & tau : *event.taus){
      for(auto genp : *event.genparticles){
	//if(abs(genp.pdgId())!=15) continue;
	if(abs(genp.pdgId())==15){
	  double tmp = deltaR(tau,genp);
	  if(tmp<dR){
	    dR = tmp;
	  }
	}
      }
    }
    if(dR<0.4){
      hist("tau_type")->Fill(0.,weight);
      if(event.taus->size() > 0){
	const auto & tau = (*event.taus)[0];
	hist("pt_real_tau1_binned")->Fill(tau.pt(), weight);
	//hist("pt_tau1_binned")->Fill(tau.pt(), weight);
      }
    }
    else{
      hist("tau_type")->Fill(1.,weight);
      if(event.taus->size() > 0){
	const auto & tau = (*event.taus)[0];
	hist("pt_fake_tau1_binned")->Fill(tau.pt(), weight);
	//hist("pt_tau1_binned")->Fill(tau.pt(), weight);
      }
    }
  }

  
  if(event.taus->size() > 0){
    const auto & tau = (*event.taus)[0];
    hist("pt_tau1_binned")->Fill(tau.pt(), weight);
    ((TH2F*)hist("pt_tau1_vs_ST"))->Fill(tau.pt(),ht+ht_lep+met, weight);
    ((TH2F*)hist("pt_top_vs_ST"))->Fill(ptTopHad, ht+ht_lep+met,weight);
    ((TH2F*)hist("pt_top_vs_pt_tau1"))->Fill(ptTopHad, tau.pt(),weight);
  }    
  
  if(!is_data){
    for(const auto & muon : *event.muons){
      double dR = 1000;
      for(auto genp : *event.genparticles){
	if(abs(genp.pdgId())==13){
	  double tmp = deltaR(muon,genp);
	  if(tmp<dR){
	    dR = tmp;
	  }
	}
      }
      if(dR<0.1){
	hist("muon_type")->Fill(0.,weight);
      }
      else{
	hist("muon_type")->Fill(1.,weight);
      }
    }
  }

  if(!is_data){
    for(const auto & electron : *event.electrons){
      double dR = 1000;
      for(auto genp : *event.genparticles){
	if(abs(genp.pdgId())==11){
	  double tmp = deltaR(electron,genp);
	  if(tmp<dR){
	    dR = tmp;
	  }
	}
      }
      if(dR<0.1){
	hist("electron_type")->Fill(0.,weight);
      }
      else{
	hist("electron_type")->Fill(1.,weight);
      }
    }
  }


  //CMSTopTags
  
  double mDiminLower = 50., mjetLower = 140., mjetUpper = 250.;
  //std::vector<TopJet>* topjets = event.topjets;
  //std::vector<TopJet> taggedtopjets;
  int N_toptaggedjets = 0;
  bool CMSTopTag = true;

  double m_disubjet_min = 0.;
  
  for(const auto & topjet : *event.topjets){

    std::vector<Jet> subjets = topjet.subjets();
    
    if(subjets.size() < 2) m_disubjet_min = 0.0;
    
    // only need to sort if subjets there are more than 3 subjets, as
    // otherwise, we use all 3 anyway.
    if(subjets.size() > 3) sort_by_pt(subjets);
    
    double m01 = 0;
    LorentzVector sum01test = subjets[0].v4()+subjets[1].v4();
    LorentzVector sum02test = subjets[0].v4()+subjets[2].v4();
    LorentzVector sum12test = subjets[1].v4()+subjets[2].v4();
    double m01pt = 0;
    double m01eta = 0;
    double m01phi = 0;
    double m01energy = 0;
    m01pt = sum01test.pt();
    m01eta = sum01test.eta();
    m01phi = sum01test.phi();
    m01energy = sum01test.energy();
    double m02pt = 0;
    double m02eta = 0;
    double m02phi = 0;
    double m02energy = 0;
    m02pt = sum02test.pt();
    m02eta = sum02test.eta();
    m02phi = sum02test.phi();
    m02energy = sum02test.energy();
    double m12pt = 0;
    double m12eta = 0;
    double m12phi = 0;
    double m12energy = 0;
    m12pt = sum12test.pt();
    m12eta = sum12test.eta();
    m12phi = sum12test.phi();
    m12energy = sum12test.energy();

    TLorentzVector sum01;
    sum01.SetPtEtaPhiE(m01pt ,m01eta ,m01phi ,m01energy );
    
    /*if(sum01.isTimelike())  */m01 = sum01.M();

    
    if(subjets.size() < 3) m_disubjet_min = m01;
    
    double m02 = 0;
    //auto sum02 = subjets[0].v4()+subjets[2].v4();
    TLorentzVector sum02;
    sum02.SetPtEtaPhiE(m02pt ,m02eta ,m02phi ,m02energy );
    /*if( sum02.isTimelike() )*/ m02 = sum02.M();
    
    double m12 = 0;
    //auto sum12 = subjets[1].v4()+subjets[2].v4();
    TLorentzVector sum12;
    sum12.SetPtEtaPhiE(m12pt ,m12eta ,m12phi ,m12energy );
    /*if( sum12.isTimelike() )*/  m12 = sum12.M();
    
    
    m_disubjet_min = std::min(m01,std::min(m02, m12));
    hist("min_mDisubjet")->Fill(m_disubjet_min, weight);
    if(m_disubjet_min < mDiminLower) CMSTopTag = false;
    
    //auto mjet = topjet.v4().M();
    TLorentzVector mjetv4;
    mjetv4.SetPtEtaPhiE(topjet.pt() ,topjet.eta() ,topjet.phi() ,topjet.energy() );

    double mjet = mjetv4.M();

    hist("M_jet")->Fill(mjet, weight);
    if(mjet < mjetLower) CMSTopTag = false;
    if(mjet > mjetUpper) CMSTopTag = false;
    
    hist("N_subjets")->Fill(subjets.size(), weight);
    if(subjets.size() < 3) CMSTopTag = false;
    
    //if (CMSTopTag) taggedtopjets.push_back(topjet); 
    if (CMSTopTag) N_toptaggedjets++; 
    
  }

  
 //int N_toptaggedjets = taggedtopjets.size();
  hist("N_TopTags")->Fill(N_toptaggedjets, weight);
  

  if(event.muons->size() > 0 && event.taus->size() > 0){
    const auto & muon = (*event.muons)[0];
    const auto & tau = (*event.taus)[0];
    if (muon.charge() == tau.charge()){
      hist("sign")->Fill(-0.5, weight);
    }
    else{
      hist("sign")->Fill(0.5, weight);
    }
  }


  hist("N")->Fill(1,weight);

  /*
  if(!is_data){
    for(const auto & gp : *event.genparticles){
      if(abs(gp.pdgId()) == 6){ 
	// now get W daughters:
	auto topd1 = gp.daughter(event.genparticles, 1);
	auto topd2 = gp.daughter(event.genparticles, 2);
	auto wd1=gp.daughter(event.genparticles, 1);
	auto wd2=gp.daughter(event.genparticles, 1);
	if(abs(topd1->pdgId())==24){
	  wd1 = topd1->daughter(event.genparticles, 1);
	  wd2 = topd1->daughter(event.genparticles, 2);
	}
	if(abs(topd2->pdgId())==24){
	  wd1 = topd2->daughter(event.genparticles, 1);
	  wd2 = topd2->daughter(event.genparticles, 2);
	}
	if(abs(wd1->pdgId())==11 || abs(wd2->pdgId())==11) n_ele++;
	if(abs(wd1->pdgId())==13 || abs(wd2->pdgId())==13 ) n_mu++;
	if(abs(wd1->pdgId())==15 || abs(wd2->pdgId())==15) n_tau++;
	for(unsigned int i =0; i<event.jets->size(); ++i){
	  Jet jet = event.jets->at(i);
	  double dR1 = deltaR(jet,*topd2);
	  double dR2 = deltaR(jet,*wd1);
	  double dR3 = deltaR(jet,*wd2);
	  if(dR1<0.4 && abs(topd2->pdgId())<7) n_j1++;
	  if(dR2<0.4 && abs(wd1->pdgId())<7) n_j2++;
	  if(dR3<0.4 && abs(wd2->pdgId())<7) n_j3++;
	}
      }
    }
  }
  */


  
  //  Jet bjet1 = bjets[0];
  /*TLorentzVector BJet1;
  BJet1.SetPtEtaPhiE(bjet1.pt() ,bjet1.eta() ,bjet1.phi() ,bjet1.energy() );
  hist("M_b1tau")->Fill((Tau1+BJet1).M(),weight);
  */

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
