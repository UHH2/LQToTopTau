#include "UHH2/LQAnalysis/include/LQAllhadHists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"


#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
//#include "TMath.h"

using namespace std;
using namespace uhh2;

LQAllhadHists::LQAllhadHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  Double_t bins[5] = {0, 200, 350, 500, 3000};
  Double_t bins2[7] = {0, 100, 200, 300, 400, 600, 1500};

  book <TH1F>("Chi2", "#chi^{2}", 100, 0, 500);

  book <TH1F>("M_t_had", "M_{t,had} total", 50, 0, 500);
  book <TH1F>("M_t_had1", "M_{t1,had}", 50, 0, 500);
  book <TH1F>("M_t_had2", "M_{t2,had}", 50, 0, 500);

  book <TH1F>("pt_t_had", "p_{T}^{t,had}", 100, 0, 1500);
  book <TH1F>("pt_t_had_rebin", "p_{T}^{t,had}", 6, bins2);
  book <TH1F>("pt_t_had1", "p_{T}^{t1,had}", 100, 0, 1500);
  book <TH1F>("pt_t_had2", "p_{T}^{t2,had}", 100, 0, 1500);
  book <TH1F>("pt_t1_pt_t2", "p_{T}^{top1} + p_{T}^{top2}", 100, 0, 3000);
  book <TH1F>("pt_t1_pt_t2_rebin", "p_{T}^{top1} + p_{T}^{top2}", 4, bins);
  book <TH1F>("M_Tau1", "M_{#tau,1}", 50, 0, 2.5);
  book <TH1F>("M_Tau2", "M_{#tau,2}", 50, 0, 2.5);
  book <TH1F>("NewM_Tau1", "M_{#tau,1} (including met)", 50, 0, 3);
  book <TH1F>("NewM_Tau2", "M_{#tau,2} (including met)", 50, 0, 2.5);
  book <TH1F>("M_LQ1", "M_{LQ,1}", 50, 0, 2000);
  book <TH1F>("M_LQ2", "M_{LQ,2}", 50, 0, 2000);
  book <TH1F>("M_LQ1_test", "M_{LQ,1}", 50, 0, 2000);
  book <TH1F>("M_LQ2_test", "M_{LQ,2}", 50, 0, 2000);
  book <TH1F>("tau1_mother", "Tau1 mother", 100, -50.5, 49.5);
  book <TH1F>("tau2_mother", "Tau2 mother", 100, -50.5, 49.5);
  book <TH1F>("deltaRtop1gentop", "DeltaR Top1 GenTop", 100, 0, 5);
  book <TH1F>("deltaRtop2gentop", "DeltaR Top2 GenTop", 100, 0, 5);
  book <TH1F>("deltaRlq1", "DeltaR LQ1", 100, 0, 5);
  book <TH1F>("deltaRlq2", "DeltaR LQ2", 100, 0, 5);
  book <TH1F>("deltaPHImettau1", "DeltaR #tau1, MET", 100, 0, 5);
  book <TH1F>("deltaPHImettau2", "DeltaR #tau2, MET", 100, 0, 5);
  book <TH1F>("deltaETAmettau1", "DeltaR #tau1, MET", 100, 0, 5);
  book <TH1F>("deltaETAmettau2", "DeltaR #tau2, MET", 100, 0, 5);

  book <TH1F>("tausign", "DiTau, 0: SS, 1: OS", 2,0,2);
  book <TH1F>("dR_tautau", "deltaR(tau1,tau2)", 100,0,5);
  book <TH1F>("m_tautau", "M(tau1,tau2)", 100, 0, 500);


  h_hadr_hyps = ctx.get_handle<std::vector<TTbarFullhadRecoHypothesis>>("HighMassHadronicTTbarFullhadReco");

  auto dataset_type = ctx.get("dataset_type");
  is_data = dataset_type == "DATA";

}


void LQAllhadHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;


  std::vector<TTbarFullhadRecoHypothesis> hadr_hyps = event.get(h_hadr_hyps); 
  const TTbarFullhadRecoHypothesis* hadr_hyp = get_best_hypothesis( hadr_hyps, "Chi2Hadronic" );

  double mTopHad1 = hadr_hyp->tophad1_v4().M();
  double mTopHad2 = hadr_hyp->tophad2_v4().M();
  double ptTopHad1 = hadr_hyp->tophad1_v4().pt();
  double ptTopHad2 = hadr_hyp->tophad2_v4().pt();

  hist("Chi2")->Fill(hadr_hyp->discriminator("Chi2Hadronic"), weight);
  hist("M_t_had1")->Fill(mTopHad1,weight);
  hist("M_t_had2")->Fill(mTopHad2,weight);
  hist("M_t_had")->Fill(mTopHad1,weight);
  hist("M_t_had")->Fill(mTopHad2,weight);
  hist("pt_t_had1")->Fill(ptTopHad1, weight);
  hist("pt_t_had2")->Fill(ptTopHad2, weight);
  hist("pt_t_had")->Fill(ptTopHad1, weight);
  hist("pt_t_had")->Fill(ptTopHad2, weight);
  hist("pt_t_had_rebin")->Fill(ptTopHad1, weight);
  hist("pt_t_had_rebin")->Fill(ptTopHad2, weight);
  hist("pt_t1_pt_t2")->Fill(ptTopHad1+ptTopHad2, weight);
  hist("pt_t1_pt_t2_rebin")->Fill(ptTopHad1+ptTopHad2, weight);

  
  double M1_LQ = hadr_hyp->LQ1_v4().M();
  double M2_LQ = hadr_hyp->LQ2_v4().M();
  hist("M_LQ1")->Fill(M1_LQ,weight);
  hist("M_LQ2")->Fill(M2_LQ,weight);
  
 
  double LQ11_hyp=0;
  double LQ12_hyp=0;
  double LQ21_hyp=0;
  double LQ22_hyp=0;
  double test_M1_LQ=0;
  double test_M2_LQ=0;
  /*
  if(event.taus->size() > 1){
    const auto & tau1 = (*event.taus)[0];
    const auto & tau2 = (*event.taus)[1];
    LQ11_hyp = (hadr_hyp->tophad1_v4() + tau1.v4()).M();
    LQ12_hyp = (hadr_hyp->tophad1_v4() + tau2.v4()).M();
    LQ21_hyp = (hadr_hyp->tophad2_v4() + tau1.v4()).M();
    LQ22_hyp = (hadr_hyp->tophad2_v4() + tau2.v4()).M();
  }


  if(fabs(LQ11_hyp-LQ22_hyp) < fabs(LQ12_hyp-LQ21_hyp)){
    M1_LQ = LQ11_hyp;
    M2_LQ = LQ22_hyp;
  }
  else{
    M1_LQ = LQ12_hyp;
    M2_LQ = LQ21_hyp;
  }

  hist("M_LQ1")->Fill(M1_LQ,weight);
  hist("M_LQ2")->Fill(M2_LQ,weight);
  */


  if(event.taus->size() > 1){
    const auto & tau1 = (*event.taus)[0];
    const auto & tau2 = (*event.taus)[1];
    hist("dR_tautau")->Fill(deltaR(tau1,tau2), weight);
    if(tau1.charge()==tau2.charge()){
      hist("tausign")->Fill(0.5, weight);
    }
    else{ hist("tausign")->Fill(1.5, weight);}
    hist("m_tautau")->Fill((tau1.v4()+tau2.v4()).M(),weight);
  }



  TVector3 v3_met1;
  TVector3 v3_tau1;
  TVector3 v3_mettau1;
  TVector3 v3_met2;
  TVector3 v3_tau2;
  TVector3 v3_mettau2;
  if(event.taus->size() > 1){
    const auto & tau1 = (*event.taus)[0];
    const auto & tau2 = (*event.taus)[1];
    v3_mettau1.SetPtEtaPhi(1,tau1.v4().eta(),tau1.v4().phi());
    v3_mettau2.SetPtEtaPhi(1,tau2.v4().eta(),tau2.v4().phi());
  }

  double deltaRtop1=10;
  double deltaRtop2=10;
  double deltaRleptoquark1=10;
  double deltaRleptoquark2=10;
  if(!is_data){
    if(event.taus->size() > 1){
      const auto & tau1 = (*event.taus)[0];
      const auto & tau2 = (*event.taus)[1];
      v3_met1.SetPtEtaPhi(event.met->pt(),0,event.met->phi());
      v3_tau1.SetPtEtaPhi(tau1.v4().pt(),tau1.v4().eta(),tau1.v4().phi());
      double scalar1 = fabs(v3_met1.Dot(v3_tau1) / v3_tau1.Mag());
      v3_mettau1.SetMag(scalar1);
      TLorentzVector mettau1_test;
      mettau1_test.SetPtEtaPhiM(v3_mettau1.Pt(),tau1.v4().eta(),tau1.v4().phi(),0);
      LorentzVector mettau1;
      mettau1.SetPt(mettau1_test.Pt());
      mettau1.SetEta(mettau1_test.Eta());
      mettau1.SetPhi(mettau1_test.Phi());
      mettau1.SetE(mettau1_test.E());
      v3_met2.SetPtEtaPhi(event.met->pt(),0,event.met->phi());
      v3_tau2.SetPtEtaPhi(tau2.v4().pt(),tau2.v4().eta(),tau2.v4().phi());
      double scalar2 = fabs(v3_met2.Dot(v3_tau2) / v3_tau2.Mag());
      v3_mettau2.SetMag(scalar2);
      TLorentzVector mettau2_test;
      mettau2_test.SetPtEtaPhiM(v3_mettau2.Pt(),tau2.v4().eta(),tau2.v4().phi(),0);
      LorentzVector mettau2;
      mettau2.SetPt(mettau2_test.Pt());
      mettau2.SetEta(mettau2_test.Eta());
      mettau2.SetPhi(mettau2_test.Phi());
      mettau2.SetE(mettau2_test.E());
      hist("NewM_Tau1")->Fill((tau1.v4()+mettau1).M(),weight);
      hist("NewM_Tau2")->Fill((tau2.v4()+mettau2).M(),weight);
      hist("M_Tau1")->Fill(tau1.v4().M(),weight);
      hist("M_Tau2")->Fill(tau2.v4().M(),weight);


      LQ11_hyp = (hadr_hyp->tophad1_v4() + tau1.v4() + mettau1).M();
      LQ12_hyp = (hadr_hyp->tophad1_v4() + tau2.v4() + mettau2).M();
      LQ21_hyp = (hadr_hyp->tophad2_v4() + tau1.v4() + mettau1).M();
      LQ22_hyp = (hadr_hyp->tophad2_v4() + tau2.v4() + mettau2).M();
       if(fabs(LQ11_hyp-LQ22_hyp) < fabs(LQ12_hyp-LQ21_hyp)){
	test_M1_LQ = LQ11_hyp;
	test_M2_LQ = LQ22_hyp;
      }
      else{
	test_M1_LQ = LQ12_hyp;
	test_M2_LQ = LQ21_hyp;
      }
      hist("M_LQ1_test")->Fill(test_M1_LQ,weight);
      hist("M_LQ2_test")->Fill(test_M2_LQ,weight);


      double deltaeta1 = event.met->v4().eta()-tau1.v4().eta();
      double deltaeta2 = event.met->v4().eta()-tau2.v4().eta();
      double deltaphi1 = event.met->v4().phi()-tau1.v4().phi();
      double deltaphi2 = event.met->v4().phi()-tau2.v4().phi();
      hist("deltaETAmettau1")->Fill(sqrt(deltaeta1*deltaeta1+deltaphi1*deltaphi1),weight);
      hist("deltaETAmettau2")->Fill(sqrt(deltaeta2*deltaeta2+deltaphi2*deltaphi2),weight);
      hist("deltaPHImettau1")->Fill(fabs(tau1.v4().phi()-event.met->phi()),weight);
      hist("deltaPHImettau2")->Fill(fabs(tau2.v4().phi()-event.met->phi()),weight);
      for(auto genp : *event.genparticles){
	double dr1 = deltaR(tau1,genp);
	double dr2 = deltaR(tau2,genp);
	if(fabs(genp.pdgId())==15 && dr1<0.4){
	  const GenParticle* taumother1 = genp.mother(event.genparticles,1);
	  hist("tau1_mother")->Fill(taumother1->pdgId(),weight);
	}
	if(fabs(genp.pdgId())==15 && dr2<0.4){
	  const GenParticle* taumother2 = genp.mother(event.genparticles,1);
	  hist("tau2_mother")->Fill(taumother2->pdgId(),weight);
	}
	if(fabs(genp.pdgId())==6){
	  double tmp1 = deltaR(genp,hadr_hyp->tophad1_v4());
	  if(tmp1<deltaRtop1){
	    deltaRtop1 = tmp1;
	  }
	  double tmp2 = deltaR(genp,hadr_hyp->tophad2_v4());
	  if(tmp2<deltaRtop2){
	    deltaRtop2 = tmp2;
	  }
	}

	if(fabs(genp.pdgId())==42){
	  double tmp1 = deltaR(genp,hadr_hyp->LQ1_v4());
	  if(tmp1<deltaRleptoquark1){
	    deltaRleptoquark1 = tmp1;
	  }
	  double tmp2 = deltaR(genp,hadr_hyp->LQ2_v4());
	  if(tmp2<deltaRleptoquark2){
	    deltaRleptoquark2 = tmp2;
	  }
	}

      }
    }
  }

  hist("deltaRtop1gentop")->Fill(deltaRtop1,weight);
  hist("deltaRtop2gentop")->Fill(deltaRtop2,weight);

  hist("deltaRlq1")->Fill(deltaRleptoquark1,weight);
  hist("deltaRlq2")->Fill(deltaRleptoquark2,weight);



}

LQAllhadHists::~LQAllhadHists(){}
