#include "UHH2/LQAnalysis/include/LQFakeTauHists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"


#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
//#include "TMath.h"

using namespace std;
using namespace uhh2;

LQFakeTauHists::LQFakeTauHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  Double_t taubins[5] = {20, 60, 120, 200, 800};
  Double_t moretaubins[8] = {20, 40, 60, 90, 120, 160, 200, 800};
  book<TH1F>("tau_type", "0 real tau, 1 fake tau", 2,-0.5,1.5);
  book<TH1F>("pt_real_tau1_binned","p_{T} real tau 1",4,taubins);
  book<TH1F>("pt_fake_tau","p_{T} fake taus",50,0,800);
  book<TH1F>("pt_fake_tau1","p_{T} fake tau 1",50,0,800);
  book<TH1F>("pt_fake_tau1_binned","p_{T} fake tau 1",4,taubins);
  book<TH1F>("pt_fake_tau1_doublebinned","p_{T} fake tau 1",7,moretaubins);
  book<TH1F>("pt_tau1_binned","p_{T} tau 1",4,taubins);

  book<TH2F>("pt_tau1_vs_iso","pt tau1 vs tau iso", 50, 0, 800 ,50, 0, 500);
  book<TH2F>("eta_tau1_vs_iso","eta tau1 vs tau iso", 50, -3, 3 ,50, 0, 500);

  book<TH1F>("faketau1_iso","tau_iso",50,0,500);
  book<TH1F>("faketau1_pt20to40_iso","tau_iso 20to40",50,0,500);
  book<TH1F>("faketau1_pt40to60_iso","tau_iso 40to60",50,0,500);
  book<TH1F>("faketau1_pt60to90_iso","tau_iso 60to90",50,0,500);
  book<TH1F>("faketau1_pt90to120_iso","tau_iso 90to120",50,0,500);
  book<TH1F>("faketau1_pt120to160_iso","tau_iso 120to160",50,0,500);
  book<TH1F>("faketau1_pt160to200_iso","tau_iso 160to200",50,0,500);
  book<TH1F>("faketau1_pt200to800_iso","tau_iso 200to800",50,0,500);

  book<TH1F>("pt_fake_gluon_tau1_binned","p_{T} fake tau 1 (gluon)",4,taubins);
  book<TH1F>("pt_fake_gluon_tau1_doublebinned","p_{T} fake tau 1 (gluon)",7,moretaubins);
  book<TH1F>("pt_fake_Up_SS_tau1_binned","p_{T} fake tau 1 (up type quark SS)",4,taubins);
  book<TH1F>("pt_fake_Up_SS_tau1_doublebinned","p_{T} fake tau 1 (up type quark SS)",7,moretaubins);
  book<TH1F>("pt_fake_Down_SS_tau1_binned","p_{T} fake tau 1 (fake down type quark SS)",4,taubins);
  book<TH1F>("pt_fake_Down_SS_tau1_doublebinned","p_{T} fake tau 1 (fake down type quark SS)",7,moretaubins);
  book<TH1F>("pt_fake_ChargeFlip_tau1_binned","p_{T} fake tau 1 (chargeflip)",4,taubins);
  book<TH1F>("pt_fake_ChargeFlip_tau1_doublebinned","p_{T} fake tau 1 (chargeflip)",7,moretaubins);
  book<TH1F>("pt_fake_others_tau1_binned","p_{T} fake tau 1 (others)",4,taubins);
  book<TH1F>("pt_fake_others_tau1_doublebinned","p_{T} fake tau 1 (others)",7,moretaubins);

  book<TH1F>("pt_fake_uquark_tau1_doublebinned","p_{T} fake tau 1 (u quarks)",7,moretaubins);
  book<TH1F>("pt_fake_cquark_tau1_doublebinned","p_{T} fake tau 1 (c quarks)",7,moretaubins);
  book<TH1F>("pt_fake_dquark_tau1_doublebinned","p_{T} fake tau 1 (d and s quarks)",7,moretaubins);
  book<TH1F>("pt_fake_bquark_tau1_doublebinned","p_{T} fake tau 1 (b quarks)",7,moretaubins);

  book<TH1F>("ptgen_fake_gluon_tau1_doublebinned","p_{T,gen} fake tau 1 (gluons)",7,moretaubins);
  book<TH1F>("ptgen_fake_uquark_tau1_doublebinned","p_{T,gen} fake tau 1 (u quarks)",7,moretaubins);
  book<TH1F>("ptgen_fake_cquark_tau1_doublebinned","p_{T,gen} fake tau 1 (c quarks)",7,moretaubins);
  book<TH1F>("ptgen_fake_dquark_tau1_doublebinned","p_{T,gen} fake tau 1 (d and s quarks)",7,moretaubins);
  book<TH1F>("ptgen_fake_bquark_tau1_doublebinned","p_{T,gen} fake tau 1 (b quarks)",7,moretaubins);


  book<TH1F>("Gluon_pT_Tau_1","Gluon_pT_Tau_1",4,0.5,4.5);
  book<TH1F>("Gluon_pT_Tau_2","Gluon_pT_Tau_2",4,0.5,4.5);
  book<TH1F>("Gluon_pT_Tau_3","Gluon_pT_Tau_3",4,0.5,4.5);
  book<TH1F>("Gluon_pT_Tau_4","Gluon_pT_Tau_4",4,0.5,4.5);

  book<TH1F>("Up_SS_pT_Tau_1","Up_SS_pT_Tau_1",4,0.5,4.5);
  book<TH1F>("Up_SS_pT_Tau_2","Up_SS_pT_Tau_2",4,0.5,4.5);
  book<TH1F>("Up_SS_pT_Tau_3","Up_SS_pT_Tau_3",4,0.5,4.5);
  book<TH1F>("Up_SS_pT_Tau_4","Up_SS_pT_Tau_4",4,0.5,4.5);

  book<TH1F>("Down_SS_pT_Tau_1","Down_SS_pT_Tau_1",4,0.5,4.5);
  book<TH1F>("Down_SS_pT_Tau_2","Down_SS_pT_Tau_2",4,0.5,4.5);
  book<TH1F>("Down_SS_pT_Tau_3","Down_SS_pT_Tau_3",4,0.5,4.5);
  book<TH1F>("Down_SS_pT_Tau_4","Down_SS_pT_Tau_4",4,0.5,4.5);

  book<TH1F>("ChargeFlip_pT_Tau_1","ChargeFlip_pT_Tau_1",4,0.5,4.5);
  book<TH1F>("ChargeFlip_pT_Tau_2","ChargeFlip_pT_Tau_2",4,0.5,4.5);
  book<TH1F>("ChargeFlip_pT_Tau_3","ChargeFlip_pT_Tau_3",4,0.5,4.5);
  book<TH1F>("ChargeFlip_pT_Tau_4","ChargeFlip_pT_Tau_4",4,0.5,4.5);


  
  book<TH1F>("dquark_pT_Tau_1","dquark_pT_Tau_1",4,0.5,4.5);
  book<TH1F>("dquark_pT_Tau_2","dquark_pT_Tau_2",4,0.5,4.5);
  book<TH1F>("dquark_pT_Tau_3","dquark_pT_Tau_3",4,0.5,4.5);
  book<TH1F>("dquark_pT_Tau_4","dquark_pT_Tau_4",4,0.5,4.5);

  book<TH1F>("bquark_pT_Tau_1","bquark_pT_Tau_1",4,0.5,4.5);
  book<TH1F>("bquark_pT_Tau_2","bquark_pT_Tau_2",4,0.5,4.5);
  book<TH1F>("bquark_pT_Tau_3","bquark_pT_Tau_3",4,0.5,4.5);
  book<TH1F>("bquark_pT_Tau_4","bquark_pT_Tau_4",4,0.5,4.5);

  


  book<TH1F>("Gluon_pT_Tau","Gluon_pT_Tau",1,1,2);
  book<TH1F>("Up_SS_pT_Tau","Up_SS_pT_Tau",1,1,2);
  book<TH1F>("Down_SS_pT_Tau","Down_SS_pT_Tau",1,1,2);
  book<TH1F>("ChargeFlip_pT_Tau","ChargeFlip_pT_Tau",1,1,2);


  auto dataset_type = ctx.get("dataset_type");
  is_data = dataset_type == "DATA";

}


void LQFakeTauHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  //hist("Weights")->Fill(1,weight);
  

  


  //const auto taus = event.taus;
  const auto & muon = (*event.muons)[0];
  //const auto & tau = (*event.taus)[0];

  /*
  for(const auto & tau : *event.taus){
    bool fake = true;
    for(auto genp : *event.genparticles){
	double dR = deltaR(tau,genp);
	if(dR<0.4 && abs(genp.pdgId())==15){
	  fake = false;
	  break;
	}
    }
  */

  /*
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
  */

  //for(const auto & tau : *event.taus){
  if(!is_data){
  for(unsigned int i =0; i<event.taus->size(); ++i){
      Tau tau = event.taus->at(i);
      bool fake = true;
      for(auto genp : *event.genparticles){
	double dR = deltaR(tau,genp);
	if(dR<0.4 && abs(genp.pdgId())==15){
	  fake = false;
	  break;
	}
      }
      if(!fake){
	//if(dR<0.4){
	hist("tau_type")->Fill(0);
	if(event.taus->size() > 0){
	  if(i==0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_real_tau1_binned")->Fill(tau.pt(), weight);
	    //hist("pt_tau1_binned")->Fill(tau.pt(), weight);
	  }
	}
      }
      if(fake){
	//else{
	hist("tau_type")->Fill(1);

	hist("pt_fake_tau")->Fill(tau.pt(),weight);

	if(i!=0) continue;

	if(event.taus->size() > 0){
	  const auto & tau = (*event.taus)[0];
	  ((TH2F*)hist("pt_tau1_vs_iso"))->Fill(tau.pt(),tau.byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
	  ((TH2F*)hist("eta_tau1_vs_iso"))->Fill(tau.eta(),tau.byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
	  hist("pt_fake_tau1")->Fill(tau.pt(), weight);
	  hist("pt_fake_tau1_binned")->Fill(tau.pt(), weight);
	  hist("pt_fake_tau1_doublebinned")->Fill(tau.pt(), weight);
	  //hist("pt_tau1_binned")->Fill(tau.pt(), weight);
	  hist("faketau1_iso")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/*/tau.pt()*/, weight);
	  if(tau.pt()<40)    hist("faketau1_pt20to40_iso")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/*/tau.pt()*/, weight);
	  if(tau.pt()>=40 && tau.pt()<60)    hist("faketau1_pt40to60_iso")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/*/tau.pt()*/, weight);
	  if(tau.pt()>=60 && tau.pt()<90)    hist("faketau1_pt60to90_iso")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/*/tau.pt()*/, weight);
	  if(tau.pt()>=90 && tau.pt()<120)    hist("faketau1_pt90to120_iso")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/*/tau.pt()*/, weight);
	  if(tau.pt()>=120 && tau.pt()<160)    hist("faketau1_pt120to160_iso")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/*/tau.pt()*/, weight);
	  if(tau.pt()>=160 && tau.pt()<200)    hist("faketau1_pt160to200_iso")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/*/tau.pt()*/, weight);
	  if(tau.pt()>=200 && tau.pt()<800)    hist("faketau1_pt200to800_iso")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/*/tau.pt()*/, weight);

	}


	double delR = 1000;
	int genpPdg = 0;
	for(auto genp_new : *event.genparticles){
	  double tmp = deltaR(tau,genp_new);
	  if(tmp<delR){
	    delR = tmp;
	    if(delR<0.4){
	      genpPdg = genp_new.pdgId();
	      if(fabs(genp_new.pdgId())==2){
		hist("ptgen_fake_uquark_tau1_doublebinned")->Fill(genp_new.pt(), weight);
	      }
	      if(fabs(genp_new.pdgId())==4){
		hist("ptgen_fake_cquark_tau1_doublebinned")->Fill(genp_new.pt(), weight);
	      }
	      if(fabs(genp_new.pdgId())==1 || fabs(genp_new.pdgId())==3){
		hist("ptgen_fake_dquark_tau1_doublebinned")->Fill(genp_new.pt(), weight);
	      }
	      if(fabs(genp_new.pdgId())==5){
		hist("ptgen_fake_bquark_tau1_doublebinned")->Fill(genp_new.pt(), weight);
	      }
	      if(genp_new.pdgId()==21){
		hist("ptgen_fake_gluon_tau1_doublebinned")->Fill(genp_new.pt(), weight);
	      }
	    }
	  }
	}
      


	//double delR = deltaR(tau,genp_new);


	//if((genp_new.pdgId()!=21) && fabs(genp_new.pdgId())!=1 && fabs(genp_new.pdgId()) !=2 && fabs(genp_new.pdgId()) !=3 && fabs(genp_new.pdgId())!=4 && fabs(genp_new.pdgId())!=5) cout << genp_new.pdgId() << endl;
	/*if(event.taus->size() > 0){
	  const auto & tau = (*event.taus)[0];
	  if((genp_new.pdgId()!=21) && fabs(genp_new.pdgId())!=1 && fabs(genp_new.pdgId()) !=2 && fabs(genp_new.pdgId()) !=3 && fabs(genp_new.pdgId())!=4 && fabs(genp_new.pdgId())!=5 && fabs(genp_new.pdgId())!=24) cout << tau.pt() << endl;
	  }*/

	if(genpPdg==21){
	  
	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_gluon_tau1_binned")->Fill(tau.pt(), weight);
	    hist("pt_fake_gluon_tau1_doublebinned")->Fill(tau.pt(), weight);
	  }
	  
	  hist("Gluon_pT_Tau")->Fill(1, weight);
	  if(event.muons->size()>0){
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Gluon_pT_Tau_1")->Fill(1, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Gluon_pT_Tau_1")->Fill(2, weight);
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Gluon_pT_Tau_1")->Fill(3, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Gluon_pT_Tau_1")->Fill(4, weight);
	    
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Gluon_pT_Tau_2")->Fill(1, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Gluon_pT_Tau_2")->Fill(2, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Gluon_pT_Tau_2")->Fill(3, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Gluon_pT_Tau_2")->Fill(4, weight);   
	    
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Gluon_pT_Tau_3")->Fill(1, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Gluon_pT_Tau_3")->Fill(2, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Gluon_pT_Tau_3")->Fill(3, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Gluon_pT_Tau_3")->Fill(4, weight);  
	    
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Gluon_pT_Tau_4")->Fill(1, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Gluon_pT_Tau_4")->Fill(2, weight);
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Gluon_pT_Tau_4")->Fill(3, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Gluon_pT_Tau_4")->Fill(4, weight);   
	  }
	}
	if (((genpPdg ==2 || genpPdg ==4) && tau.charge()==1) || ((genpPdg ==-2 || genpPdg ==-4) && tau.charge()==-1)){

	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_Up_SS_tau1_binned")->Fill(tau.pt(), weight);
	    hist("pt_fake_Up_SS_tau1_doublebinned")->Fill(tau.pt(), weight);
	  }

	  hist("Up_SS_pT_Tau")->Fill(1, weight);
	  if(event.muons->size()>0){
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Up_SS_pT_Tau_1")->Fill(1, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Up_SS_pT_Tau_1")->Fill(2, weight);
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Up_SS_pT_Tau_1")->Fill(3, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Up_SS_pT_Tau_1")->Fill(4, weight);   
	    
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Up_SS_pT_Tau_2")->Fill(1, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Up_SS_pT_Tau_2")->Fill(2, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Up_SS_pT_Tau_2")->Fill(3, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Up_SS_pT_Tau_2")->Fill(4, weight);   
            
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Up_SS_pT_Tau_3")->Fill(1, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Up_SS_pT_Tau_3")->Fill(2, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Up_SS_pT_Tau_3")->Fill(3, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Up_SS_pT_Tau_3")->Fill(4, weight);   
                           
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Up_SS_pT_Tau_4")->Fill(1, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Up_SS_pT_Tau_4")->Fill(2, weight);
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Up_SS_pT_Tau_4")->Fill(3, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Up_SS_pT_Tau_4")->Fill(4, weight);   
	  }
	}
                     
	if (((genpPdg ==1 || genpPdg ==3 || genpPdg ==5) && tau.charge()==1) || ((genpPdg ==-1 || genpPdg ==-3 || genpPdg ==-5) && tau.charge()==-1)){

	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_Down_SS_tau1_binned")->Fill(tau.pt(), weight);
	    hist("pt_fake_Down_SS_tau1_doublebinned")->Fill(tau.pt(), weight);
	  }

	  hist("Down_SS_pT_Tau")->Fill(1, weight);
	  if(event.muons->size()>0){
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Down_SS_pT_Tau_1")->Fill(1, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Down_SS_pT_Tau_1")->Fill(2, weight);
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Down_SS_pT_Tau_1")->Fill(3, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Down_SS_pT_Tau_1")->Fill(4, weight);   
                           
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Down_SS_pT_Tau_2")->Fill(1, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Down_SS_pT_Tau_2")->Fill(2, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Down_SS_pT_Tau_2")->Fill(3, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Down_SS_pT_Tau_2")->Fill(4, weight);   
                           
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Down_SS_pT_Tau_3")->Fill(1, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Down_SS_pT_Tau_3")->Fill(2, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Down_SS_pT_Tau_3")->Fill(3, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Down_SS_pT_Tau_3")->Fill(4, weight);   
                           
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("Down_SS_pT_Tau_4")->Fill(1, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("Down_SS_pT_Tau_4")->Fill(2, weight);
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("Down_SS_pT_Tau_4")->Fill(3, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("Down_SS_pT_Tau_4")->Fill(4, weight);   
	  }
	}
	if (((genpPdg ==1 || genpPdg ==3 || genpPdg ==5) && tau.charge()==-1) || ((genpPdg ==-1 || genpPdg ==-3 || genpPdg ==-5) && tau.charge()==1) || ((genpPdg ==2 || genpPdg ==4) && tau.charge()==-1) || ((genpPdg ==-2 || genpPdg ==-4) && tau.charge()==1)){

	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_ChargeFlip_tau1_binned")->Fill(tau.pt(), weight);
	    hist("pt_fake_ChargeFlip_tau1_doublebinned")->Fill(tau.pt(), weight);
	  }

	  hist("ChargeFlip_pT_Tau")->Fill(1, weight);
	  if(event.muons->size()>0){
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()==tau.charge())   hist("ChargeFlip_pT_Tau_1")->Fill(1, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("ChargeFlip_pT_Tau_1")->Fill(2, weight);
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("ChargeFlip_pT_Tau_1")->Fill(3, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("ChargeFlip_pT_Tau_1")->Fill(4, weight);   
            
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()==tau.charge())   hist("ChargeFlip_pT_Tau_2")->Fill(1, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("ChargeFlip_pT_Tau_2")->Fill(2, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("ChargeFlip_pT_Tau_2")->Fill(3, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("ChargeFlip_pT_Tau_2")->Fill(4, weight);   
            
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("ChargeFlip_pT_Tau_3")->Fill(1, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("ChargeFlip_pT_Tau_3")->Fill(2, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("ChargeFlip_pT_Tau_3")->Fill(3, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("ChargeFlip_pT_Tau_3")->Fill(4, weight);   
            
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("ChargeFlip_pT_Tau_4")->Fill(1, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("ChargeFlip_pT_Tau_4")->Fill(2, weight);
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("ChargeFlip_pT_Tau_4")->Fill(3, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("ChargeFlip_pT_Tau_4")->Fill(4, weight);   
	  }
	}


	if ((genpPdg ==2 && tau.charge()==1) || (genpPdg ==-2 && tau.charge()==-1) || (genpPdg ==2 && tau.charge()==-1) || (genpPdg ==-2 && tau.charge()==1)){
	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_uquark_tau1_doublebinned")->Fill(tau.pt(), weight);
	  }
	}




	if ((genpPdg ==4 && tau.charge()==1) || (genpPdg ==-4 && tau.charge()==-1) || (genpPdg ==4 && tau.charge()==-1) || (genpPdg ==-4 && tau.charge()==1)){
	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_cquark_tau1_doublebinned")->Fill(tau.pt(), weight);
	  }

	}



	if (((genpPdg ==1 || genpPdg ==3) && tau.charge()==1) || ((genpPdg ==-1 || genpPdg ==-3) && tau.charge()==-1) || ((genpPdg ==1 || genpPdg ==3) && tau.charge()==-1) || ((genpPdg ==-1 || genpPdg ==-3) && tau.charge()==1)){

	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_dquark_tau1_doublebinned")->Fill(tau.pt(), weight);
	  }


	  if(event.muons->size()>0){
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()==tau.charge())   hist("dquark_pT_Tau_1")->Fill(1, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("dquark_pT_Tau_1")->Fill(2, weight);
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("dquark_pT_Tau_1")->Fill(3, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("dquark_pT_Tau_1")->Fill(4, weight);   
            
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()==tau.charge())   hist("dquark_pT_Tau_2")->Fill(1, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("dquark_pT_Tau_2")->Fill(2, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("dquark_pT_Tau_2")->Fill(3, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("dquark_pT_Tau_2")->Fill(4, weight);   
            
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("dquark_pT_Tau_3")->Fill(1, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("dquark_pT_Tau_3")->Fill(2, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("dquark_pT_Tau_3")->Fill(3, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("dquark_pT_Tau_3")->Fill(4, weight);   
            
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("dquark_pT_Tau_4")->Fill(1, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("dquark_pT_Tau_4")->Fill(2, weight);
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("dquark_pT_Tau_4")->Fill(3, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("dquark_pT_Tau_4")->Fill(4, weight); 
	  }
	}

	if ((genpPdg ==5 && tau.charge()==1) || (genpPdg ==-5 && tau.charge()==-1) || (genpPdg ==5 && tau.charge()==-1) || (genpPdg ==-5 && tau.charge()==1)){

	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_bquark_tau1_doublebinned")->Fill(tau.pt(), weight);
	  }

	  if(event.muons->size()>0){
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()==tau.charge())   hist("bquark_pT_Tau_1")->Fill(1, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("bquark_pT_Tau_1")->Fill(2, weight);
	    if (tau.pt() < 60 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("bquark_pT_Tau_1")->Fill(3, weight);
	    if (tau.pt() < 60 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("bquark_pT_Tau_1")->Fill(4, weight);   
            
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()==tau.charge())   hist("bquark_pT_Tau_2")->Fill(1, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("bquark_pT_Tau_2")->Fill(2, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("bquark_pT_Tau_2")->Fill(3, weight);
	    if (tau.pt() >= 60 && tau.pt()<120 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("bquark_pT_Tau_2")->Fill(4, weight);   
            
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("bquark_pT_Tau_3")->Fill(1, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("bquark_pT_Tau_3")->Fill(2, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("bquark_pT_Tau_3")->Fill(3, weight);
	    if (tau.pt() >= 120 && tau.pt()<200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("bquark_pT_Tau_3")->Fill(4, weight);   
            
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()==tau.charge())   hist("bquark_pT_Tau_4")->Fill(1, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()==tau.charge())   hist("bquark_pT_Tau_4")->Fill(2, weight);
	    if (tau.pt() >= 200 && muon.charge()==1 && muon.charge()!=tau.charge())   hist("bquark_pT_Tau_4")->Fill(3, weight);
	    if (tau.pt() >= 200 && muon.charge()==-1 && muon.charge()!=tau.charge())   hist("bquark_pT_Tau_4")->Fill(4, weight); 
	  }
	}



	if((genpPdg!=21) && fabs(genpPdg)!=1 && fabs(genpPdg) !=2 && fabs(genpPdg) !=3 && fabs(genpPdg)!=4 && fabs(genpPdg)!=5){
	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_others_tau1_binned")->Fill(tau.pt(), weight);
	    hist("pt_fake_others_tau1_doublebinned")->Fill(tau.pt(), weight);
	  }
	}


	    
	  
      }
    }
	
  }

  
  if(event.taus->size() > 0){
    const auto & tau = (*event.taus)[0];
    hist("pt_tau1_binned")->Fill(tau.pt(), weight);
  }    
  





 
  

}

LQFakeTauHists::~LQFakeTauHists(){}
