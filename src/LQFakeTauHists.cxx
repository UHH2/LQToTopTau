#include "UHH2/LQAnalysis/include/LQFakeTauHists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"


#include <math.h>
#include "TH1F.h"
#include <iostream>
//#include "TMath.h"

using namespace std;
using namespace uhh2;

LQFakeTauHists::LQFakeTauHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  Double_t taubins[5] = {20, 60, 120, 200, 800};
  book<TH1F>("tau_type", "0 real tau, 1 fake tau", 2,-0.5,1.5);
  book<TH1F>("pt_real_tau1_binned","p_{T} real tau 1",4,taubins);
  book<TH1F>("pt_fake_tau1_binned","p_{T} fake tau 1",4,taubins);
  book<TH1F>("pt_tau1_binned","p_{T} tau 1",4,taubins);

  book<TH1F>("pt_fake_gluon_tau1_binned","p_{T} fake gluon tau 1",4,taubins);
  book<TH1F>("pt_fake_Up_SS_tau1_binned","p_{T} fake up type quark SS tau 1",4,taubins);
  book<TH1F>("pt_fake_Down_SS_tau1_binned","p_{T} fake down type quark SS tau 1",4,taubins);
  book<TH1F>("pt_fake_ChargeFlip_tau1_binned","p_{T} fake chargeflip tau 1",4,taubins);
  book<TH1F>("pt_fake_WBoson_tau1_binned","p_{T} fake W-Boson tau 1",4,taubins);

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

  book<TH1F>("Gluon_pT_Tau","Gluon_pT_Tau",1,1,2);
  book<TH1F>("Up_SS_pT_Tau","Up_SS_pT_Tau",1,1,2);
  book<TH1F>("Down_SS_pT_Tau","Down_SS_pT_Tau",1,1,2);
  book<TH1F>("ChargeFlip_pT_Tau","ChargeFlip_pT_Tau",1,1,2);


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
  const auto & tau = (*event.taus)[0];

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
    hist("tau_type")->Fill(0);
    if(event.taus->size() > 0){
      const auto & tau = (*event.taus)[0];
      hist("pt_real_tau1_binned")->Fill(tau.pt(), weight);
      //hist("pt_tau1_binned")->Fill(tau.pt(), weight);
    }
  }
  else{
    hist("tau_type")->Fill(1);
    if(event.taus->size() > 0){
      const auto & tau = (*event.taus)[0];
      hist("pt_fake_tau1_binned")->Fill(tau.pt(), weight);
      //hist("pt_tau1_binned")->Fill(tau.pt(), weight);
    }
    for(auto genp_new : *event.genparticles){
      double delR = deltaR(tau,genp_new);
      if(delR<0.4){
	//if((genp_new.pdgId()!=21) && fabs(genp_new.pdgId())!=1 && fabs(genp_new.pdgId()) !=2 && fabs(genp_new.pdgId()) !=3 && fabs(genp_new.pdgId())!=4 && fabs(genp_new.pdgId())!=5) cout << genp_new.pdgId() << endl;
	
	if(genp_new.pdgId()==21){
	  
	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_gluon_tau1_binned")->Fill(tau.pt(), weight);
	  }
	  
	  hist("Gluon_pT_Tau")->Fill(1, weight);
	  
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
	if (((genp_new.pdgId() ==2 || genp_new.pdgId() ==4) && tau.charge()==1) || ((genp_new.pdgId() ==-2 || genp_new.pdgId() ==-4) && tau.charge()==-1)){

	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_Up_SS_tau1_binned")->Fill(tau.pt(), weight);
	  }

	  hist("Up_SS_pT_Tau")->Fill(1, weight);

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
                     
	if (((genp_new.pdgId() ==1 || genp_new.pdgId() ==3 || genp_new.pdgId() ==5) && tau.charge()==1) || ((genp_new.pdgId() ==-1 || genp_new.pdgId() ==-3 || genp_new.pdgId() ==-5) && tau.charge()==-1)){

	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_Down_SS_tau1_binned")->Fill(tau.pt(), weight);
	  }

	  hist("Down_SS_pT_Tau")->Fill(1, weight);

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
	if (((genp_new.pdgId() ==1 || genp_new.pdgId() ==3 || genp_new.pdgId() ==5) && tau.charge()==-1) || ((genp_new.pdgId() ==-1 || genp_new.pdgId() ==-3 || genp_new.pdgId() ==-5) && tau.charge()==1) || ((genp_new.pdgId() ==2 || genp_new.pdgId() ==4) && tau.charge()==-1) || ((genp_new.pdgId() ==-2 || genp_new.pdgId() ==-4) && tau.charge()==1)){

	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_ChargeFlip_tau1_binned")->Fill(tau.pt(), weight);
	  }

	  hist("ChargeFlip_pT_Tau")->Fill(1, weight);

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


	if(abs(genp_new.pdgId())==24){
	  
	  if(event.taus->size() > 0){
	    const auto & tau = (*event.taus)[0];
	    hist("pt_fake_WBoson_tau1_binned")->Fill(tau.pt(), weight);
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
