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
  book<TH1F>("ST", "H_{T}", 50, 0, 5000);
  book<TH1F>("ST_binned", "H_{T}", 8, bins);
  book<TH1F>("ST_testbinned", "H_{T}", 32, 800,4000);

  book<TH1F>("isolation","#tau iso [GeV]",20,0,2);
  book<TH1F>("isolation_sideband","#tau iso [GeV]",15,0,300);
  double isobins[101];
  isobins[0]=0.5;
  isobins[1]=1.5;
  for(int i=2; i<101; i++){
    isobins[i]=i*4;
  }
  book<TH1F>("isolation_fullrange","#tau iso [GeV]",100,isobins);

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
  book<TH1F>("M_btau", "M_{b#tau}", 20, 0,1000);

  book<TH1F>("DeltaR_btau", "DeltaR_{b#tau}", 50, 0,5);
  book<TH1F>("pt_b", "pt bjetM", 50,0,1000);
  book<TH1F>("pt_b_binned", "pt bjetM", 4,bbins);


  book<TH1F>("MT", "M_{T}(mu,missing ET)", 50,0,800);
  //book<TH1F>("Weights", "weights", 2,0.5,2.5);
  book<TH1F>("eta_tilde", "|#tilde{#eta}|", 8,0,2.4);
  book<TH1F>("eta_tilde2", "|#tilde{#eta}|", 24,0,2.4);
  double eta_bins[3]={0,0.9,2.5};
  book<TH1F>("eta_tilde_bin1", "|#tilde{#eta}|", 2,eta_bins);
  book<TH1F>("eta_tilde_bin2", "|#tilde{#eta}|", 2,eta_bins);
  book<TH1F>("muon_type", "0 real muon, 1 fake muon", 2,-0.5,1.5);
  book<TH1F>("electron_type", "0 real electron, 1 fake electrron", 2,-0.5,1.5);
  book<TH1F>("tau_type", "0 real tau, 1 fake tau", 2,-0.5,1.5);

  book<TH1F>("pt_real_tau1_binned","p_{T} real tau 1",4,taubins);
  book<TH1F>("pt_fake_tau1_binned","p_{T} fake tau 1",4,taubins);
  book<TH1F>("pt_tau1_binned","p_{T} tau 1",4,taubins);
  book<TH1F>("pt_tau_binned","p_{T} taus",4,taubins);

  book<TH1F>("M_jet", "M_{Jet}", 100, 0, 2000);
  book<TH1F>("N_subjets", "N_{Subjets} in a Topjet", 11, -0.5, 10.5);
  book<TH1F>("min_mDisubjet", "Min(m_{ij})", 50, 0, 1000);
  book<TH1F>("N_TopTags", "Number of CMSTopTags",6 ,-0.5, 5.5 );

  book<TH1F>("pt_thad", "P_{T}^{top,had} [GeV/c]", 20,0,1200);
  book<TH1F>("pt_tlep", "P_{T}^{top,lep} [GeV/c]", 20,0,1200);
  book<TH1F>("pt_tcom", "P_{T}^{top,combined} [GeV/c]", 20,0,1200);

  book<TH1F>("sign", "#mu#tau sign", 2,-1,1);

  book<TH1F>("N", "counting exp.", 1,0.5,1.5);


  double pttoprebins[5]={0,100,200,300,1200};
  double pttopbins[4]={0,100,200,1200};
  double testbin[5] = {0,120,200,380,1200};

  book<TH1F>("M_tophad_own", "M^{top,had} [GeV/c^{2}]", 50, 0, 500 ) ;
  book<TH1F>("Pt_tophad_own", "P_{T}^{top,had} [GeV/c]", 60, 0, 1200 ) ;
  book<TH1F>("Pt_tophad_own_rebin", "P_{T}^{top,had} [GeV/c]", 12, 0, 1200 ) ;
  book<TH1F>("Pt_tophad_own_binned", "P_{T}^{top,had} [GeV/c]", 4, pttoprebins ) ;
  book<TH1F>("Pt_tophad_own_binned2", "P_{T}^{top,had} [GeV/c]", 3, pttopbins ) ;
  book<TH1F>("Pt_tophad_own_test", "P_{T}^{top,had} [GeV/c]", 4, testbin ) ;
  book <TH1F>("Chi2", "#chi^{2}", 100, 0, 200);
  book <TH1F>("top_recjets", "number of jets used for top hypothesis", 6, 0.5, 6.5);

  /*
  book<TH2F>("pt_tau1_vs_ST","pt tau 1 vs ST", 50, 0, 800 ,100, 0, 2500);
  book<TH2F>("pt_tau1_vs_MET","pt tau 1 vs MET", 50, 0, 800 ,50, 0, 1000);
  */

  //h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>("HighMassTTbarReconstruction");
  //m_discriminator_name = "Chi2";


  h_hadr_hyps = ctx.get_handle<std::vector<TTbarFullhadRecoHypothesis>>("HighMassHadronicTTbarFullhadReco");

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


  for(const auto & tau : *event.taus){
    hist("isolation")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
    hist("isolation_sideband")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
    hist("isolation_fullrange")->Fill(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
    hist("pt_tau_binned")->Fill(tau.pt(),weight);
  }
  


  
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

  
  if(muons->size()>0){
    if(taus->size()>0){
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
    hist("MT")->Fill(sqrt(2*muon.pt()*event.met->pt()* (1-cos(event.met->phi()-muon.phi())) ), weight);
  }

  const auto jets = event.jets;
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
	}
      }
    }
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
    }
  }

  /*
  if(event.taus->size() > 0){
    const auto & tau = (*event.taus)[0];
    hist("pt_tau1_binned")->Fill(tau.pt(), weight);
    ((TH2F*)hist("pt_tau1_vs_ST"))->Fill(tau.pt(),ht+ht_lep+met, weight);
    ((TH2F*)hist("pt_tau1_vs_MET"))->Fill(tau.pt(),met, weight);
  }    
  */
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
	hist("muon_type")->Fill(0);
      }
      else{
	hist("muon_type")->Fill(1);
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
	hist("electron_type")->Fill(0);
      }
      else{
	hist("electron_type")->Fill(1);
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
