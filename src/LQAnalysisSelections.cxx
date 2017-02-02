#include "UHH2/LQAnalysis/include/LQAnalysisSelections.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ObjectIdUtils.h"



#include <stdexcept>

using namespace uhh2;
using namespace std;


InvMass2MuVeto::InvMass2MuVeto(double m_min_, double m_max_):m_min(m_min_), m_max(m_max_){}
bool InvMass2MuVeto::passes(const Event & event){
  const auto muons = event.muons;
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
	  if(InvMass > m_min && InvMass < m_max) return false;
	}
    }
  return true;
}

OppositeSignCut::OppositeSignCut(){
}

bool OppositeSignCut::passes(const Event & event)
{
  for(const auto & muon : *event.muons)
    {
      for(const auto & tau : *event.taus)
	{
	  if (muon.charge() == tau.charge()) return false;
	}     
    }
  return true;
}

SameSignCut::SameSignCut(){
}

bool SameSignCut::passes(const Event & event)
{
  for(const auto & muon : *event.muons)
    {
      for(const auto & tau : *event.taus)
	{
	  if (muon.charge() != tau.charge()) return false;
	}     
    }
  return true;
}



SameSignCutLeadingLep::SameSignCutLeadingLep(){
}
bool SameSignCutLeadingLep::passes(const Event & event)
{
  if(event.muons->size() > 0 && event.taus->size() > 0){
    const auto & muon = (*event.muons)[0];
    const auto & tau = (*event.taus)[0];
    if (muon.charge() == tau.charge()) return true;
  }
  return false;
}


EleTauSameSignCut::EleTauSameSignCut(){
}
bool EleTauSameSignCut::passes(const Event & event)
{
  if(event.electrons->size() > 0 && event.taus->size() > 0){
    const auto & ele = (*event.electrons)[0];
    const auto & tau = (*event.taus)[0];
    if (ele.charge() == tau.charge()) return true;
  }
  return false;
}

JetTauCleaning::JetTauCleaning(){
}
bool JetTauCleaning::passes(const Event & event)
{
  /*
  for(unsigned int i=0; i<event.jets->size(); i++){
    Jet jet = event.jets->at(i);
    for(const auto & tau : *event.taus){
      if(event.jets->size()>0){
	if(deltaR(tau,jet)<0.4){
	  event.jets->erase(event.jets->begin()+i);
	  i--;
	}
      }
    }
  }*/
  for(const auto & tau : *event.taus){
    for(unsigned int i=0; i<event.jets->size(); i++){
      if(event.jets->size()>0){
	Jet jet = event.jets->at(i);
	if(deltaR(tau,jet)<0.4){
	  event.jets->erase(event.jets->begin()+i);
	  i--;
	}
      }
    }
  }
  return true;
}

GetFakeTaus::GetFakeTaus(){
}
bool GetFakeTaus::passes(const Event & event)
{
  for(unsigned int i =0; i<event.taus->size(); ++i){
    Tau tau = event.taus->at(i);
    for(auto genp : *event.genparticles){
      double dR = deltaR(tau,genp);
      if(dR<0.4 && abs(genp.pdgId())==15){
	return false;
      }
    }
  }
  return true;
}

GetRealTaus::GetRealTaus(){
}
bool GetRealTaus::passes(const Event & event)
{
  double dR = 1000;
  for(const auto & tau : *event.taus){
    for(auto genp : *event.genparticles){
      if(abs(genp.pdgId())==15){
	double tmp = deltaR(tau,genp);
	if(tmp<dR){
	  dR = tmp;
	}
      }
    }
  }
  if(dR<0.4){
  }
  else{
    return false;
  }
  return true;
}

IsoNonIso::IsoNonIso(){
}
bool IsoNonIso::passes(const Event & event){
  double N_iso = 0;
  double N_noniso = 0;
  for(const auto & tau : *event.taus){
    if(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()<1.5){
      N_iso += 1;
    }
    if(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()>=1.5){
      N_noniso += 1;
    }
  }
  if(N_iso>0 && N_noniso>0) return false;
  return true;
  }



MbtauSelection::MbtauSelection(double minMbtau, double maxMbtau): minMbtau_(minMbtau),maxMbtau_(maxMbtau) {}
bool MbtauSelection::passes(const Event & event){
  const auto jets = event.jets;
  vector<Jet> bjets;
  for (unsigned int i =0; i<jets->size(); ++i) {
    if(jets->at(i).btag_combinedSecondaryVertex()>0.679) {
      bjets.push_back(jets->at(i));
    }
  }
  //bool pass = false;
  const auto taus = event.taus;
  for (unsigned int i=0; i<=(*taus).size(); ++i){
    if((*taus).size()>i){
      Tau tau = (*taus)[i];
      TLorentzVector Tau;
      Tau.SetPtEtaPhiE(tau.pt() ,tau.eta() ,tau.phi() ,tau.energy() );
      for (unsigned int i =0; i<=bjets.size(); ++i) {
	if (bjets.size()> i) {
	  Jet bjet = bjets[i];
	  TLorentzVector BJet;
	  BJet.SetPtEtaPhiE(bjet.pt() ,bjet.eta() ,bjet.phi() ,bjet.energy() );
	  double invmass=(Tau+BJet).M();
	  //return pass = invmass > minMbtau_ && (maxMbtau_ < 0 || invmass < maxMbtau_);
	  if(invmass < minMbtau_ && (maxMbtau_ < 0 || invmass > maxMbtau_)) return false;
	}
      }
    }
  }
  //return false;
  return true;
}


NJetCut::NJetCut(int nmin_, int nmax_, double ptmin_, double etamax_): nmin(nmin_), nmax(nmax_), ptmin(ptmin_), etamax(etamax_){}
bool NJetCut::passes(const Event & event){
  int nparticle=0;
  for(auto & jet : *event.jets) {
    if (jet.pt() > ptmin && fabs(jet.eta()<etamax)) nparticle++;
  }
  return nparticle >= nmin;
}



METCut::METCut(double min_met_, double max_met_): min_met(min_met_), max_met(max_met_){}
bool METCut::passes(const Event & event){
  double MET = event.met->pt();
  bool pass = false;
  pass = MET > min_met && (max_met < 0 || MET < max_met);
  return pass;
}

MTCut::MTCut(double min_mt_, double max_mt_): min_mt(min_mt_), max_mt(max_mt_){}
bool MTCut::passes(const Event & event){
  double MT_muMET =0.;
  for(const auto & muon : *event.muons){
     MT_muMET = sqrt(2*muon.pt()*event.met->pt()* (1-cos(deltaPhi(*event.met,muon)) ));
    break;
  }
  bool pass = false;
  pass = MT_muMET > min_mt && (max_mt < 0 || MT_muMET < max_mt);
  return pass;
}


HtSelection::HtSelection(double ht_min_, double ht_max_):ht_min(ht_min_), ht_max(ht_max_){}
bool HtSelection::passes(const Event & event){
  auto met = event.met->pt();

  double ht = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;
  for(const auto & jet : *event.jets){
    ht_jets += jet.pt();
  }
  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }
  for(const auto & tau : *event.taus){
    ht_lep += tau.pt();
  }
  ht = ht_lep + ht_jets + met;

  bool pass = false;
  pass = ht > ht_min && (ht_max < 0 || ht < ht_max);
  return pass;
}

PtLeadingJetSelection::PtLeadingJetSelection(double pt_min_, double pt_max_):pt_min(pt_min_), pt_max(pt_max_){}
bool PtLeadingJetSelection::passes(const Event & event){

  bool pass = true;
  double pt_leadingjet = event.jets->at(0).pt();
  pass = pt_leadingjet >= pt_min && (pt_leadingjet <= pt_max || pt_max < 0);
  return pass;
}

// see https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation53XReReco
/*
float btagging::csv_threshold(const csv_wp & wp){
    using namespace btagging;
    switch(wp){
        case csv_wp::loose: return 0.244f;
        case csv_wp::medium: return 0.679f;
        case csv_wp::tight: return 0.898f;
    }
    // This should never happen; even if, the coompiler should warn in the switch.
    // But to avoid a compiler warning that no value is returned, include this line:
    throw invalid_argument("unknown working point given to btagging::csv_threshold");
}

NBTagSelection::NBTagSelection(int nmin_, int nmax_, btagging::csv_wp wp): nmin(nmin_), nmax(nmax_), min_csv(btagging::csv_threshold(wp)){}

bool NBTagSelection::passes(const Event & event){
    int nbtag = 0;
    for(const Jet & j : *event.jets){
        if(j.btag_combinedSecondaryVertex() >= min_csv) ++nbtag;
    }
    return nbtag >= nmin && (nmax < 0 || nbtag <= nmax);
}
*/

ElectronIso::ElectronIso(double iso_):iso(iso_){}
bool ElectronIso::operator()(const Electron & electron, const uhh2::Event &) const {
if(electron.relIso()>iso) return false;
return true;
}

LQMatch::LQMatch() {}    
bool LQMatch::passes(const Event & event){
  int n_ele = 0, n_mu = 0, n_tau = 0, sum = 0, n_j1=0, n_j2=0, n_j3=0, n_j4=0, n_j5=0, n_j6=0;
  assert(event.genparticles);

  for(const auto & gp : *event.genparticles){
    if(gp.pdgId() == 6){ 
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
      if(abs(wd1->pdgId())==13 || abs(wd2->pdgId())==13) n_mu++;
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
    if(gp.pdgId() == -6){ 
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
      if(abs(wd1->pdgId())==13 || abs(wd2->pdgId())==13) n_mu++;
      if(abs(wd1->pdgId())==15 || abs(wd2->pdgId())==15) n_tau++;
      for(unsigned int i =0; i<event.jets->size(); ++i){
	Jet jet = event.jets->at(i);
	double dR1 = deltaR(jet,*topd2);
	double dR2 = deltaR(jet,*wd1);
	double dR3 = deltaR(jet,*wd2);
	if(dR1<0.4 && abs(topd2->pdgId())<7) n_j4++;
	if(dR2<0.4 && abs(wd1->pdgId())<7) n_j5++;
	if(dR3<0.4 && abs(wd2->pdgId())<7) n_j6++;
      }
    }
  }

  sum = n_ele + n_mu + n_tau;
  if(sum>1) return false;
  if(! ((n_j1>0 && n_j2>0 && n_j3>0) || (n_j4>0 && n_j5>0 && n_j6>0)) ) return false;

  /*
  LorentzVector bv4;
  LorentzVector q1v4;
  LorentzVector q2v4;
  float drminb=99, drminq1=99, drminq2=99;
  int n1=0, n2=0;

  if(n_j1>0 && n_j2>0 && n_j3>0){
    for(const auto & gp : *event.genparticles){
      if(gp.pdgId() == 6){ 
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
	if(abs(wd1->pdgId())==11 || abs(wd2->pdgId())==11) n1++;
	if(abs(wd1->pdgId())==13 || abs(wd2->pdgId())==13) n1++;
	if(abs(wd1->pdgId())==15 || abs(wd2->pdgId())==15) n1++;
	if(n1==0){
	  for(unsigned int i =0; i<event.jets->size(); ++i){
	    Jet jet = event.jets->at(i);
	    double dR1 = deltaR(jet,*topd2);
	    double dR2 = deltaR(jet,*wd1);
	    double dR3 = deltaR(jet,*wd2);
	    if(dR1<0.4 && abs(topd2->pdgId())<7 && dR1<drminb){
	      bv4=jet.v4();
	      drminb=dR1;
	    }
	    if(dR2<0.4 && abs(wd1->pdgId())<7 && dR2<drminq1){
	      q1v4=jet.v4();
	      drminq1=dR2;
	    }
	    if(dR3<0.4 && abs(wd2->pdgId())<7 && dR3<drminq2){
	      q2v4=jet.v4();
	      drminq2=dR3;
	    }
	  }
	}
      }
    }
    cout << "1" << endl;
    if(bv4==q1v4 && bv4!=q2v4) cout << "matchable mass: " <<  (bv4+q2v4).M() << endl;
    cout << "2" << endl;
    if(bv4==q2v4 && bv4!=q1v4) cout << "matchable mass: " <<  (bv4+q1v4).M() << endl;
    cout << "3" << endl;
    if(bv4==q1v4 && bv4==q2v4) cout << "matchable mass: " <<  bv4.M() << endl;
    cout << "4" << endl;
    if(bv4!=q1v4 && bv4!=q2v4 && q1v4!=q2v4){ cout << "matchable mass: " << (bv4+q1v4+q2v4).M() << endl;}
  }

  drminb=99;
  drminq1=99;
  drminq2=99;
  LorentzVector abv4;
  LorentzVector aq1v4;
  LorentzVector aq2v4;
  if(n_j4>0 && n_j5>0 && n_j6>0){
    for(const auto & gp : *event.genparticles){
      if(gp.pdgId() == 6){ 
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
	if(abs(wd1->pdgId())==11 || abs(wd2->pdgId())==11) n2++;
	if(abs(wd1->pdgId())==13 || abs(wd2->pdgId())==13) n2++;
	if(abs(wd1->pdgId())==15 || abs(wd2->pdgId())==15) n2++;
	if(n2==0){
	  for(unsigned int i =0; i<event.jets->size(); ++i){
	    Jet jet = event.jets->at(i);
	    double dR1 = deltaR(jet,*topd2);
	    double dR2 = deltaR(jet,*wd1);
	    double dR3 = deltaR(jet,*wd2);
	    if(dR1<0.4 && abs(topd2->pdgId())<7 && dR1<drminb){
	      abv4=jet.v4();
	      drminb=dR1;
	    }
	    if(dR2<0.4 && abs(wd1->pdgId())<7 && dR2<drminq1){
	      aq1v4=jet.v4();
	      drminq1=dR2;
	    }
	    if(dR3<0.4 && abs(wd2->pdgId())<7 && dR3<drminq2){
	      aq2v4=jet.v4();
	      drminq2=dR3;
	    }
	  }
	}
      }
    }
    cout << "5" << endl;
    if(abv4==aq1v4 && abv4!=aq2v4) cout << "matchable mass: " <<  (abv4+aq2v4).M() << endl;
    cout << "6" << endl;
    if(abv4==aq2v4 && abv4!=aq1v4) cout << "matchable mass: " <<  (abv4+aq1v4).M() << endl;
    cout << "7" << endl;
    if(abv4==aq1v4 && abv4==aq2v4) cout << "matchable mass: " <<  abv4.M() << endl;
    cout << "8" << endl;
    if(abv4!=aq1v4 && abv4!=aq2v4 && aq1v4!=aq2v4){ cout << "matchable mass: " << (abv4+aq1v4+aq2v4).M() << endl;}
  }
  */

  return true;
}



