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



SameSignCut::SameSignCut(){
}

 bool SameSignCut::passes(const Event & event)
 {
   for(const auto & muon : *event.muons)
     {
       for(const auto & tau : *event.taus)
	 {
	   if (muon.charge() == tau.charge()) return true;
	 }     
     }
   return false;
 }

EleTauSameSignCut::EleTauSameSignCut(){
}

 bool EleTauSameSignCut::passes(const Event & event)
 {
   for(const auto & electron : *event.electrons)
     {
       for(const auto & tau : *event.taus)
	 {
	   if (electron.charge() == tau.charge()) return true;
	 }     
     }
   return false;
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
  bool pass = false;
  double MET = event.met->pt();
  pass = MET > min_met && (max_met < 0 || MET < max_met);
  return pass;
}



// see https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation53XReReco
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


ElectronIso::ElectronIso(double iso_):iso(iso_){}
bool ElectronIso::operator()(const Electron & electron, const uhh2::Event &) const {
if(electron.relIso()>iso) return false;
return true;
}


