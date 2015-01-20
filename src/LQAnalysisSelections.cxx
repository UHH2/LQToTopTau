#include "UHH2/LQAnalysis/include/LQAnalysisSelections.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ObjectIdUtils.h"



#include <stdexcept>

using namespace uhh2;
using namespace std;
/*
NJetSelection::NJetSelection(int nmin_, int nmax_): nmin(nmin_), nmax(nmax_){}
bool NJetSelection::passes(const Event & event){
    int njets = event.jets->size();
    return njets >= nmin && (nmax < 0 || njets <= nmax);
}
*/
/*
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
*/

METCut::METCut(double min_met_, double max_met_): min_met(min_met_), max_met(max_met_){}

bool METCut::passes(const Event & event)
{
  return event.met->pt() > min_met && event.met->pt() < max_met;
  /*double MET = event.met->pt();
  if( MET < min_met) return false;
  if( MET > max_met) return false;
  return true;*/ 
}



