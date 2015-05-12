#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TauHists.h"
#include "UHH2/common/include/TauIds.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/LQAnalysis/include/LQAnalysisSelections.h"
#include "UHH2/LQAnalysis/include/LQAnalysisHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/EventVariables.h"

using namespace std;
using namespace uhh2;

/** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
 * 
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class LQAnalysisOldMuModule: public AnalysisModule {
public:
    
    explicit LQAnalysisOldMuModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
    
  std::unique_ptr<MuonCleaner> muoncleaner_iso;
  std::unique_ptr<ElectronCleaner> electroncleaner_iso;
  std::unique_ptr<JetCleaner> jetcleaner;
  std::unique_ptr<MuonIDKinematic> muonidkinematic;
  std::unique_ptr<MuonCleaner> muoncleaner;
  std::unique_ptr<TauCleaner> taucleaner;
  std::unique_ptr<ElectronCleaner> electroncleaner;
  
  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  std::unique_ptr<Selection> njet_sel, TwoBTagL, BTagM, BTagT, ntau_sel, ele_sel, muon_sel, isomuon_sel, leadingjet_sel, leadingjet300_sel, secondjet100_sel;
  std::vector<std::unique_ptr<Selection> > fullhad_sel;

  std::vector<std::unique_ptr<AnalysisModule>> pre_modules;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> lq_PreSelection, electron_PreSelection, muon_PreSelection, jet_PreSelection, tau_PreSelection, event_PreSelection;
  std::unique_ptr<Hists> h_lq_LeadingJet150, h_tau_LeadingJet150, h_mu_LeadingJet150, h_ele_LeadingJet150, h_jet_LeadingJet150, h_event_LeadingJet150;
  std::unique_ptr<Hists> h_lq_threejets, h_tau_threejets, h_mu_threejets, h_ele_threejets, h_jet_threejets, h_event_threejets;
  std::unique_ptr<Hists> h_lq_st1000, h_tau_st1000, h_mu_st1000, h_ele_st1000, h_jet_st1000, h_event_st1000;
  std::unique_ptr<Hists> h_lq_met50, h_tau_met50, h_mu_met50, h_ele_met50, h_jet_met50, h_event_met50;
  std::unique_ptr<Hists> h_lq_mediumtau, h_tau_mediumtau, h_mu_mediumtau, h_ele_mediumtau, h_jet_mediumtau, h_event_mediumtau;
  std::unique_ptr<Hists> h_lq_MuIso, h_tau_MuIso, h_mu_MuIso, h_ele_MuIso, h_jet_MuIso, h_event_MuIso;
  std::unique_ptr<Hists> h_lq_Mareike, h_tau_Mareike, h_mu_Mareike, h_ele_Mareike, h_jet_Mareike, h_event_Mareike;
  std::unique_ptr<Hists> h_ele_full, h_tau_full, h_event_full, h_lq_full, h_jet_full, h_muon_full;


  JetId BTagLoose, BTagMedium, BTagTight;
  MuonId MuIso;
  ElectronId EleIso;

};


LQAnalysisOldMuModule::LQAnalysisOldMuModule(Context & ctx){
    // In the constructor, the typical tasks are to create
    // other modules like cleaners (1), selections (2) and Hist classes (3).
    // But you can do more and e.g. access the configuration, as shown below.
    
    cout << "Hello World from LQAnalysisOldMuModule!" << endl;
    
    // If needed, access the configuration of the module here, e.g.:
    string testvalue = ctx.get("TestKey", "<not set>");
    cout << "TestKey in the configuration was: " << testvalue << endl;
    
    // If running in SFrame, the keys "dataset_version", "dataset_type" and "dataset_lumi"
    // are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
    
    // 1. setup other modules.
    jetcleaner.reset(new JetCleaner(30.0, 2.4));
    muonidkinematic.reset(new MuonIDKinematic(30.0,3.0));
    muoncleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), MuonIDKinematic(30.0, 2.1))));
    electroncleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_PHYS14_25ns_medium, PtEtaCut(20.0, 2.5))));
    taucleaner.reset(new TauCleaner(AndId<Tau>(TauIDMedium(), PtEtaCut(20.0, 2.1))));
    BTagLoose = CSVBTag(CSVBTag::WP_LOOSE);
    BTagMedium = CSVBTag(CSVBTag::WP_MEDIUM);
    BTagTight = CSVBTag(CSVBTag::WP_TIGHT);
    MuIso = MuonIso(0.12);
    EleIso = ElectronIso(0.12);
    muoncleaner_iso.reset(new MuonCleaner(AndId<Muon>(MuIso, MuonIDKinematic(30.0, 2.4))));
    electroncleaner_iso.reset(new ElectronCleaner(AndId<Electron>(EleIso, PtEtaCut(30.0, 2.5))));

    pre_modules.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));

    // 2. set up selections:
    njet_sel.reset(new NJetSelection(3,-1));
    leadingjet_sel.reset(new NJetCut(1,-1,150,5.0));
    leadingjet300_sel.reset(new NJetCut(1,-1,300,5.0));
    secondjet100_sel.reset(new NJetCut(2,-1,100,5.0));
    ntau_sel.reset(new NTauSelection(1,-1));
    TwoBTagL.reset(new NJetSelection(2,999,BTagLoose));
    BTagM.reset(new NJetSelection(1,999,BTagMedium));
    BTagT.reset(new NJetSelection(1,999,BTagTight));
    ele_sel.reset(new NElectronSelection(0,0));
    muon_sel.reset(new NMuonSelection(1,-1));
    isomuon_sel.reset(new NMuonSelection(1,-1,MuIso));

    int n_cuts = 4;
    fullhad_sel.resize(n_cuts);
    fullhad_sel[0].reset(new NJetSelection(3,-1));
    fullhad_sel[1].reset(new NTauSelection(1));
    fullhad_sel[2].reset(new NJetSelection(1,999,BTagMedium));
    fullhad_sel[3].reset(new NJetCut(2,-1,50,3.0));    
    //fullhad_sel[4].reset(new METCut(185,-1));

    // 3. Set up Hists classes:
    electron_PreSelection.reset(new ElectronHists(ctx, "LQMod_Electrons_PreSel"));
    muon_PreSelection.reset(new MuonHists(ctx, "LQMod_Muons_PreSel"));
    tau_PreSelection.reset(new TauHists(ctx, "LQMod_Taus_PreSel"));
    jet_PreSelection.reset(new JetHists(ctx, "LQMod_Jets_PreSel"));
    event_PreSelection.reset(new EventHists(ctx, "LQMod_Events_PreSel"));
    lq_PreSelection.reset(new LQAnalysisHists(ctx, "LQMod_LQ_PreSel"));

    h_lq_LeadingJet150.reset(new LQAnalysisHists(ctx, "LQMod_LQ_LeadingJet150"));
    h_tau_LeadingJet150.reset(new TauHists(ctx, "LQMod_Taus_LeadingJet150"));
    h_mu_LeadingJet150.reset(new MuonHists(ctx, "LQMod_Muons_LeadingJet150"));
    h_ele_LeadingJet150.reset(new ElectronHists(ctx, "LQMod_Electrons_LeadingJet150"));
    h_jet_LeadingJet150.reset(new JetHists(ctx, "LQMod_Jets_LeadingJet150"));
    h_event_LeadingJet150.reset(new EventHists(ctx, "LQMod_Events_LeadingJet150"));

    h_lq_MuIso.reset(new LQAnalysisHists(ctx, "LQMod_LQ_MuIso"));
    h_tau_MuIso.reset(new TauHists(ctx, "LQMod_Taus_MuIso"));
    h_mu_MuIso.reset(new MuonHists(ctx, "LQMod_Muons_MuIso"));
    h_ele_MuIso.reset(new ElectronHists(ctx, "LQMod_Electrons_MuIso"));
    h_jet_MuIso.reset(new JetHists(ctx, "LQMod_Jets_MuIso"));
    h_event_MuIso.reset(new EventHists(ctx, "LQMod_Events_MuIso"));

    h_lq_threejets.reset(new LQAnalysisHists(ctx, "LQMod_LQ_threejets"));
    h_tau_threejets.reset(new TauHists(ctx, "LQMod_Taus_threejets"));
    h_mu_threejets.reset(new MuonHists(ctx, "LQMod_Muons_threejets"));
    h_ele_threejets.reset(new ElectronHists(ctx, "LQMod_Electrons_threejets"));
    h_jet_threejets.reset(new JetHists(ctx, "LQMod_Jets_threejets"));
    h_event_threejets.reset(new EventHists(ctx, "LQMod_Events_threejets"));

    h_lq_st1000.reset(new LQAnalysisHists(ctx, "LQMod_LQ_st1000"));
    h_tau_st1000.reset(new TauHists(ctx, "LQMod_Taus_st1000"));
    h_mu_st1000.reset(new MuonHists(ctx, "LQMod_Muons_st1000"));
    h_ele_st1000.reset(new ElectronHists(ctx, "LQMod_Electrons_st1000"));
    h_jet_st1000.reset(new JetHists(ctx, "LQMod_Jets_st1000"));
    h_event_st1000.reset(new EventHists(ctx, "LQMod_Events_st1000"));

    h_lq_met50.reset(new LQAnalysisHists(ctx, "LQMod_LQ_met50"));
    h_tau_met50.reset(new TauHists(ctx, "LQMod_Taus_met50"));
    h_mu_met50.reset(new MuonHists(ctx, "LQMod_Muons_met50"));
    h_ele_met50.reset(new ElectronHists(ctx, "LQMod_Electrons_met50"));
    h_jet_met50.reset(new JetHists(ctx, "LQMod_Jets_met50"));
    h_event_met50.reset(new EventHists(ctx, "LQMod_Events_met50"));

    h_lq_mediumtau.reset(new LQAnalysisHists(ctx, "LQMod_LQ_mediumtau"));
    h_tau_mediumtau.reset(new TauHists(ctx, "LQMod_Taus_mediumtau"));
    h_mu_mediumtau.reset(new MuonHists(ctx, "LQMod_Muons_mediumtau"));
    h_ele_mediumtau.reset(new ElectronHists(ctx, "LQMod_Electrons_mediumtau"));
    h_jet_mediumtau.reset(new JetHists(ctx, "LQMod_Jets_mediumtau"));
    h_event_mediumtau.reset(new EventHists(ctx, "LQMod_Events_mediumtau"));

    h_lq_Mareike.reset(new LQAnalysisHists(ctx, "LQMod_LQ_Mareike"));
    h_tau_Mareike.reset(new TauHists(ctx, "LQMod_Taus_Mareike"));
    h_mu_Mareike.reset(new MuonHists(ctx, "LQMod_Muons_Mareike"));
    h_ele_Mareike.reset(new ElectronHists(ctx, "LQMod_Electrons_Mareike"));
    h_jet_Mareike.reset(new JetHists(ctx, "LQMod_Jets_Mareike"));
    h_event_Mareike.reset(new EventHists(ctx, "LQMod_Events_Mareike"));

    h_ele_full.reset(new ElectronHists(ctx, "ele_fullselection"));
    h_muon_full.reset(new MuonHists(ctx, "muon_fullselection"));
    h_jet_full.reset(new JetHists(ctx, "jet_fullselection"));
    h_event_full.reset(new EventHists(ctx, "event_fullselection"));
    h_tau_full.reset(new TauHists(ctx, "tau_fullselection"));
    h_lq_full.reset(new LQAnalysisHists(ctx, "lq_fullselection"));




}


bool LQAnalysisOldMuModule::process(Event & event) {
    // This is the main procedure, called for each event. Typically,
    // do some pre-processing by calling the modules' process method
    // of the modules constructed in the constructor (1).
    // Then, test whether the event passes some selection and -- if yes --
    // use it to fill the histograms (2).
    // Finally, decide whether or not to keep the event in the output (3);
    // this is controlled by the return value of this method: If it
    // returns true, the event is kept; if it returns false, the event
    // is thrown away.
    
    //cout << "LQAnalysisOldMuModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
    
    // 1. run all modules; here: only jet cleaning.
    /*jetcleaner->process(event);
    muoncleaner->process(event);
    electroncleaner->process(event);
    taucleaner->process(event);*/

    for (auto & mod : pre_modules) {
      mod->process(event);
    }
    
    
    /*
    const auto jets = event.jets;
    if(jets->size() > 0){
      const auto & jet = (*jets)[0];
      if(jet.pt()<150) return false;
    }
    if(jets->size() > 1){
      const auto & jet = (*jets)[1];
      if(jet.pt()<50) return false;
    }
    if(jets->size() > 2){
      const auto & jet = (*jets)[2];
      if(jet.pt()<30) return false;
    }
    */

    // 2. test selections and fill histograms

    muoncleaner_iso->process(event);
    electroncleaner_iso->process(event);
    if(!isomuon_sel->passes(event)) return false;


    electron_PreSelection->fill(event);
    muon_PreSelection->fill(event);
    tau_PreSelection->fill(event);
    jet_PreSelection->fill(event);
    event_PreSelection->fill(event);
    lq_PreSelection->fill(event);

    // define met and st
    auto met = event.met->pt();
    double ht = 0.0;
    for(const auto & jet : *event.jets){
      ht += jet.pt();
    }
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
    double st=0.0;
    st = ht + ht_lep + met;

    /*
    const auto jets = event.jets;
    if(jets->size() > 0){
      const auto & jet = (*jets)[0];
      if(jet.pt()<150) return false;
    }
    */
    // muon isolation
    /*for (const auto & muon : *event.muons){
      if(muon.relIso()>0.12) return false;
      }*/

    //if (met<50) return false;
    //if (st<1000) return false;


    if(isomuon_sel->passes(event)){
      h_lq_MuIso->fill(event);
      h_tau_MuIso->fill(event);
      h_mu_MuIso->fill(event);
      h_ele_MuIso->fill(event);
      h_jet_MuIso->fill(event);
      h_event_MuIso->fill(event);
    }

    //if(!isomuon_sel->passes(event)) return false;

    const auto jets = event.jets;
    if(jets->size() > 0){
      const auto & jet = (*jets)[0];
      if(jet.pt()>150){
      h_lq_LeadingJet150->fill(event);
      h_tau_LeadingJet150->fill(event);
      h_mu_LeadingJet150->fill(event);
      h_ele_LeadingJet150->fill(event);
      h_jet_LeadingJet150->fill(event);
      h_event_LeadingJet150->fill(event);
      }
    }

    if (njet_sel->passes(event)){
      h_lq_threejets->fill(event);
      h_tau_threejets->fill(event);
      h_mu_threejets->fill(event);
      h_ele_threejets->fill(event);
      h_jet_threejets->fill(event);
      h_event_threejets->fill(event);
    }

    if(st>1000){
      h_lq_st1000->fill(event);
      h_tau_st1000->fill(event);
      h_mu_st1000->fill(event);
      h_ele_st1000->fill(event);
      h_jet_st1000->fill(event);
      h_event_st1000->fill(event);
    }
    
    if(met>50){
      h_lq_met50->fill(event);
      h_tau_met50->fill(event);
      h_mu_met50->fill(event);
      h_ele_met50->fill(event);
      h_jet_met50->fill(event);
      h_event_met50->fill(event);
    }

    taucleaner->process(event);
    for(unsigned int i=0; i<event.jets->size(); i++){
      Jet jet = event.jets->at(i);
      for(const auto & tau : *event.taus){
	if(deltaR(tau,jet)<0.4){
	  event.jets->erase(event.jets->begin()+i);
	  i--;
	}
      }
    }
    if(!ntau_sel->passes(event)) return false;
    h_lq_mediumtau->fill(event);
    h_tau_mediumtau->fill(event);
    h_mu_mediumtau->fill(event);
    h_ele_mediumtau->fill(event);
    h_jet_mediumtau->fill(event);
    h_event_mediumtau->fill(event);
    

    if(jets->size() > 0){
      const auto & jet = (*jets)[0];
      if(jet.pt()<150) return false;
    }
    if (!njet_sel->passes(event)) return false;
    if(st<1000) return false;
    if(met<50) return false;

    //if(!isomuon_sel->passes(event)) return false;
    h_lq_Mareike->fill(event);
    h_tau_Mareike->fill(event);
    h_mu_Mareike->fill(event);
    h_ele_Mareike->fill(event);
    h_jet_Mareike->fill(event);
    h_event_Mareike->fill(event);



    bool complete_selection = true;
    std::vector<bool> v_accept(fullhad_sel.size());
    for (unsigned i=0; i<fullhad_sel.size(); ++i) {
      bool accept = fullhad_sel[i]->passes(event);
      v_accept[i] = accept;
      if (!accept) {
	complete_selection = false;
      }
    }
    if (complete_selection) {
      h_ele_full->fill(event);
      h_muon_full->fill(event);
      h_jet_full->fill(event);
      h_tau_full->fill(event);
      h_event_full->fill(event);
      h_lq_full->fill(event);
    }
    
    //if (!njet_sel->passes(event)) return false;


    //if (!BTagM->passes(event)) return false;





    
    // 3. decide whether or not to keep the current event in the output:
    return complete_selection;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisOldMuModule)
