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
class LQAnalysisMuModule: public AnalysisModule {
public:
    
    explicit LQAnalysisMuModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
    
    
  std::unique_ptr<JetCleaner> jetcleaner;
  std::unique_ptr<MuonIDKinematic> muonidkinematic;
  std::unique_ptr<MuonCleaner> muoncleaner;
  std::unique_ptr<TauCleaner> taucleaner;
  std::unique_ptr<ElectronCleaner> electroncleaner;
  
  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  std::unique_ptr<Selection> fivejet_sel, BTagM, BTagT, ntau_sel;
  std::vector<std::unique_ptr<Selection> > fullhad_sel;

  std::vector<std::unique_ptr<AnalysisModule>> pre_modules;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> lq_PreSelection, electron_PreSelection, muon_PreSelection, jet_PreSelection, tau_PreSelection, event_PreSelection;
  std::unique_ptr<Hists> h_ele_full, h_tau_full, h_event_full, h_lq_full, h_jet_full, h_muon_full;
  std::unique_ptr<Hists> ele_HT1000, muon_HT1000, tau_HT1000, event_HT1000, jet_HT1000, lq_HT1000;
  std::unique_ptr<Hists> ele_HT1300, muon_HT1300, tau_HT1300, event_HT1300, jet_HT1300, lq_HT1300;
  std::unique_ptr<Hists> ele_HT1800, muon_HT1800, tau_HT1800, event_HT1800, jet_HT1800, lq_HT1800;

  JetId BTagMedium, BTagTight;


};


LQAnalysisMuModule::LQAnalysisMuModule(Context & ctx){
    // In the constructor, the typical tasks are to create
    // other modules like cleaners (1), selections (2) and Hist classes (3).
    // But you can do more and e.g. access the configuration, as shown below.
    
    cout << "Hello World from LQAnalysisMuModule!" << endl;
    
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
    muoncleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), MuonIDKinematic(20.0, 2.1))));
    electroncleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_CSA14_50ns_medium, PtEtaCut(20.0, 2.5))));
    taucleaner.reset(new TauCleaner(AndId<Tau>(TauIDMedium(), PtEtaCut(20.0, 2.1))));
    BTagMedium = CSVBTag(CSVBTag::WP_MEDIUM);
    BTagTight = CSVBTag(CSVBTag::WP_TIGHT);

    pre_modules.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));

    // 2. set up selections:
    fivejet_sel.reset(new NJetSelection(3,-1));
    ntau_sel.reset(new NTauSelection(1,-1));
    BTagM.reset(new NJetSelection(1,999,BTagMedium));
    BTagT.reset(new NJetSelection(1,999,BTagTight));

    int n_cuts = 4;
    fullhad_sel.resize(n_cuts);
    fullhad_sel[0].reset(new NJetSelection(5,-1));
    fullhad_sel[1].reset(new NTauSelection(1));
    fullhad_sel[2].reset(new NJetSelection(1,999,BTagMedium));
    fullhad_sel[3].reset(new NJetCut(2,-1,50,3.0));    
    //fullhad_sel[4].reset(new METCut(185,-1));

    // 3. Set up Hists classes:
    electron_PreSelection.reset(new ElectronHists(ctx, "Ele_PreSel"));
    muon_PreSelection.reset(new MuonHists(ctx, "Mu_PreSel"));
    tau_PreSelection.reset(new TauHists(ctx, "Tau_PreSel"));
    jet_PreSelection.reset(new JetHists(ctx, "Jet_PreSel"));
    event_PreSelection.reset(new EventHists(ctx, "Event_PreSel"));
    lq_PreSelection.reset(new LQAnalysisHists(ctx, "LQhists_PreSel"));

    h_ele_full.reset(new ElectronHists(ctx, "ele_fullselection"));
    h_muon_full.reset(new MuonHists(ctx, "muon_fullselection"));
    h_jet_full.reset(new JetHists(ctx, "jet_fullselection"));
    h_event_full.reset(new EventHists(ctx, "event_fullselection"));
    h_tau_full.reset(new TauHists(ctx, "tau_fullselection"));
    h_lq_full.reset(new LQAnalysisHists(ctx, "lq_fullselection"));

    ele_HT1000.reset(new ElectronHists(ctx, "ele_HT1000"));
    muon_HT1000.reset(new MuonHists(ctx, "muon_HT1000"));
    jet_HT1000.reset(new JetHists(ctx, "jet_HT1000"));
    event_HT1000.reset(new EventHists(ctx, "event_HT1000"));
    tau_HT1000.reset(new TauHists(ctx, "tau_HT1000"));
    lq_HT1000.reset(new LQAnalysisHists(ctx, "lq_HT1000"));

    ele_HT1300.reset(new ElectronHists(ctx, "ele_HT1300"));
    muon_HT1300.reset(new MuonHists(ctx, "muon_HT1300"));
    jet_HT1300.reset(new JetHists(ctx, "jet_HT1300"));
    event_HT1300.reset(new EventHists(ctx, "event_HT1300"));
    tau_HT1300.reset(new TauHists(ctx, "tau_HT1300"));
    lq_HT1300.reset(new LQAnalysisHists(ctx, "lq_HT1300"));

    ele_HT1800.reset(new ElectronHists(ctx, "ele_HT1800"));
    muon_HT1800.reset(new MuonHists(ctx, "muon_HT1800"));
    jet_HT1800.reset(new JetHists(ctx, "jet_HT1800"));
    event_HT1800.reset(new EventHists(ctx, "event_HT1800"));
    tau_HT1800.reset(new TauHists(ctx, "tau_HT1800"));
    lq_HT1800.reset(new LQAnalysisHists(ctx, "lq_HT1800"));


}


bool LQAnalysisMuModule::process(Event & event) {
    // This is the main procedure, called for each event. Typically,
    // do some pre-processing by calling the modules' process method
    // of the modules constructed in the constructor (1).
    // Then, test whether the event passes some selection and -- if yes --
    // use it to fill the histograms (2).
    // Finally, decide whether or not to keep the event in the output (3);
    // this is controlled by the return value of this method: If it
    // returns true, the event is kept; if it returns false, the event
    // is thrown away.
    
    //cout << "LQAnalysisMuModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
    
    // 1. run all modules; here: only jet cleaning.
    jetcleaner->process(event);
    muoncleaner->process(event);
    electroncleaner->process(event);
    taucleaner->process(event);

    for (auto & mod : pre_modules) {
      mod->process(event);
    }
    
    // 2. test selections and fill histograms
    

    electron_PreSelection->fill(event);
    muon_PreSelection->fill(event);
    tau_PreSelection->fill(event);
    jet_PreSelection->fill(event);
    event_PreSelection->fill(event);
    lq_PreSelection->fill(event);

    auto met = event.met->pt();
    if (met<50) return false;

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

    if (ht<1000) return false;


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
    
    if (!fivejet_sel->passes(event)) return false;
    if (!BTagM->passes(event)) return false;

    ele_HT1000->fill(event);
    muon_HT1000->fill(event);
    jet_HT1000->fill(event);
    event_HT1000->fill(event);
    tau_HT1000->fill(event);
    lq_HT1000->fill(event);


    if (st<1300) return false;
    ele_HT1300->fill(event);
    muon_HT1300->fill(event);
    jet_HT1300->fill(event);
    event_HT1300->fill(event);
    tau_HT1300->fill(event);
    lq_HT1300->fill(event);


    if (st<1800) return false;
    ele_HT1800->fill(event);
    muon_HT1800->fill(event);
    jet_HT1800->fill(event);
    event_HT1800->fill(event);
    tau_HT1800->fill(event);
    lq_HT1800->fill(event);




    
    // 3. decide whether or not to keep the current event in the output:
    return complete_selection;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisMuModule)