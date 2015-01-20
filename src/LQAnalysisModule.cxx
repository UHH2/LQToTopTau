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

using namespace std;
using namespace uhh2;

/** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
 * 
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class LQAnalysisModule: public AnalysisModule {
public:
    
    explicit LQAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
    
    
  std::unique_ptr<JetCleaner> jetcleaner;
  std::unique_ptr<MuonIDKinematic> muonidkinematic;
  std::unique_ptr<MuonCleaner> muoncleaner;
  std::unique_ptr<TauCleaner> taucleaner;
  std::unique_ptr<ElectronCleaner> electroncleaner;
  
  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  std::unique_ptr<Selection> njet_sel, bsel, ntau_sel;
  std::vector<std::unique_ptr<Selection> > fullhad_sel;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_nocuts, h_njet, h_ntau, h_bsel, h_ele, h_event, h_muon, h_jet, h_tau, h_ele_full, h_tau_full, h_event_full, h_jet_full, h_muon_full;

  JetId BTagMedium;


};


LQAnalysisModule::LQAnalysisModule(Context & ctx){
    // In the constructor, the typical tasks are to create
    // other modules like cleaners (1), selections (2) and Hist classes (3).
    // But you can do more and e.g. access the configuration, as shown below.
    
    cout << "Hello World from LQAnalysisModule!" << endl;
    
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
    muoncleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), MuonIDKinematic(20.0, 2.4))));
    electroncleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_CSA14_50ns_medium, PtEtaCut(20.0, 2.5))));
    taucleaner.reset(new TauCleaner(AndId<Tau>(TauIDMedium(), PtEtaCut(20.0, 2.1))));
    BTagMedium = CSVBTag(CSVBTag::WP_MEDIUM);


    // 2. set up selections:
    njet_sel.reset(new NJetSelection(2,-1));
    ntau_sel.reset(new NTauSelection(1,-1));
    //bsel.reset(new NBTagSelection(1));

    int n_cuts = 3;
    fullhad_sel.resize(n_cuts);
    fullhad_sel[0].reset(new NJetSelection(5,-1));
    fullhad_sel[1].reset(new NTauSelection(1));
    fullhad_sel[2].reset(new NJetSelection(1,999,BTagMedium));
    //fullhad_sel[2].reset(new METCut(100,-1));
    

    // 3. Set up Hists classes:
    h_nocuts.reset(new LQAnalysisHists(ctx, "NoCuts"));
    h_njet.reset(new LQAnalysisHists(ctx, "Njet"));
    h_ntau.reset(new TauHists(ctx, "Ntau"));
    h_bsel.reset(new TauHists(ctx, "Bsel"));
    h_ele.reset(new ElectronHists(ctx, "ele_nocuts"));
    h_muon.reset(new MuonHists(ctx, "muon_nocuts"));
    h_tau.reset(new TauHists(ctx, "tau_nocuts"));
    h_jet.reset(new JetHists(ctx, "jet_nocuts"));
    h_event.reset(new EventHists(ctx, "event_nocuts"));
    h_ele_full.reset(new ElectronHists(ctx, "ele_fullselection"));
    h_muon_full.reset(new MuonHists(ctx, "muon_fullselection"));
    h_jet_full.reset(new JetHists(ctx, "jet_fullselection"));
    h_event_full.reset(new EventHists(ctx, "event_fullselection"));
    h_tau_full.reset(new TauHists(ctx, "tau_fullselection"));


}


bool LQAnalysisModule::process(Event & event) {
    // This is the main procedure, called for each event. Typically,
    // do some pre-processing by calling the modules' process method
    // of the modules constructed in the constructor (1).
    // Then, test whether the event passes some selection and -- if yes --
    // use it to fill the histograms (2).
    // Finally, decide whether or not to keep the event in the output (3);
    // this is controlled by the return value of this method: If it
    // returns true, the event is kept; if it returns false, the event
    // is thrown away.
    
    //cout << "LQAnalysisModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
    
    // 1. run all modules; here: only jet cleaning.
    jetcleaner->process(event);
<<<<<<< HEAD
    muoncleaner->process(event);
    electroncleaner->process(event);
    taucleaner->process(event);

=======
    
    //auto met = event.met->pt();
>>>>>>> 82d9b8e4d08df0845834203bafdd380fc5743e12
    // 2. test selections and fill histograms
    
    h_nocuts->fill(event);
    
    bool njet_selection = njet_sel->passes(event);
    if(njet_selection){
        h_njet->fill(event);
    }

    bool ntau_selection = ntau_sel->passes(event);
    if(ntau_selection){
        h_ntau->fill(event);
    }


    /*
    bool bjet_selection = bsel->passes(event);
    if(bjet_selection){
      h_bsel->fill(event);
    }
    */

    h_ele->fill(event);
    h_event->fill(event);
    h_muon->fill(event);
    h_tau->fill(event);
    h_jet->fill(event);


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
    }
    





    
    // 3. decide whether or not to keep the current event in the output:
    return njet_selection && complete_selection;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisModule)
