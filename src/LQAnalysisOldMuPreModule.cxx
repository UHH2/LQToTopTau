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
#include "UHH2/common/include/JetCorrections.h"

using namespace std;
using namespace uhh2;

/** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
 * 
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class LQAnalysisOldMuPreModule: public AnalysisModule {
public:
    
    explicit LQAnalysisOldMuPreModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
    
    
  std::unique_ptr<JetCleaner> jetcleaner;
  std::unique_ptr<MuonIDKinematic> muonidkinematic;
  std::unique_ptr<MuonCleaner> muoncleaner;
  std::unique_ptr<TauCleaner> taucleaner;
  std::unique_ptr<ElectronCleaner> electroncleaner;
  std::unique_ptr<JetLeptonCleaner> jetleptoncleaner;
  
  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  std::unique_ptr<Selection> njet_sel, bsel, ntau_sel, nmuon_sel;
  std::vector<std::unique_ptr<Selection> > fullhad_sel;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_lq_nocut, h_tau_nocut, h_mu_nocut, h_ele_nocut, h_jet_nocut, h_event_nocut;
  std::unique_ptr<Hists> h_lq_cleaner, h_tau_cleaner, h_mu_cleaner, h_ele_cleaner, h_jet_cleaner, h_event_cleaner;
  std::unique_ptr<Hists> h_lq_beforeST, h_tau_beforeST, h_mu_beforeST, h_ele_beforeST, h_jet_beforeST, h_event_beforeST;
  std::unique_ptr<Hists> h_lq_PreSel, h_tau_PreSel, h_mu_PreSel, h_ele_PreSel, h_jet_PreSel, h_event_PreSel;

  JetId BTagMedium;


};


LQAnalysisOldMuPreModule::LQAnalysisOldMuPreModule(Context & ctx){
    // In the constructor, the typical tasks are to create
    // other modules like cleaners (1), selections (2) and Hist classes (3).
    // But you can do more and e.g. access the configuration, as shown below.
    
    cout << "Hello World from LQAnalysisOldMuPreModule!" << endl;
    
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
    jetcleaner.reset(new JetCleaner(30.0, 2.5));
    muonidkinematic.reset(new MuonIDKinematic(30.0,3.0));
    muoncleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), MuonIDKinematic(30.0, 2.1))));
    electroncleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_PHYS14_25ns_medium, PtEtaCut(30.0, 2.5))));
    taucleaner.reset(new TauCleaner(AndId<Tau>(TauIDDecayModeFinding(), PtEtaCut(20.0, 2.1))));
    BTagMedium = CSVBTag(CSVBTag::WP_MEDIUM);
    //jetleptoncleaner.reset(new JetLeptonCleaner(AndId<Muon>(MuonIDTight(), MuonIDKinematic(20.0, 2.4)),AndId<Electron>(ElectronID_CSA14_50ns_medium, PtEtaCut(30.0, 2.5))));
    jetleptoncleaner.reset(new JetLeptonCleaner(JERFiles::PHYS14_L123_MC));

    // 2. set up selections:
    njet_sel.reset(new NJetCut(2,-1,50,2.5));
    ntau_sel.reset(new NTauSelection(1,-1));
    nmuon_sel.reset(new NMuonSelection(1,-1));
    bsel.reset(new NBTagSelection(1,-1));

    int n_cuts = 4;
    fullhad_sel.resize(n_cuts);
    fullhad_sel[0].reset(new NJetSelection(2,-1));
    fullhad_sel[1].reset(new NTauSelection(1));
    //fullhad_sel[2].reset(new NJetSelection(1,999,BTagMedium));
    fullhad_sel[2].reset(new NJetCut(2,-1,50,3.0));
    fullhad_sel[3].reset(new NMuonSelection(1));
    //fullhad_sel[4].reset(new METCut(100,-1));



    // 3. Set up Hists classes:
    h_lq_nocut.reset(new LQAnalysisHists(ctx, "LQPreMod_LQ_NoCuts"));
    h_tau_nocut.reset(new TauHists(ctx, "LQPreMod_Taus_NoCuts"));
    h_mu_nocut.reset(new MuonHists(ctx, "LQPreMod_Muons_NoCuts"));
    h_ele_nocut.reset(new ElectronHists(ctx, "LQPreMod_Electrons_NoCuts"));
    h_jet_nocut.reset(new JetHists(ctx, "LQPreMod_Jets_NoCuts"));
    h_event_nocut.reset(new EventHists(ctx, "LQPreMod_Events_NoCuts"));

    h_lq_cleaner.reset(new LQAnalysisHists(ctx, "LQPreMod_LQ_Cleaner"));
    h_tau_cleaner.reset(new TauHists(ctx, "LQPreMod_Taus_Cleaner"));
    h_mu_cleaner.reset(new MuonHists(ctx, "LQPreMod_Muons_Cleaner"));
    h_ele_cleaner.reset(new ElectronHists(ctx, "LQPreMod_Electrons_Cleaner"));
    h_jet_cleaner.reset(new JetHists(ctx, "LQPreMod_Jets_Cleaner"));
    h_event_cleaner.reset(new EventHists(ctx, "LQPreMod_Events_Cleaner"));

    h_lq_beforeST.reset(new LQAnalysisHists(ctx, "LQPreMod_LQ_beforeST"));
    h_tau_beforeST.reset(new TauHists(ctx, "LQPreMod_Taus_beforeST"));
    h_mu_beforeST.reset(new MuonHists(ctx, "LQPreMod_Muons_beforeST"));
    h_ele_beforeST.reset(new ElectronHists(ctx, "LQPreMod_Electrons_beforeST"));
    h_jet_beforeST.reset(new JetHists(ctx, "LQPreMod_Jets_beforeST"));
    h_event_beforeST.reset(new EventHists(ctx, "LQPreMod_Events_beforeST"));

    h_lq_PreSel.reset(new LQAnalysisHists(ctx, "LQPreMod_LQ_PreSel"));
    h_tau_PreSel.reset(new TauHists(ctx, "LQPreMod_Taus_PreSel"));
    h_mu_PreSel.reset(new MuonHists(ctx, "LQPreMod_Muons_PreSel"));
    h_ele_PreSel.reset(new ElectronHists(ctx, "LQPreMod_Electrons_PreSel"));
    h_jet_PreSel.reset(new JetHists(ctx, "LQPreMod_Jets_PreSel"));
    h_event_PreSel.reset(new EventHists(ctx, "LQPreMod_Events_PreSel"));



   

}


bool LQAnalysisOldMuPreModule::process(Event & event) {
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





  // fill hists before all selection cuts
    h_lq_nocut->fill(event);
    h_tau_nocut->fill(event);
    h_mu_nocut->fill(event);
    h_ele_nocut->fill(event);
    h_jet_nocut->fill(event);
    h_event_nocut->fill(event);

    //cleaning modules
    taucleaner->process(event);
    muoncleaner->process(event);
    electroncleaner->process(event);
    jetleptoncleaner->process(event);

    for(unsigned int i=0; i<event.jets->size(); i++){
      Jet jet = event.jets->at(i);
      for(const auto & tau : *event.taus){
	if(deltaR(tau,jet)<0.4){
	  event.jets->erase(event.jets->begin()+i);
	}
      }
    }

    jetcleaner->process(event);

    // fill hists after cleaning modules 
    h_lq_cleaner->fill(event);
    h_tau_cleaner->fill(event);
    h_mu_cleaner->fill(event);
    h_ele_cleaner->fill(event);
    h_jet_cleaner->fill(event);
    h_event_cleaner->fill(event);

    // define met
    auto met = event.met->pt();

    // define ht and st
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



    // 2. test selections and fill histograms
    bool complete_selection = true;
    std::vector<bool> v_accept(fullhad_sel.size());
    for (unsigned i=0; i<fullhad_sel.size(); ++i) {
      bool accept = fullhad_sel[i]->passes(event);
      v_accept[i] = accept;
      if (!accept) {
	complete_selection = false;
      }
    }


    if(!njet_sel->passes(event)) return false;
    if(!ntau_sel->passes(event)) return false;
    if(!nmuon_sel->passes(event)) return false;

    h_lq_beforeST->fill(event);
    h_tau_beforeST->fill(event);
    h_mu_beforeST->fill(event);
    h_ele_beforeST->fill(event);
    h_jet_beforeST->fill(event);
    h_event_beforeST->fill(event);

    if(st<350) return false;

    h_lq_PreSel->fill(event);
    h_tau_PreSel->fill(event);
    h_mu_PreSel->fill(event);
    h_ele_PreSel->fill(event);
    h_jet_PreSel->fill(event);
    h_event_PreSel->fill(event);

    
    // 3. decide whether or not to keep the current event in the output:
    return complete_selection;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisOldMuPreModule)
