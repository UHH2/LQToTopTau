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



#include "TH1F.h"


using namespace std;
using namespace uhh2;

/** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
 * 
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class LQAnalysisPreModule: public AnalysisModule {
public:
    
    explicit LQAnalysisPreModule(Context & ctx);
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
  std::unique_ptr<Selection> njet_sel, bsel, ntau_sel, nmuon_sel, nele_sel;
  std::vector<std::unique_ptr<Selection> > fullhad_sel;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  //std::unique_ptr<Hists> h_njet, h_ntau, h_bsel, h_ele, h_event, h_muon, h_jet, h_tau, h_ele_full, h_tau_full, h_event_full, h_jet_full, h_muon_full, h_test1, h_test2;

  std::unique_ptr<Hists> h_lq_nocut, h_tau_nocut, h_mu_nocut, h_ele_nocut, h_jet_nocut, h_event_nocut;
  std::unique_ptr<Hists> h_lq_cleaner, h_tau_cleaner, h_mu_cleaner, h_ele_cleaner, h_jet_cleaner, h_event_cleaner;
  std::unique_ptr<Hists> h_lq_MET50Only, h_tau_MET50Only, h_mu_MET50Only, h_ele_MET50Only, h_jet_MET50Only, h_event_MET50Only;
  std::unique_ptr<Hists> h_lq_TauOnly, h_tau_TauOnly, h_mu_TauOnly, h_ele_TauOnly, h_jet_TauOnly, h_event_TauOnly;
  std::unique_ptr<Hists> h_lq_TwoJetsOnly, h_tau_TwoJetsOnly, h_mu_TwoJetsOnly, h_ele_TwoJetsOnly, h_jet_TwoJetsOnly, h_event_TwoJetsOnly;
  std::unique_ptr<Hists> h_lq_ST350Only, h_tau_ST350Only, h_mu_ST350Only, h_ele_ST350Only, h_jet_ST350Only, h_event_ST350Only;
  std::unique_ptr<Hists> h_lq_PreSel, h_tau_PreSel, h_mu_PreSel, h_ele_PreSel, h_jet_PreSel, h_event_PreSel;


  //std::unique_ptr<Hists> h_before, h_after;

  JetId BTagMedium;
  MuonId muid;
  ElectronId eleid;


};


LQAnalysisPreModule::LQAnalysisPreModule(Context & ctx){
    // In the constructor, the typical tasks are to create
    // other modules like cleaners (1), selections (2) and Hist classes (3).
    // But you can do more and e.g. access the configuration, as shown below.
    
    cout << "Hello World from LQAnalysisPreModule!" << endl;
    
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
    muoncleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), MuonIDKinematic(30.0, 2.4))));
    electroncleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_PHYS14_25ns_medium, PtEtaCut(30.0, 2.5))));
    taucleaner.reset(new TauCleaner(AndId<Tau>(TauIDMedium(), PtEtaCut(20.0, 2.1))));
    //taucleaner.reset(new TauCleaner(AndId<Tau>(TauIDMedium(), PtEtaCut(0.0, 6.0))));
    BTagMedium = CSVBTag(CSVBTag::WP_MEDIUM);
    jetleptoncleaner.reset(new JetLeptonCleaner(JERFiles::PHYS14_L123_MC));
    muid = MuonIDTight();
    eleid = ElectronID_PHYS14_25ns_medium;

    // 2. set up selections:
    njet_sel.reset(new NJetSelection(2,-1));
    ntau_sel.reset(new NTauSelection(1,-1));
    nmuon_sel.reset(new NMuonSelection(0,0));
    nele_sel.reset(new NElectronSelection(0,0));
    bsel.reset(new NBTagSelection(1,-1));

    int n_cuts = 3;
    fullhad_sel.resize(n_cuts);
    fullhad_sel[0].reset(new NJetSelection(2,-1));
    fullhad_sel[1].reset(new NTauSelection(1));
    //fullhad_sel[2].reset(new NJetSelection(1,999,BTagMedium));
    fullhad_sel[2].reset(new NJetCut(2,-1,50,3.0));
    //fullhad_sel[4].reset(new METCut(100,-1));

    // 3. Set up Hists classes:
    //h_nocuts.reset(new LQAnalysisHists(ctx, "NoCuts"));
    /*h_before.reset(new LQAnalysisHists(ctx, "bfore"));
      h_after.reset(new LQAnalysisHists(ctx, "after"));
    h_test1.reset(new JetHists(ctx, "jet_test1"));
    h_test2.reset(new JetHists(ctx, "jet_test2"));*/

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

    h_lq_MET50Only.reset(new LQAnalysisHists(ctx, "LQPreMod_LQ_MET50Only"));
    h_tau_MET50Only.reset(new TauHists(ctx, "LQPreMod_Taus_MET50Only"));
    h_mu_MET50Only.reset(new MuonHists(ctx, "LQPreMod_Muons_MET50Only"));
    h_ele_MET50Only.reset(new ElectronHists(ctx, "LQPreMod_Electrons_MET50Only"));
    h_jet_MET50Only.reset(new JetHists(ctx, "LQPreMod_Jets_MET50Only"));
    h_event_MET50Only.reset(new EventHists(ctx, "LQPreMod_Events_MET50Only"));

    h_lq_TauOnly.reset(new LQAnalysisHists(ctx, "LQPreMod_LQ_TauOnly"));
    h_tau_TauOnly.reset(new TauHists(ctx, "LQPreMod_Taus_TauOnly"));
    h_mu_TauOnly.reset(new MuonHists(ctx, "LQPreMod_Muons_TauOnly"));
    h_ele_TauOnly.reset(new ElectronHists(ctx, "LQPreMod_Electrons_TauOnly"));
    h_jet_TauOnly.reset(new JetHists(ctx, "LQPreMod_Jets_TauOnly"));
    h_event_TauOnly.reset(new EventHists(ctx, "LQPreMod_Events_TauOnly"));

    h_lq_TwoJetsOnly.reset(new LQAnalysisHists(ctx, "LQPreMod_LQ_TwoJetsOnly"));
    h_tau_TwoJetsOnly.reset(new TauHists(ctx, "LQPreMod_Taus_TwoJetsOnly"));
    h_mu_TwoJetsOnly.reset(new MuonHists(ctx, "LQPreMod_Muons_TwoJetsOnly"));
    h_ele_TwoJetsOnly.reset(new ElectronHists(ctx, "LQPreMod_Electrons_TwoJetsOnly"));
    h_jet_TwoJetsOnly.reset(new JetHists(ctx, "LQPreMod_Jets_TwoJetsOnly"));
    h_event_TwoJetsOnly.reset(new EventHists(ctx, "LQPreMod_Events_TwoJetsOnly"));

    h_lq_ST350Only.reset(new LQAnalysisHists(ctx, "LQPreMod_LQ_ST350Only"));
    h_tau_ST350Only.reset(new TauHists(ctx, "LQPreMod_Taus_ST350Only"));
    h_mu_ST350Only.reset(new MuonHists(ctx, "LQPreMod_Muons_ST350Only"));
    h_ele_ST350Only.reset(new ElectronHists(ctx, "LQPreMod_Electrons_ST3500Only"));
    h_jet_ST350Only.reset(new JetHists(ctx, "LQPreMod_Jets_ST350Only"));
    h_event_ST350Only.reset(new EventHists(ctx, "LQPreMod_Events_ST350Only"));

    h_lq_PreSel.reset(new LQAnalysisHists(ctx, "LQPreMod_LQ_PreSel"));
    h_tau_PreSel.reset(new TauHists(ctx, "LQPreMod_Taus_PreSel"));
    h_mu_PreSel.reset(new MuonHists(ctx, "LQPreMod_Muons_PreSel"));
    h_ele_PreSel.reset(new ElectronHists(ctx, "LQPreMod_Electrons_PreSel"));
    h_jet_PreSel.reset(new JetHists(ctx, "LQPreMod_Jets_PreSel"));
    h_event_PreSel.reset(new EventHists(ctx, "LQPreMod_Events_PreSel"));





}


bool LQAnalysisPreModule::process(Event & event) {
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

  //h_test1->fill(event);    

  // fill hists before all selection cuts
    h_lq_nocut->fill(event);
    h_tau_nocut->fill(event);
    h_mu_nocut->fill(event);
    h_ele_nocut->fill(event);
    h_jet_nocut->fill(event);
    h_event_nocut->fill(event);


   // 1. run all modules; here: only jet cleaning.



    //if(!nele_sel->passes(event)) return false;
    muoncleaner->process(event);
    electroncleaner->process(event);
    jetleptoncleaner->process(event);

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

    jetcleaner->process(event);


    //h_test2->fill(event);
    /*    
    for(const auto & tau : *event.taus){
      if(!tau.decayModeFinding() || !tau.byLooseCombinedIsolationDeltaBetaCorr3Hits()) return false;
    }

    h_before->fill(event);
    taucleaner->process(event);
    h_after->fill(event);
    */


    // fill hists after cleaning modules 
    h_lq_cleaner->fill(event);
    h_tau_cleaner->fill(event);
    h_mu_cleaner->fill(event);
    h_ele_cleaner->fill(event);
    h_jet_cleaner->fill(event);
    h_event_cleaner->fill(event);

    // 2. test selections and fill histograms
    
    //h_nocuts->fill(event);


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


    if(!nmuon_sel->passes(event)) return false;
    if(!nele_sel->passes(event)) return false;
    
 

    if(!njet_sel->passes(event)) return false;
    const auto jets = event.jets;
    if(jets->size() > 0){
      const auto & jet = (*jets)[0];
      if(jet.pt()<50) return false;
    }
    if(jets->size() > 1){
      const auto & jet = (*jets)[1];
      if(jet.pt()<50) return false;
    }

    h_lq_TwoJetsOnly->fill(event);
    h_tau_TwoJetsOnly->fill(event);
    h_mu_TwoJetsOnly->fill(event);
    h_ele_TwoJetsOnly->fill(event);
    h_jet_TwoJetsOnly->fill(event);
    h_event_TwoJetsOnly->fill(event);
    
    if(st<500) return false;

    h_lq_ST350Only->fill(event);
    h_tau_ST350Only->fill(event);
    h_mu_ST350Only->fill(event);
    h_ele_ST350Only->fill(event);
    h_jet_ST350Only->fill(event);
    h_event_ST350Only->fill(event);
    
    if(met<50) return false;

    h_lq_MET50Only->fill(event);
    h_tau_MET50Only->fill(event);
    h_mu_MET50Only->fill(event);
    h_ele_MET50Only->fill(event);
    h_jet_MET50Only->fill(event);
    h_event_MET50Only->fill(event);
    
    if(!ntau_sel->passes(event)) return false;
 
    h_lq_TauOnly->fill(event);
    h_tau_TauOnly->fill(event);
    h_mu_TauOnly->fill(event);
    h_ele_TauOnly->fill(event);
    h_jet_TauOnly->fill(event);
    h_event_TauOnly->fill(event);
    
   


    bool complete_selection = true;
    std::vector<bool> v_accept(fullhad_sel.size());
    for (unsigned i=0; i<fullhad_sel.size(); ++i) {
      bool accept = fullhad_sel[i]->passes(event);
      v_accept[i] = accept;
      if (!accept) {
	complete_selection = false;
      }
    }
    /*if (complete_selection) {
      h_ele_full->fill(event);
      h_muon_full->fill(event);
      h_jet_full->fill(event);
      h_tau_full->fill(event);
      h_event_full->fill(event);
      }*/
    

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
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisPreModule)
