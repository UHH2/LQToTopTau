#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/common/include/CommonModules.h"
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
#include "UHH2/LQAnalysis/include/LQAnalysisPreHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/PrintingModules.h"
//#include "UHH2/common/test/TestJetLeptonCleaner.cpp"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/TauUncerts.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/EventVariables.h"


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
  std::string region_;

  std::unique_ptr<JetCleaner> jetcleaner;
  std::unique_ptr<MuonIDKinematic> muonidkinematic;
  std::unique_ptr<TauCleaner> taucleaner;
  std::unique_ptr<JetLeptonCleaner> jetleptoncleaner;

  
  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  std::unique_ptr<Selection> njet_sel, ntau_sel, nmuon_sel, nele_sel, lumi_sel, trigger_sel, JetTauCleaner, metcut, stcut;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_lq_nocut, h_tau_nocut, h_mu_nocut, h_ele_nocut, h_jet_nocut, h_event_nocut, h_lumi_nocut;
  std::unique_ptr<Hists> h_lq_trigger, h_tau_trigger, h_mu_trigger, h_ele_trigger, h_jet_trigger, h_event_trigger, h_lumi_trigger;
  std::unique_ptr<Hists> h_lq_cleaner, h_tau_cleaner, h_mu_cleaner, h_ele_cleaner, h_jet_cleaner, h_event_cleaner;
  std::unique_ptr<Hists> h_lq_MET50Only, h_tau_MET50Only, h_mu_MET50Only, h_ele_MET50Only, h_jet_MET50Only, h_event_MET50Only;
  std::unique_ptr<Hists> h_lq_MuonOnly, h_tau_MuonOnly, h_mu_MuonOnly, h_ele_MuonOnly, h_jet_MuonOnly, h_event_MuonOnly;
  std::unique_ptr<Hists> h_lq_TauOnly, h_tau_TauOnly, h_mu_TauOnly, h_ele_TauOnly, h_jet_TauOnly, h_event_TauOnly;
  std::unique_ptr<Hists> h_lq_TwoJetsOnly, h_tau_TwoJetsOnly, h_mu_TwoJetsOnly, h_ele_TwoJetsOnly, h_jet_TwoJetsOnly, h_event_TwoJetsOnly;
  std::unique_ptr<Hists> h_lq_ST350Only, h_tau_ST350Only, h_mu_ST350Only, h_ele_ST350Only, h_jet_ST350Only, h_event_ST350Only;
  std::unique_ptr<Hists> h_lq_PreSel, h_tau_PreSel, h_mu_PreSel, h_ele_PreSel, h_jet_PreSel, h_event_PreSel, h_lumi_PreSel;


  std::unique_ptr<CommonModules> common;

  JetId jet_id;
  PtEtaCut* pteta;
  MuonId MuIso, MuId;
  ElectronId EleIso, EleId;
  TauId TauonId;

  bool is_data;

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

  auto dataset_type = ctx.get("dataset_type");
  is_data = dataset_type == "DATA";

  // 1. setup other modules.
  EleId = AndId<Electron>(ElectronID_Spring15_25ns_medium, PtEtaCut(30.0, 2.5));
  MuId = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4),MuonIso(0.15));

  //switch tauId depending on whether we're in the signal or control region
  region_ = ctx.get("Region");
  if(region_!="SR" && region_!="CR")
    throw std::runtime_error("undefined argument for 'region' (must be 'SR' or 'CR'): "+region_);
  if(region_ == "SR"){
    TauonId = AndId<Tau>(TauIDMedium(), PtEtaCut(20.0, 2.1));
  }
  if(region_ == "CR"){
    TauonId = AndId<Tau>(TauIDDecayModeFinding(), PtEtaCut(20.0, 2.1));
  }

  pteta = new PtEtaCut(30, 2.5);
  jet_id = *pteta;
  jetcleaner.reset(new JetCleaner(ctx, jet_id));
  common.reset(new CommonModules());
  //common->disable_mcpileupreweight();
  //common->disable_metfilters();
  //common->disable_pvfilter();
  common->disable_lumisel();
  common->disable_jersmear();
  common->switch_jetlepcleaner(true);
  common->set_electron_id(EleId);
  common->set_muon_id(MuId);
  common->set_tau_id(TauonId);
  common->init(ctx);


  // 2. set up selections:
  njet_sel.reset(new NJetSelection(2,-1));
  ntau_sel.reset(new NTauSelection(2,-1));
  nmuon_sel.reset(new NMuonSelection(0,0));
  nele_sel.reset(new NElectronSelection(0,0));
  metcut.reset(new METCut(50,-1));
  stcut.reset(new HtSelection(300,-1));
  trigger_sel.reset(new TriggerSelection("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v*"));
  //trigger_sel.reset(new TriggerSelection("HLT_AK4PFJet30_v*"));
  lumi_sel.reset(new LumiSelection(ctx));
  JetTauCleaner.reset(new JetTauCleaning());

  // 3. Set up Hists classes:
  h_lq_nocut.reset(new LQAnalysisPreHists(ctx, "LQPreMod_LQ_NoCuts"));
  h_tau_nocut.reset(new TauHists(ctx, "LQPreMod_Taus_NoCuts"));
  h_mu_nocut.reset(new MuonHists(ctx, "LQPreMod_Muons_NoCuts"));
  h_ele_nocut.reset(new ElectronHists(ctx, "LQPreMod_Electrons_NoCuts"));
  h_jet_nocut.reset(new JetHists(ctx, "LQPreMod_Jets_NoCuts"));
  h_event_nocut.reset(new EventHists(ctx, "LQPreMod_Events_NoCuts"));
  h_lumi_nocut.reset(new LuminosityHists(ctx, "LQPreMod_Lumi_NoCuts"));

  h_lq_trigger.reset(new LQAnalysisPreHists(ctx, "LQPreMod_LQ_Trigger"));
  h_tau_trigger.reset(new TauHists(ctx, "LQPreMod_Taus_Trigger"));
  h_mu_trigger.reset(new MuonHists(ctx, "LQPreMod_Muons_Trigger"));
  h_ele_trigger.reset(new ElectronHists(ctx, "LQPreMod_Electrons_Trigger"));
  h_jet_trigger.reset(new JetHists(ctx, "LQPreMod_Jets_Trigger"));
  h_event_trigger.reset(new EventHists(ctx, "LQPreMod_Events_Trigger"));
  h_lumi_trigger.reset(new LuminosityHists(ctx, "LQPreMod_Lumi_Trigger"));

  h_lq_cleaner.reset(new LQAnalysisPreHists(ctx, "LQPreMod_LQ_Cleaner"));
  h_tau_cleaner.reset(new TauHists(ctx, "LQPreMod_Taus_Cleaner"));
  h_mu_cleaner.reset(new MuonHists(ctx, "LQPreMod_Muons_Cleaner"));
  h_ele_cleaner.reset(new ElectronHists(ctx, "LQPreMod_Electrons_Cleaner"));
  h_jet_cleaner.reset(new JetHists(ctx, "LQPreMod_Jets_Cleaner"));
  h_event_cleaner.reset(new EventHists(ctx, "LQPreMod_Events_Cleaner"));

  h_lq_MET50Only.reset(new LQAnalysisPreHists(ctx, "LQPreMod_LQ_MET50Only"));
  h_tau_MET50Only.reset(new TauHists(ctx, "LQPreMod_Taus_MET50Only"));
  h_mu_MET50Only.reset(new MuonHists(ctx, "LQPreMod_Muons_MET50Only"));
  h_ele_MET50Only.reset(new ElectronHists(ctx, "LQPreMod_Electrons_MET50Only"));
  h_jet_MET50Only.reset(new JetHists(ctx, "LQPreMod_Jets_MET50Only"));
  h_event_MET50Only.reset(new EventHists(ctx, "LQPreMod_Events_MET50Only"));

  h_lq_MuonOnly.reset(new LQAnalysisPreHists(ctx, "LQPreMod_LQ_MuonOnly"));
  h_tau_MuonOnly.reset(new TauHists(ctx, "LQPreMod_Taus_MuonOnly"));
  h_mu_MuonOnly.reset(new MuonHists(ctx, "LQPreMod_Muons_MuonOnly"));
  h_ele_MuonOnly.reset(new ElectronHists(ctx, "LQPreMod_Electrons_MuonOnly"));
  h_jet_MuonOnly.reset(new JetHists(ctx, "LQPreMod_Jets_MuonOnly"));
  h_event_MuonOnly.reset(new EventHists(ctx, "LQPreMod_Events_MuonOnly"));

  h_lq_TauOnly.reset(new LQAnalysisPreHists(ctx, "LQPreMod_LQ_TauOnly"));
  h_tau_TauOnly.reset(new TauHists(ctx, "LQPreMod_Taus_TauOnly"));
  h_mu_TauOnly.reset(new MuonHists(ctx, "LQPreMod_Muons_TauOnly"));
  h_ele_TauOnly.reset(new ElectronHists(ctx, "LQPreMod_Electrons_TauOnly"));
  h_jet_TauOnly.reset(new JetHists(ctx, "LQPreMod_Jets_TauOnly"));
  h_event_TauOnly.reset(new EventHists(ctx, "LQPreMod_Events_TauOnly"));

  h_lq_TwoJetsOnly.reset(new LQAnalysisPreHists(ctx, "LQPreMod_LQ_TwoJetsOnly"));
  h_tau_TwoJetsOnly.reset(new TauHists(ctx, "LQPreMod_Taus_TwoJetsOnly"));
  h_mu_TwoJetsOnly.reset(new MuonHists(ctx, "LQPreMod_Muons_TwoJetsOnly"));
  h_ele_TwoJetsOnly.reset(new ElectronHists(ctx, "LQPreMod_Electrons_TwoJetsOnly"));
  h_jet_TwoJetsOnly.reset(new JetHists(ctx, "LQPreMod_Jets_TwoJetsOnly"));
  h_event_TwoJetsOnly.reset(new EventHists(ctx, "LQPreMod_Events_TwoJetsOnly"));

  h_lq_ST350Only.reset(new LQAnalysisPreHists(ctx, "LQPreMod_LQ_ST350Only"));
  h_tau_ST350Only.reset(new TauHists(ctx, "LQPreMod_Taus_ST350Only"));
  h_mu_ST350Only.reset(new MuonHists(ctx, "LQPreMod_Muons_ST350Only"));
  h_ele_ST350Only.reset(new ElectronHists(ctx, "LQPreMod_Electrons_ST3500Only"));
  h_jet_ST350Only.reset(new JetHists(ctx, "LQPreMod_Jets_ST350Only"));
  h_event_ST350Only.reset(new EventHists(ctx, "LQPreMod_Events_ST350Only"));

  h_lq_PreSel.reset(new LQAnalysisPreHists(ctx, "LQPreMod_LQ_PreSel"));
  h_tau_PreSel.reset(new TauHists(ctx, "LQPreMod_Taus_PreSel"));
  h_mu_PreSel.reset(new MuonHists(ctx, "LQPreMod_Muons_PreSel"));
  h_ele_PreSel.reset(new ElectronHists(ctx, "LQPreMod_Electrons_PreSel"));
  h_jet_PreSel.reset(new JetHists(ctx, "LQPreMod_Jets_PreSel"));
  h_event_PreSel.reset(new EventHists(ctx, "LQPreMod_Events_PreSel"));
  h_lumi_PreSel.reset(new LuminosityHists(ctx, "LQPreMod_Lumi_PreSel"));

   

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
    

  if(is_data){
    if(!lumi_sel->passes(event)) return false;
  }
  
  // fill hists before all selection cuts
  h_lq_nocut->fill(event);
  h_tau_nocut->fill(event);
  h_mu_nocut->fill(event);
  h_ele_nocut->fill(event);
  h_jet_nocut->fill(event);
  h_event_nocut->fill(event);
  h_lumi_nocut->fill(event);
  
  //print all trigger names
  /*
  for (unsigned int i=0; i<event.get_current_triggernames().size();i++)
    cout<< event.get_current_triggernames()[i]<<"\n";
  cout<<"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
  */

  //trigger sel + hists
  //if(!trigger_sel->passes(event)) return false;  
  h_lq_trigger->fill(event);
  h_tau_trigger->fill(event);
  h_mu_trigger->fill(event);
  h_ele_trigger->fill(event);
  h_jet_trigger->fill(event);
  h_event_trigger->fill(event);
  h_lumi_trigger->fill(event);

  //cleaning modules
  bool pass_common = common->process(event);
  if(!pass_common) return false;

  /*
  taucleaner->process(event);
  */

  if(!JetTauCleaner->passes(event)) return false;
  jetcleaner->process(event);
 
  // fill hists after cleaning modules   
  h_lq_cleaner->fill(event);
  h_tau_cleaner->fill(event);
  h_mu_cleaner->fill(event);
  h_ele_cleaner->fill(event);
  h_jet_cleaner->fill(event);
  h_event_cleaner->fill(event);

  if(!ntau_sel->passes(event)) return false;
  h_lq_TauOnly->fill(event);
  h_tau_TauOnly->fill(event);
  h_mu_TauOnly->fill(event);
  h_ele_TauOnly->fill(event);
  h_jet_TauOnly->fill(event);
  h_event_TauOnly->fill(event);

  if(!nmuon_sel->passes(event)) return false;
  if(!nele_sel->passes(event)) return false;
  h_lq_MuonOnly->fill(event);
  h_tau_MuonOnly->fill(event);
  h_mu_MuonOnly->fill(event);
  h_ele_MuonOnly->fill(event);
  h_jet_MuonOnly->fill(event);
  h_event_MuonOnly->fill(event);
    
  if(!njet_sel->passes(event)) return false;
  const auto jets = event.jets;
  if(jets->size() > 0){
    const auto & jet1 = (*jets)[0];
    const auto & jet2 = (*jets)[1];
    if(jet1.pt()<50) return false;
    if(jet2.pt()<50) return false;
  }
  h_lq_TwoJetsOnly->fill(event);
  h_tau_TwoJetsOnly->fill(event);
  h_mu_TwoJetsOnly->fill(event);
  h_ele_TwoJetsOnly->fill(event);
  h_jet_TwoJetsOnly->fill(event);
  h_event_TwoJetsOnly->fill(event);
    
  if(!stcut->passes(event)) return false;
  h_lq_ST350Only->fill(event);
  h_tau_ST350Only->fill(event);
  h_mu_ST350Only->fill(event);
  h_ele_ST350Only->fill(event);
  h_jet_ST350Only->fill(event);
  h_event_ST350Only->fill(event);

  if(!metcut->passes(event)) return false;
  h_lq_MET50Only->fill(event);
  h_tau_MET50Only->fill(event);
  h_mu_MET50Only->fill(event);
  h_ele_MET50Only->fill(event);
  h_jet_MET50Only->fill(event);
  h_event_MET50Only->fill(event);

  h_lq_PreSel->fill(event);
  h_tau_PreSel->fill(event);
  h_mu_PreSel->fill(event);
  h_ele_PreSel->fill(event);
  h_jet_PreSel->fill(event);
  h_event_PreSel->fill(event);
  h_lumi_PreSel->fill(event);
    
  // 3. decide whether or not to keep the current event in the output:
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisPreModule)
