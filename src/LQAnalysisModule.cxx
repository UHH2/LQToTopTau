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
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/LQAnalysis/include/LQAnalysisSelections.h"
#include "UHH2/LQAnalysis/include/LQAnalysisHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/EventVariables.h"

#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQAnalysis/include/TTbarFullhadRecoHypothesisDiscriminators.h"
#include "UHH2/LQAnalysis/include/TTbarFullhadRecoHypothesis.h"
#include "UHH2/LQAnalysis/include/TTbarFullhadReco.h"
#include "UHH2/LQAnalysis/include/LQAllhadHists.h"


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
    
   
  
  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  std::unique_ptr<Selection> njet_sel, twojet_sel, fourjet_sel, TwoBTagL, TwoBTagM, BTagL, BTagM, BTagT, toptag, ntau_sel, mbtau_sel;

  std::vector<std::unique_ptr<AnalysisModule>> pre_modules;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> lq_PreSelection, electron_PreSelection, muon_PreSelection, jet_PreSelection, topjet_PreSelection, tau_PreSelection, event_PreSelection, hypothesis_PreSelection;
  std::unique_ptr<Hists> h_lq_ThreeJets, h_tau_ThreeJets, h_mu_ThreeJets, h_ele_ThreeJets, h_jet_ThreeJets, h_topjet_ThreeJets, h_event_ThreeJets, h_hypothesis_ThreeJets;
  std::unique_ptr<Hists> h_lq_ST400, h_tau_ST400, h_mu_ST400, h_ele_ST400, h_jet_ST400, h_topjet_ST400, h_event_ST400, h_hypothesis_ST400;
  std::unique_ptr<Hists> h_lq_ST500, h_tau_ST500, h_mu_ST500, h_ele_ST500, h_jet_ST500, h_topjet_ST500, h_event_ST500, h_hypothesis_ST500;
  std::unique_ptr<Hists> h_lq_ST600, h_tau_ST600, h_mu_ST600, h_ele_ST600, h_jet_ST600, h_topjet_ST600, h_event_ST600, h_hypothesis_ST600;

  JetId BTagLoose, BTagMedium, BTagTight;
  TopJetId CMSTopTagger;
  MuonId MuIso;

  //std::unique_ptr<AnalysisModule> jetlepcleantest;
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<JetCleaner> jetcleaner;
  std::string Sys_PU;

  bool is_data;
  TauId TauonId;

  std::vector<std::unique_ptr<AnalysisModule>> recomodules;
  std::unique_ptr<AnalysisModule> ttgenprod;
  Event::Handle<TTbarGen> h_ttbargen;
  Event::Handle<std::vector<TTbarFullhadRecoHypothesis>> h_hadr_hyps;

};


LQAnalysisModule::LQAnalysisModule(Context & ctx){   
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
  jetcleaner.reset(new JetCleaner(ctx, 30.0, 2.4));
  BTagLoose = CSVBTag(CSVBTag::WP_LOOSE);
  BTagMedium = CSVBTag(CSVBTag::WP_MEDIUM);
  BTagTight = CSVBTag(CSVBTag::WP_TIGHT);
  CMSTopTagger = CMSTopTag(50,140,250);

  TauonId = AndId<Tau>(TauIDMedium(), PtEtaCut(30.0, 2.1));

  Sys_PU = ctx.get("puvariation");

  common.reset(new CommonModules());
  //common->disable_mcpileupreweight();
  common->disable_lumisel();
  //common->disable_jersmear();
  //common->disable_jec();
  //common->disable_metfilters();
  //common->disable_pvfilter();
  //common->set_electron_id(EleId);
  //common->set_muon_id(MuId);
  //common->set_tau_id(TauonId);
  //common->init(ctx);
  common->init(ctx,Sys_PU);


  // ttbar GEN
  ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
  // ttbar RECO hypotheses
  h_hadr_hyps = ctx.get_handle<std::vector<TTbarFullhadRecoHypothesis>>("HighMassHadronicTTbarFullhadReco");
  recomodules.emplace_back(new HighMassHadronicTTbarReco(ctx));
  recomodules.emplace_back(new TTbarFullhadRecoChi2Discriminator(ctx,"HighMassHadronicTTbarFullhadReco"));

  pre_modules.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));

  // 2. set up selections:
  twojet_sel.reset(new NJetSelection(2,-1));
  njet_sel.reset(new NJetSelection(3,-1));
  fourjet_sel.reset(new NJetSelection(4,-1));
  ntau_sel.reset(new NTauSelection(2,-1));
  TwoBTagL.reset(new NJetSelection(2,999,BTagLoose));
  TwoBTagM.reset(new NJetSelection(2,999,BTagMedium));
  BTagL.reset(new NJetSelection(1,999,BTagLoose));
  BTagM.reset(new NJetSelection(1,999,BTagMedium));
  BTagT.reset(new NJetSelection(1,999,BTagTight));
  mbtau_sel.reset(new MbtauSelection(150,-1));
  toptag.reset(new NTopJetSelection(1,999,CMSTopTagger));


  // 3. Set up Hists classes:
  electron_PreSelection.reset(new ElectronHists(ctx, "LQMod_Electrons_PreSel"));
  muon_PreSelection.reset(new MuonHists(ctx, "LQMod_Muons_PreSel"));
  tau_PreSelection.reset(new TauHists(ctx, "LQMod_Taus_PreSel"));
  jet_PreSelection.reset(new JetHists(ctx, "LQMod_Jets_PreSel"));
  topjet_PreSelection.reset(new TopJetHists(ctx, "LQMod_TopJets_PreSel"));
  event_PreSelection.reset(new EventHists(ctx, "LQMod_Events_PreSel"));
  lq_PreSelection.reset(new LQAnalysisHists(ctx, "LQMod_LQ_PreSel"));
  hypothesis_PreSelection.reset(new LQAllhadHists(ctx, "LQMod_Hypothesis_PreSel"));

  h_ele_ThreeJets.reset(new ElectronHists(ctx, "LQMod_Electrons_ThreeJets"));
  h_mu_ThreeJets.reset(new MuonHists(ctx, "LQMod_Muons_ThreeJets"));
  h_jet_ThreeJets.reset(new JetHists(ctx, "LQMod_Jets_ThreeJets"));
  h_topjet_ThreeJets.reset(new TopJetHists(ctx, "LQMod_TopJets_ThreeJets"));
  h_event_ThreeJets.reset(new EventHists(ctx, "LQMod_Events_ThreeJets"));
  h_tau_ThreeJets.reset(new TauHists(ctx, "LQMod_Taus_ThreeJets"));
  h_lq_ThreeJets.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ThreeJets"));
  h_hypothesis_ThreeJets.reset(new LQAllhadHists(ctx, "LQMod_Hypothesis_ThreeJets"));

  h_ele_ST400.reset(new ElectronHists(ctx, "LQMod_Electrons_ST400"));
  h_mu_ST400.reset(new MuonHists(ctx, "LQMod_Muons_ST400"));
  h_jet_ST400.reset(new JetHists(ctx, "LQMod_Jets_ST400"));
  h_topjet_ST400.reset(new TopJetHists(ctx, "LQMod_TopJets_ST400"));
  h_event_ST400.reset(new EventHists(ctx, "LQMod_Events_ST400"));
  h_tau_ST400.reset(new TauHists(ctx, "LQMod_Taus_ST400"));
  h_lq_ST400.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST400"));
  h_hypothesis_ST400.reset(new LQAllhadHists(ctx, "LQMod_Hypothesis_ST400"));

  h_ele_ST500.reset(new ElectronHists(ctx, "LQMod_Electrons_ST500"));
  h_mu_ST500.reset(new MuonHists(ctx, "LQMod_Muons_ST500"));
  h_jet_ST500.reset(new JetHists(ctx, "LQMod_Jets_ST500"));
  h_topjet_ST500.reset(new TopJetHists(ctx, "LQMod_TopJets_ST500"));
  h_event_ST500.reset(new EventHists(ctx, "LQMod_Events_ST500"));
  h_tau_ST500.reset(new TauHists(ctx, "LQMod_Taus_ST500"));
  h_lq_ST500.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST500"));
  h_hypothesis_ST500.reset(new LQAllhadHists(ctx, "LQMod_Hypothesis_ST500"));

  h_ele_ST600.reset(new ElectronHists(ctx, "LQMod_Electrons_ST600"));
  h_mu_ST600.reset(new MuonHists(ctx, "LQMod_Muons_ST600"));
  h_jet_ST600.reset(new JetHists(ctx, "LQMod_Jets_ST600"));
  h_topjet_ST600.reset(new TopJetHists(ctx, "LQMod_TopJets_ST600"));
  h_event_ST600.reset(new EventHists(ctx, "LQMod_Events_ST600"));
  h_tau_ST600.reset(new TauHists(ctx, "LQMod_Taus_ST600"));
  h_lq_ST600.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST600"));
  h_hypothesis_ST600.reset(new LQAllhadHists(ctx, "LQMod_Hypothesis_ST600"));


}


bool LQAnalysisModule::process(Event & event) {
  //cout << "LQAnalysisModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
    

  // 1. run all modules; here: only jet cleaning.
  bool pass_common = common->process(event);
  if(!pass_common) return false;
  jetcleaner->process(event);
  // HT calculator
  for (auto & mod : pre_modules) {
    mod->process(event);
  }


  if(!ntau_sel->passes(event)) return false;
  for (auto & m : recomodules) {
    m->process(event);
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

  if (!twojet_sel->passes(event)) return false;
  electron_PreSelection->fill(event);
  muon_PreSelection->fill(event);
  tau_PreSelection->fill(event);
  topjet_PreSelection->fill(event);
  jet_PreSelection->fill(event);
  event_PreSelection->fill(event);
  lq_PreSelection->fill(event);
  hypothesis_PreSelection->fill(event);

 

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

  //const auto taus = event.taus;
  const auto & tau1 = (*event.taus)[0];
  const auto & tau2 = (*event.taus)[1];
  if(tau1.pt()<50) return false;
  if(tau2.pt()<50) return false;

  if (!njet_sel->passes(event)) return false;
  h_lq_ThreeJets->fill(event);
  h_tau_ThreeJets->fill(event);
  h_mu_ThreeJets->fill(event);
  h_ele_ThreeJets->fill(event);
  h_jet_ThreeJets->fill(event);
  h_topjet_ThreeJets->fill(event);
  h_event_ThreeJets->fill(event);
  h_hypothesis_ThreeJets->fill(event);
 
 

  /*
  for(const auto & tau : *event.taus){
    if(tau.pt()<60) return false;
  }
  */
  
  if(!fourjet_sel->passes(event)) return false; // five jets!!

  //if(met<150) return false;
  //if(!BTagT->passes(event)) return false;
  if(!BTagM->passes(event)) return false;
  //if(!TwoBTagM->passes(event)) return false;
    
  //if(!BTagM->passes(event)) return false;
  //if(!BTagL->passes(event)) return false;
  //if(!toptag->passes(event)) return false;

  if(st<300) return false;
  h_lq_ST400->fill(event);
  h_tau_ST400->fill(event);
  h_mu_ST400->fill(event);
  h_ele_ST400->fill(event);
  h_jet_ST400->fill(event);
  h_topjet_ST400->fill(event);
  h_event_ST400->fill(event);
  h_hypothesis_ST400->fill(event);

  if(st<500) return false;
  h_lq_ST500->fill(event);
  h_tau_ST500->fill(event);
  h_mu_ST500->fill(event);
  h_ele_ST500->fill(event);
  h_jet_ST500->fill(event);
  h_topjet_ST500->fill(event);
  h_event_ST500->fill(event);
  h_hypothesis_ST500->fill(event);

  if(st<600) return false;
  h_lq_ST600->fill(event);
  h_tau_ST600->fill(event);
  h_mu_ST600->fill(event);
  h_ele_ST600->fill(event);
  h_jet_ST600->fill(event);
  h_topjet_ST600->fill(event);
  h_event_ST600->fill(event);
  h_hypothesis_ST600->fill(event);



  /*
    vector<Jet> bjets;
    for (unsigned int i =0; i<jets->size(); ++i) {
    if(jets->at(i).btag_combinedSecondaryVertex()>0.244) {
    bjets.push_back(jets->at(i));
    }
    }
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
    if((Tau+BJet).M()<200) return false;
    }
    }
    }
    }
  */

 

  //print all trigger names    
  /*for (unsigned int i=0; i<event.get_current_triggernames().size();i++)
    cout<< event.get_current_triggernames()[i]<<"\n";
    cout<<"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
  */

    


    
  // 3. decide whether or not to keep the current event in the output:
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisModule)
