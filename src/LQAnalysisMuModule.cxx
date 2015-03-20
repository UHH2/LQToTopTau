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
  std::unique_ptr<Selection> njet_sel, TwoBTagL, BTagM, BTagT, ntau_sel, ele_sel, muon_sel, onemuon_sel, twomuons_sel, isomuon_sel, leadingjet_sel, leadingjet300_sel, secondjet100_sel, samesign_sel;
  std::vector<std::unique_ptr<Selection> > fullhad_sel;

  std::vector<std::unique_ptr<AnalysisModule>> pre_modules;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> lq_PreSelection, electron_PreSelection, muon_PreSelection, jet_PreSelection, tau_PreSelection, event_PreSelection;
  std::unique_ptr<Hists> h_lq_LeadingJet150, h_tau_LeadingJet150, h_mu_LeadingJet150, h_ele_LeadingJet150, h_jet_LeadingJet150, h_event_LeadingJet150;
  std::unique_ptr<Hists> h_lq_ThreeJets, h_tau_ThreeJets, h_mu_ThreeJets, h_ele_ThreeJets, h_jet_ThreeJets, h_event_ThreeJets;
  std::unique_ptr<Hists> h_lq_ST400, h_tau_ST400, h_mu_ST400, h_ele_ST400, h_jet_ST400, h_event_ST400;
  std::unique_ptr<Hists> h_lq_ST500, h_tau_ST500, h_mu_ST500, h_ele_ST500, h_jet_ST500, h_event_ST500;
  std::unique_ptr<Hists> h_lq_ST700, h_tau_ST700, h_mu_ST700, h_ele_ST700, h_jet_ST700, h_event_ST700;
  std::unique_ptr<Hists> h_lq_ST900, h_tau_ST900, h_mu_ST900, h_ele_ST900, h_jet_ST900, h_event_ST900;
  std::unique_ptr<Hists> h_lq_Mareike, h_tau_Mareike, h_mu_Mareike, h_ele_Mareike, h_jet_Mareike, h_event_Mareike;
  std::unique_ptr<Hists> h_lq_FactorTwo, h_tau_FactorTwo, h_mu_FactorTwo, h_ele_FactorTwo, h_jet_FactorTwo, h_event_FactorTwo;
  std::unique_ptr<Hists> h_ele_full, h_tau_full, h_event_full, h_lq_full, h_jet_full, h_muon_full;
  std::unique_ptr<Hists> ele_HT1000, muon_HT1000, tau_HT1000, event_HT1000, jet_HT1000, lq_HT1000;
  std::unique_ptr<Hists> ele_HT1100MET100_BTagM, muon_HT1100MET100_BTagM, tau_HT1100MET100_BTagM, event_HT1100MET100_BTagM, jet_HT1100MET100_BTagM, lq_HT1100MET100_BTagM;
  std::unique_ptr<Hists> ele_HT1100MET100, muon_HT1100MET100, tau_HT1100MET100, event_HT1100MET100, jet_HT1100MET100, lq_HT1100MET100;
  std::unique_ptr<Hists> ele_HT1000_BTagM, muon_HT1000_BTagM, tau_HT1000_BTagM, event_HT1000_BTagM, jet_HT1000_BTagM, lq_HT1000_BTagM;
  std::unique_ptr<Hists> ele_HT1000_BTagT, muon_HT1000_BTagT, tau_HT1000_BTagT, event_HT1000_BTagT, jet_HT1000_BTagT, lq_HT1000_BTagT;
  std::unique_ptr<Hists> ele_HT1000_TwoBTagL, muon_HT1000_TwoBTagL, tau_HT1000_TwoBTagL, event_HT1000_TwoBTagL, jet_HT1000_TwoBTagL, lq_HT1000_TwoBTagL;
  std::unique_ptr<Hists> ele_PreSel_SameSign, muon_PreSel_SameSign, tau_PreSel_SameSign, event_PreSel_SameSign, jet_PreSel_SameSign, lq_PreSel_SameSign;

  JetId BTagLoose, BTagMedium, BTagTight;
  MuonId MuIso;


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
    muoncleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), MuonIDKinematic(30.0, 2.1))));
    electroncleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_PHYS14_25ns_medium, PtEtaCut(20.0, 2.5))));
    taucleaner.reset(new TauCleaner(AndId<Tau>(TauIDMedium(), PtEtaCut(20.0, 2.1))));
    BTagLoose = CSVBTag(CSVBTag::WP_LOOSE);
    BTagMedium = CSVBTag(CSVBTag::WP_MEDIUM);
    BTagTight = CSVBTag(CSVBTag::WP_TIGHT);
    MuIso = MuonIso(0.12);


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
    onemuon_sel.reset(new NMuonSelection(1,1));
    twomuons_sel.reset(new NMuonSelection(2,2));
    isomuon_sel.reset(new NMuonSelection(1,-1,MuIso));
    samesign_sel.reset(new SameSignCut());

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

    h_ele_ThreeJets.reset(new ElectronHists(ctx, "LQMod_Electrons_ThreeJets"));
    h_mu_ThreeJets.reset(new MuonHists(ctx, "LQMod_Muons_ThreeJets"));
    h_jet_ThreeJets.reset(new JetHists(ctx, "LQMod_Jets_ThreeJets"));
    h_event_ThreeJets.reset(new EventHists(ctx, "LQMod_Events_ThreeJets"));
    h_tau_ThreeJets.reset(new TauHists(ctx, "LQMod_Taus_ThreeJets"));
    h_lq_ThreeJets.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ThreeJets"));

    h_lq_LeadingJet150.reset(new LQAnalysisHists(ctx, "LQMod_LQ_LeadingJet150"));
    h_tau_LeadingJet150.reset(new TauHists(ctx, "LQMod_Taus_LeadingJet150"));
    h_mu_LeadingJet150.reset(new MuonHists(ctx, "LQMod_Muons_LeadingJet150"));
    h_ele_LeadingJet150.reset(new ElectronHists(ctx, "LQMod_Electrons_LeadingJet150"));
    h_jet_LeadingJet150.reset(new JetHists(ctx, "LQMod_Jets_LeadingJet150"));
    h_event_LeadingJet150.reset(new EventHists(ctx, "LQMod_Events_LeadingJet150"));

    h_ele_ST400.reset(new ElectronHists(ctx, "LQMod_Electrons_ST400"));
    h_mu_ST400.reset(new MuonHists(ctx, "LQMod_Muons_ST400"));
    h_jet_ST400.reset(new JetHists(ctx, "LQMod_Jets_ST400"));
    h_event_ST400.reset(new EventHists(ctx, "LQMod_Events_ST400"));
    h_tau_ST400.reset(new TauHists(ctx, "LQMod_Taus_ST400"));
    h_lq_ST400.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST400"));

    h_ele_ST500.reset(new ElectronHists(ctx, "LQMod_Electrons_ST500"));
    h_mu_ST500.reset(new MuonHists(ctx, "LQMod_Muons_ST500"));
    h_jet_ST500.reset(new JetHists(ctx, "LQMod_Jets_ST500"));
    h_event_ST500.reset(new EventHists(ctx, "LQMod_Events_ST500"));
    h_tau_ST500.reset(new TauHists(ctx, "LQMod_Taus_ST500"));
    h_lq_ST500.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST500"));

    h_ele_ST700.reset(new ElectronHists(ctx, "LQMod_Electrons_ST700"));
    h_mu_ST700.reset(new MuonHists(ctx, "LQMod_Muons_ST700"));
    h_jet_ST700.reset(new JetHists(ctx, "LQMod_Jets_ST700"));
    h_event_ST700.reset(new EventHists(ctx, "LQMod_Events_ST700"));
    h_tau_ST700.reset(new TauHists(ctx, "LQMod_Taus_ST700"));
    h_lq_ST700.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST700"));

    h_ele_ST900.reset(new ElectronHists(ctx, "LQMod_Electrons_ST900"));
    h_mu_ST900.reset(new MuonHists(ctx, "LQMod_Muons_ST900"));
    h_jet_ST900.reset(new JetHists(ctx, "LQMod_Jets_ST900"));
    h_event_ST900.reset(new EventHists(ctx, "LQMod_Events_ST900"));
    h_tau_ST900.reset(new TauHists(ctx, "LQMod_Taus_ST900"));
    h_lq_ST900.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST900"));

    h_lq_Mareike.reset(new LQAnalysisHists(ctx, "LQMod_LQ_Mareike"));
    h_tau_Mareike.reset(new TauHists(ctx, "LQMod_Taus_Mareike"));
    h_mu_Mareike.reset(new MuonHists(ctx, "LQMod_Muons_Mareike"));
    h_ele_Mareike.reset(new ElectronHists(ctx, "LQMod_Electrons_Mareike"));
    h_jet_Mareike.reset(new JetHists(ctx, "LQMod_Jets_Mareike"));
    h_event_Mareike.reset(new EventHists(ctx, "LQMod_Events_Mareike"));

    h_lq_FactorTwo.reset(new LQAnalysisHists(ctx, "LQMod_LQ_FactorTwo"));
    h_tau_FactorTwo.reset(new TauHists(ctx, "LQMod_Taus_FactorTwo"));
    h_mu_FactorTwo.reset(new MuonHists(ctx, "LQMod_Muons_FactorTwo"));
    h_ele_FactorTwo.reset(new ElectronHists(ctx, "LQMod_Electrons_FactorTwo"));
    h_jet_FactorTwo.reset(new JetHists(ctx, "LQMod_Jets_FactorTwo"));
    h_event_FactorTwo.reset(new EventHists(ctx, "LQMod_Events_FactorTwo"));

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

    ele_HT1100MET100_BTagM.reset(new ElectronHists(ctx, "ele_HT1100MET100_BTagM"));
    muon_HT1100MET100_BTagM.reset(new MuonHists(ctx, "muon_HT1100MET100_BTagM"));
    jet_HT1100MET100_BTagM.reset(new JetHists(ctx, "jet_HT1100MET100_BTagM"));
    event_HT1100MET100_BTagM.reset(new EventHists(ctx, "event_HT1100MET100_BTagM"));
    tau_HT1100MET100_BTagM.reset(new TauHists(ctx, "tau_HT1100MET100_BTagM"));
    lq_HT1100MET100_BTagM.reset(new LQAnalysisHists(ctx, "lq_HT1100MET100_BTagM"));

    ele_HT1100MET100.reset(new ElectronHists(ctx, "ele_HT1100MET100"));
    muon_HT1100MET100.reset(new MuonHists(ctx, "muon_HT1100MET100"));
    jet_HT1100MET100.reset(new JetHists(ctx, "jet_HT1100MET100"));
    event_HT1100MET100.reset(new EventHists(ctx, "event_HT1100MET100"));
    tau_HT1100MET100.reset(new TauHists(ctx, "tau_HT1100MET100"));
    lq_HT1100MET100.reset(new LQAnalysisHists(ctx, "lq_HT1100MET100"));

    ele_HT1000_BTagM.reset(new ElectronHists(ctx, "ele_HT1000_BTagM"));
    muon_HT1000_BTagM.reset(new MuonHists(ctx, "muon_HT1000_BTagM"));
    jet_HT1000_BTagM.reset(new JetHists(ctx, "jet_HT1000_BTagM"));
    event_HT1000_BTagM.reset(new EventHists(ctx, "event_HT1000_BTagM"));
    tau_HT1000_BTagM.reset(new TauHists(ctx, "tau_HT1000_BTagM"));
    lq_HT1000_BTagM.reset(new LQAnalysisHists(ctx, "lq_HT1000_BTagM"));

    ele_HT1000_BTagT.reset(new ElectronHists(ctx, "ele_HT1000_BTagT"));
    muon_HT1000_BTagT.reset(new MuonHists(ctx, "muon_HT1000_BTagT"));
    jet_HT1000_BTagT.reset(new JetHists(ctx, "jet_HT1000_BTagT"));
    event_HT1000_BTagT.reset(new EventHists(ctx, "event_HT1000_BTagT"));
    tau_HT1000_BTagT.reset(new TauHists(ctx, "tau_HT1000_BTagT"));
    lq_HT1000_BTagT.reset(new LQAnalysisHists(ctx, "lq_HT1000_BTagT"));

    ele_HT1000_TwoBTagL.reset(new ElectronHists(ctx, "ele_HT1000_TwoBTagL"));
    muon_HT1000_TwoBTagL.reset(new MuonHists(ctx, "muon_HT1000_TwoBTagL"));
    jet_HT1000_TwoBTagL.reset(new JetHists(ctx, "jet_HT1000_TwoBTagL"));
    event_HT1000_TwoBTagL.reset(new EventHists(ctx, "event_HT1000_TwoBTagL"));
    tau_HT1000_TwoBTagL.reset(new TauHists(ctx, "tau_HT1000_TwoBTagL"));
    lq_HT1000_TwoBTagL.reset(new LQAnalysisHists(ctx, "lq_HT1000_TwoBTagL"));

    ele_PreSel_SameSign.reset(new ElectronHists(ctx, "ele_PreSel_SameSign"));
    muon_PreSel_SameSign.reset(new MuonHists(ctx, "muon_PreSel_SameSign"));
    jet_PreSel_SameSign.reset(new JetHists(ctx, "jet_PreSel_SameSign"));
    event_PreSel_SameSign.reset(new EventHists(ctx, "event_PreSel_SameSign"));
    tau_PreSel_SameSign.reset(new TauHists(ctx, "tau_PreSel_SameSign"));
    lq_PreSel_SameSign.reset(new LQAnalysisHists(ctx, "lq_PreSel_SameSign"));

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
    
    if(!muon_sel->passes(event)) return false;
    
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

    electron_PreSelection->fill(event);
    muon_PreSelection->fill(event);
    tau_PreSelection->fill(event);
    jet_PreSelection->fill(event);
    event_PreSelection->fill(event);
    lq_PreSelection->fill(event);

    ele_PreSel_SameSign->fill(event);
    muon_PreSel_SameSign->fill(event);
    jet_PreSel_SameSign->fill(event);
    event_PreSel_SameSign->fill(event);
    tau_PreSel_SameSign->fill(event);
    lq_PreSel_SameSign->fill(event);

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


    if (!njet_sel->passes(event)) return false;
    h_lq_ThreeJets->fill(event);
    h_tau_ThreeJets->fill(event);
    h_mu_ThreeJets->fill(event);
    h_ele_ThreeJets->fill(event);
    h_jet_ThreeJets->fill(event);
    h_event_ThreeJets->fill(event);

    const auto jets = event.jets;
    if(jets->size() > 0){
      const auto & jet = (*jets)[0];
      if(jet.pt()<150) return false;
    }
    if(!leadingjet_sel->passes(event)) return false;
    h_lq_LeadingJet150->fill(event);
    h_tau_LeadingJet150->fill(event);
    h_mu_LeadingJet150->fill(event);
    h_ele_LeadingJet150->fill(event);
    h_jet_LeadingJet150->fill(event);
    h_event_LeadingJet150->fill(event);



    if(st<400) return false;
    h_lq_ST400->fill(event);
    h_tau_ST400->fill(event);
    h_mu_ST400->fill(event);
    h_ele_ST400->fill(event);
    h_jet_ST400->fill(event);
    h_event_ST400->fill(event);

    if(st<500) return false;
    h_lq_ST500->fill(event);
    h_tau_ST500->fill(event);
    h_mu_ST500->fill(event);
    h_ele_ST500->fill(event);
    h_jet_ST500->fill(event);
    h_event_ST500->fill(event);

    if(st<700) return false;
    h_lq_ST700->fill(event);
    h_tau_ST700->fill(event);
    h_mu_ST700->fill(event);
    h_ele_ST700->fill(event);
    h_jet_ST700->fill(event);
    h_event_ST700->fill(event);

    if(st<900) return false;
    h_lq_ST900->fill(event);
    h_tau_ST900->fill(event);
    h_mu_ST900->fill(event);
    h_ele_ST900->fill(event);
    h_jet_ST900->fill(event);
    h_event_ST900->fill(event);

    if(st<1000) return false;
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
    ele_HT1000->fill(event);
    muon_HT1000->fill(event);
    jet_HT1000->fill(event);
    event_HT1000->fill(event);
    tau_HT1000->fill(event);
    lq_HT1000->fill(event);


    //if (!BTagM->passes(event)) return false;



    if (BTagM->passes(event)){
      ele_HT1000_BTagM->fill(event);
      muon_HT1000_BTagM->fill(event);
      jet_HT1000_BTagM->fill(event);
      event_HT1000_BTagM->fill(event);
      tau_HT1000_BTagM->fill(event);
      lq_HT1000_BTagM->fill(event);
    }    

    if (TwoBTagL->passes(event)){
      ele_HT1000_TwoBTagL->fill(event);
      muon_HT1000_TwoBTagL->fill(event);
      jet_HT1000_TwoBTagL->fill(event);
      event_HT1000_TwoBTagL->fill(event);
      tau_HT1000_TwoBTagL->fill(event);
      lq_HT1000_TwoBTagL->fill(event);
    }
    
    if (BTagT->passes(event)){
      ele_HT1000_BTagT->fill(event);
      muon_HT1000_BTagT->fill(event);
      jet_HT1000_BTagT->fill(event);
      event_HT1000_BTagT->fill(event);
      tau_HT1000_BTagT->fill(event);
      lq_HT1000_BTagT->fill(event);
    }


    if(st<1100) return false;
    if(met<100) return false;
    //if(!leadingjet300_sel->passes(event)) return false;
    //if(!secondjet100_sel->passes(event)) return false;

    ele_HT1100MET100->fill(event);
    muon_HT1100MET100->fill(event);
    jet_HT1100MET100->fill(event);
    event_HT1100MET100->fill(event);
    tau_HT1100MET100->fill(event);
    lq_HT1100MET100->fill(event);
    

    if(BTagM->passes(event)){
      ele_HT1100MET100_BTagM->fill(event);
      muon_HT1100MET100_BTagM->fill(event);
      jet_HT1100MET100_BTagM->fill(event);
      event_HT1100MET100_BTagM->fill(event);
      tau_HT1100MET100_BTagM->fill(event);
      lq_HT1100MET100_BTagM->fill(event);
    }

    h_lq_FactorTwo->fill(event);
    h_tau_FactorTwo->fill(event);
    h_mu_FactorTwo->fill(event);
    h_ele_FactorTwo->fill(event);
    h_jet_FactorTwo->fill(event);
    h_event_FactorTwo->fill(event);


    
    // 3. decide whether or not to keep the current event in the output:
    return complete_selection;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisMuModule)
