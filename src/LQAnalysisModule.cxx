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
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/LQAnalysis/include/LQAnalysisSelections.h"
#include "UHH2/LQAnalysis/include/LQAnalysisHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/EventVariables.h"

#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/PrintingModules.h"

#include "TH1F.h"

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
  std::unique_ptr<Selection> njet_sel, fivejet_sel, TwoBTagL, BTagM, BTagT, TwoBTagM, TwoBTagT, TopTag, ntau_sel, muon_sel, ele_sel;
  //  std::vector<std::unique_ptr<Selection> > fullhad_sel;

  std::vector<std::unique_ptr<AnalysisModule>> pre_modules;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> lq_PreSelection, electron_PreSelection, muon_PreSelection, jet_PreSelection, tau_PreSelection, event_PreSelection;
  std::unique_ptr<Hists> lq_Selection, electron_Selection, muon_Selection, jet_Selection, tau_Selection, event_Selection;
  std::unique_ptr<Hists> lq_TwoBTagL, electron_TwoBTagL, muon_TwoBTagL, jet_TwoBTagL, tau_TwoBTagL, event_TwoBTagL;
  std::unique_ptr<Hists> lq_BTagM, electron_BTagM, muon_BTagM, jet_BTagM, tau_BTagM, event_BTagM;
  std::unique_ptr<Hists> lq_TwoBTagM, electron_TwoBTagM, muon_TwoBTagM, jet_TwoBTagM, tau_TwoBTagM, event_TwoBTagM;
  std::unique_ptr<Hists> lq_BTagT, electron_BTagT, muon_BTagT, jet_BTagT, tau_BTagT, event_BTagT;

  JetId BTagLoose, BTagMedium, BTagTight;
  MuonId muid;
  ElectronId eleid;
  TopJetId cmstoptag;

  std::unique_ptr<AnalysisModule> printer;
  std::unique_ptr<AnalysisModule> ttgenprod;
  Event::Handle<TTbarGen> h_ttbargen;

  TH1F* ttbardecay = new TH1F("TTbardecay","ttbardecay",11,-0.5,10.5);


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
    muoncleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), MuonIDKinematic(30.0, 2.1))));
    electroncleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_PHYS14_25ns_medium, PtEtaCut(30.0, 2.5))));
    taucleaner.reset(new TauCleaner(AndId<Tau>(TauIDMedium(), PtEtaCut(20.0, 2.1))));

    BTagLoose = CSVBTag(CSVBTag::WP_LOOSE);
    BTagMedium = CSVBTag(CSVBTag::WP_MEDIUM);
    BTagTight = CSVBTag(CSVBTag::WP_TIGHT);
    muid = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4));
    eleid = (AndId<Electron>(ElectronID_PHYS14_25ns_medium, PtEtaCut(20.0, 2.5)));
    cmstoptag = CMSTopTag(50,140,250);

    pre_modules.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));

    printer.reset(new GenParticlesPrinter(ctx));
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    
    // 2. set up selections:
    njet_sel.reset(new NJetSelection(6,-1));
    fivejet_sel.reset(new NJetSelection(5,-1));
    ntau_sel.reset(new NTauSelection(1,-1));
    TwoBTagL.reset(new NJetSelection(2,999,BTagLoose));
    BTagM.reset(new NJetSelection(1,999,BTagMedium));
    TwoBTagM.reset(new NJetSelection(2,999,BTagMedium));
    BTagT.reset(new NJetSelection(1,999,BTagTight));
    TwoBTagT.reset(new NJetSelection(2,999,BTagTight));
    TopTag.reset(new NTopJetSelection(1,999,cmstoptag));

    muon_sel.reset(new NMuonSelection(1,-1,muid));
    ele_sel.reset(new NElectronSelection(1,-1,eleid));



    // 3. Set up Hists classes:
    electron_PreSelection.reset(new ElectronHists(ctx, "LQMod_Electrons_PreSel"));
    muon_PreSelection.reset(new MuonHists(ctx, "LQMod_Muons_PreSel"));
    tau_PreSelection.reset(new TauHists(ctx, "LQMod_Taus_PreSel"));
    jet_PreSelection.reset(new JetHists(ctx, "LQMod_Jets_PreSel"));
    event_PreSelection.reset(new EventHists(ctx, "LQMod_Events_PreSel"));
    lq_PreSelection.reset(new LQAnalysisHists(ctx, "LQMod_LQ_PreSel"));

    electron_Selection.reset(new ElectronHists(ctx, "LQMod_Electrons_Sel"));
    muon_Selection.reset(new MuonHists(ctx, "LQMod_Muons_Sel"));
    tau_Selection.reset(new TauHists(ctx, "LQMod_Taus_Sel"));
    jet_Selection.reset(new JetHists(ctx, "LQMod_Jets_Sel"));
    event_Selection.reset(new EventHists(ctx, "LQMod_Events_Sel"));
    lq_Selection.reset(new LQAnalysisHists(ctx, "LQMod_LQ_Sel"));

    electron_TwoBTagL.reset(new ElectronHists(ctx, "LQMod_Electrons_TwoBTagL"));
    muon_TwoBTagL.reset(new MuonHists(ctx, "LQMod_Muons_TwoBTagL"));
    tau_TwoBTagL.reset(new TauHists(ctx, "LQMod_Taus_TwoBTagL"));
    jet_TwoBTagL.reset(new JetHists(ctx, "LQMod_Jets_TwoBTagL"));
    event_TwoBTagL.reset(new EventHists(ctx, "LQMod_Events_TwoBTagL"));
    lq_TwoBTagL.reset(new LQAnalysisHists(ctx, "LQMod_LQ_TwoBTagL"));

    electron_BTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_BTagM"));
    muon_BTagM.reset(new MuonHists(ctx, "LQMod_Muons_BTagM"));
    tau_BTagM.reset(new TauHists(ctx, "LQMod_Taus_BTagM"));
    jet_BTagM.reset(new JetHists(ctx, "LQMod_Jets_BTagM"));
    event_BTagM.reset(new EventHists(ctx, "LQMod_Events_BTagM"));
    lq_BTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_BTagM"));

    electron_TwoBTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_TwoBTagM"));
    muon_TwoBTagM.reset(new MuonHists(ctx, "LQMod_Muons_TwoBTagM"));
    tau_TwoBTagM.reset(new TauHists(ctx, "LQMod_Taus_TwoBTagM"));
    jet_TwoBTagM.reset(new JetHists(ctx, "LQMod_Jets_TwoBTagM"));
    event_TwoBTagM.reset(new EventHists(ctx, "LQMod_Events_TwoBTagM"));
    lq_TwoBTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_TwoBTagM"));

    electron_BTagT.reset(new ElectronHists(ctx, "LQMod_Electrons_BTagT"));
    muon_BTagT.reset(new MuonHists(ctx, "LQMod_Muons_BTagT"));
    tau_BTagT.reset(new TauHists(ctx, "LQMod_Taus_BTagT"));
    jet_BTagT.reset(new JetHists(ctx, "LQMod_Jets_BTagT"));
    event_BTagT.reset(new EventHists(ctx, "LQMod_Events_BTagT"));
    lq_BTagT.reset(new LQAnalysisHists(ctx, "LQMod_LQ_BTagT"));



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
    //jetcleaner->process(event);
    //muoncleaner->process(event);  // overlap with semileptonic channel
    //electroncleaner->process(event);  // overlap with semileptonic channel
    //taucleaner->process(event);

    for (auto & mod : pre_modules) {
      mod->process(event);
    }
    
    // 2. test selections and fill histograms
    
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

    // fill preselection hists
    electron_PreSelection->fill(event);
    muon_PreSelection->fill(event);
    tau_PreSelection->fill(event);
    jet_PreSelection->fill(event);
    event_PreSelection->fill(event);
    lq_PreSelection->fill(event);

    // declare met and st cut
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
 



    // leading jet pt > 100 GeV
    const auto jets = event.jets;
    if(jets->size() > 0){
      const auto & jet = (*jets)[0];
      if(jet.pt()<100) return false;
    }


    if (!njet_sel->passes(event)) return false;
    if(ht<800) return false;
    if(met<150) return false;


    electron_Selection->fill(event);
    muon_Selection->fill(event);
    tau_Selection->fill(event);
    jet_Selection->fill(event);
    event_Selection->fill(event);
    lq_Selection->fill(event);


    if(TwoBTagL->passes(event)){
      electron_TwoBTagL->fill(event);
      muon_TwoBTagL->fill(event);
      tau_TwoBTagL->fill(event);
      jet_TwoBTagL->fill(event);
      event_TwoBTagL->fill(event);
      lq_TwoBTagL->fill(event);
    }

    if(BTagM->passes(event)){
      electron_BTagM->fill(event);
      muon_BTagM->fill(event);
      tau_BTagM->fill(event);
      jet_BTagM->fill(event);
      event_BTagM->fill(event);
      lq_BTagM->fill(event);
    }

    if(TwoBTagM->passes(event)){
      electron_TwoBTagM->fill(event);
      muon_TwoBTagM->fill(event);
      tau_TwoBTagM->fill(event);
      jet_TwoBTagM->fill(event);
      event_TwoBTagM->fill(event);
      lq_TwoBTagM->fill(event);
    }

    if(BTagT->passes(event)){
      electron_BTagT->fill(event);
      muon_BTagT->fill(event);
      tau_BTagT->fill(event);
      jet_BTagT->fill(event);
      event_BTagT->fill(event);
      lq_BTagT->fill(event);
    }

    /*
    //printer->process(event);
    ttgenprod->process(event);
    const auto & ttbargen = event.get(h_ttbargen);
    //cout << "Decay channel is " << int(ttbargen.DecayChannel()) << endl;

    ttbardecay->Fill(ttbargen.DecayChannel());
    */

    
    // 3. decide whether or not to keep the current event in the output:
    return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisModule)
