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
  std::unique_ptr<Selection> njet_sel, fivejet_sel, BTagM, BTagT, TwoBTagM, TwoBTagT, TopTag, ntau_sel, muon_sel, ele_sel;
  std::vector<std::unique_ptr<Selection> > fullhad_sel;

  std::vector<std::unique_ptr<AnalysisModule>> pre_modules;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> lq_PreSelection, electron_PreSelection, muon_PreSelection, jet_PreSelection, tau_PreSelection, event_PreSelection;
  std::unique_ptr<Hists> h_ele_full, h_tau_full, h_event_full, h_lq_full, h_jet_full, h_muon_full;
  std::unique_ptr<Hists> ele_threejets, muon_threejets, tau_threejets, event_threejets, jet_threejets, lq_threejets;
  std::unique_ptr<Hists> ele_leadingjet, muon_leadingjet, tau_leadingjet, event_leadingjet, jet_leadingjet, lq_leadingjet;
  std::unique_ptr<Hists> ele_MET120, muon_MET120, tau_MET120, event_MET120, jet_MET120, lq_MET120;
  std::unique_ptr<Hists> ele_MET200, muon_MET200, tau_MET200, event_MET200, jet_MET200, lq_MET200;
  std::unique_ptr<Hists> ele_MET200HT1300, muon_MET200HT1300, tau_MET200HT1300, event_MET200HT1300, jet_MET200HT1300, lq_MET200HT1300;
  std::unique_ptr<Hists> ele_MET200HT1300_BTagT, muon_MET200HT1300_BTagT, tau_MET200HT1300_BTagT, event_MET200HT1300_BTagT, jet_MET200HT1300_BTagT, lq_MET200HT1300_BTagT;
  std::unique_ptr<Hists> ele_MET200HT1300_TwoBTagM, muon_MET200HT1300_TwoBTagM, tau_MET200HT1300_TwoBTagM, event_MET200HT1300_TwoBTagM, jet_MET200HT1300_TwoBTagM, lq_MET200HT1300_TwoBTagM;
  std::unique_ptr<Hists> ele_MET200LeadingTau150, muon_MET200LeadingTau150, tau_MET200LeadingTau150, event_MET200LeadingTau150, jet_MET200LeadingTau150, lq_MET200LeadingTau150;
  std::unique_ptr<Hists> ele_MET200LeadingTau150_BTagM, muon_MET200LeadingTau150_BTagM, tau_MET200LeadingTau150_BTagM, event_MET200LeadingTau150_BTagM, jet_MET200LeadingTau150_BTagM, lq_MET200LeadingTau150_BTagM;
  std::unique_ptr<Hists> ele_MET200LeadingTau150_TwoBTagM, muon_MET200LeadingTau150_TwoBTagM, tau_MET200LeadingTau150_TwoBTagM, event_MET200LeadingTau150_TwoBTagM, jet_MET200LeadingTau150_TwoBTagM, lq_MET200LeadingTau150_TwoBTagM;
  std::unique_ptr<Hists> ele_MET200LeadingTau150_BTagT, muon_MET200LeadingTau150_BTagT, tau_MET200LeadingTau150_BTagT, event_MET200LeadingTau150_BTagT, jet_MET200LeadingTau150_BTagT, lq_MET200LeadingTau150_BTagT;
  std::unique_ptr<Hists> ele_HT1000, muon_HT1000, tau_HT1000, event_HT1000, jet_HT1000, lq_HT1000;
  std::unique_ptr<Hists> ele_HT1000MET120, muon_HT1000MET120, tau_HT1000MET120, event_HT1000MET120, jet_HT1000MET120, lq_HT1000MET120;
  std::unique_ptr<Hists> ele_HT1000MET120_BTagT, muon_HT1000MET120_BTagT, tau_HT1000MET120_BTagT, event_HT1000MET120_BTagT, jet_HT1000MET120_BTagT, lq_HT1000MET120_BTagT;
  std::unique_ptr<Hists> ele_HT1000MET150, muon_HT1000MET150, tau_HT1000MET150, event_HT1000MET150, jet_HT1000MET150, lq_HT1000MET150;
  std::unique_ptr<Hists> ele_HT1500MET150, muon_HT1500MET150, tau_HT1500MET150, event_HT1500MET150, jet_HT1500MET150, lq_HT1500MET150;
  std::unique_ptr<Hists> ele_HT1500MET150_BTagM, muon_HT1500MET150_BTagM, tau_HT1500MET150_BTagM, event_HT1500MET150_BTagM, jet_HT1500MET150_BTagM, lq_HT1500MET150_BTagM;
  std::unique_ptr<Hists> ele_HT1500MET150_TwoBTagM, muon_HT1500MET150_TwoBTagM, tau_HT1500MET150_TwoBTagM, event_HT1500MET150_TwoBTagM, jet_HT1500MET150_TwoBTagM, lq_HT1500MET150_TwoBTagM;
  std::unique_ptr<Hists> ele_HT1600MET150_TwoBTagM, muon_HT1600MET150_TwoBTagM, tau_HT1600MET150_TwoBTagM, event_HT1600MET150_TwoBTagM, jet_HT1600MET150_TwoBTagM, lq_HT1600MET150_TwoBTagM;
  std::unique_ptr<Hists> ele_HT1800MET150_TwoBTagM, muon_HT1800MET150_TwoBTagM, tau_HT1800MET150_TwoBTagM, event_HT1800MET150_TwoBTagM, jet_HT1800MET150_TwoBTagM, lq_HT1800MET150_TwoBTagM;
  std::unique_ptr<Hists> ele_HT1000MET150_BTagM, muon_HT1000MET150_BTagM, tau_HT1000MET150_BTagM, event_HT1000MET150_BTagM, jet_HT1000MET150_BTagM, lq_HT1000MET150_BTagM;
  std::unique_ptr<Hists> ele_HT1000MET150_TwoBTagM, muon_HT1000MET150_TwoBTagM, tau_HT1000MET150_TwoBTagM, event_HT1000MET150_TwoBTagM, jet_HT1000MET150_TwoBTagM, lq_HT1000MET150_TwoBTagM;

  JetId BTagMedium, BTagTight;
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
    njet_sel.reset(new NJetSelection(5,-1));
    fivejet_sel.reset(new NJetSelection(5,-1));
    ntau_sel.reset(new NTauSelection(1,-1));
    BTagM.reset(new NJetSelection(1,999,BTagMedium));
    TwoBTagM.reset(new NJetSelection(2,999,BTagMedium));
    BTagT.reset(new NJetSelection(1,999,BTagTight));
    TwoBTagT.reset(new NJetSelection(2,999,BTagTight));
    TopTag.reset(new NTopJetSelection(1,999,cmstoptag));

    muon_sel.reset(new NMuonSelection(1,-1,muid));
    ele_sel.reset(new NElectronSelection(1,-1,eleid));

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
    lq_PreSelection.reset(new LQAnalysisHists(ctx, "LQ_PreSel"));

    h_ele_full.reset(new ElectronHists(ctx, "LQMod_Electrons_fullselection"));
    h_muon_full.reset(new MuonHists(ctx, "LQMod_Muons_fullselection"));
    h_jet_full.reset(new JetHists(ctx, "LQMod_Jets_fullselection"));
    h_event_full.reset(new EventHists(ctx, "LQMod_Events_fullselection"));
    h_tau_full.reset(new TauHists(ctx, "LQMod_Taus_fullselection"));
    h_lq_full.reset(new LQAnalysisHists(ctx, "LQMod_LQ_fullselection"));

    ele_threejets.reset(new ElectronHists(ctx, "LQMod_Electrons_ThreeJets"));
    muon_threejets.reset(new MuonHists(ctx, "LQMod_Muons_ThreeJets"));
    jet_threejets.reset(new JetHists(ctx, "LQMod_Jets_ThreeJets"));
    event_threejets.reset(new EventHists(ctx, "LQMod_Events_ThreeJets"));
    tau_threejets.reset(new TauHists(ctx, "LQMod_Taus_ThreeJets"));
    lq_threejets.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ThreeJets"));

    ele_leadingjet.reset(new ElectronHists(ctx, "LQMod_Electrons_LeadingJet"));
    muon_leadingjet.reset(new MuonHists(ctx, "LQMod_Muons_LeadingJet"));
    jet_leadingjet.reset(new JetHists(ctx, "LQMod_Jets_LeadingJet"));
    event_leadingjet.reset(new EventHists(ctx, "LQMod_Events_LeadingJet"));
    tau_leadingjet.reset(new TauHists(ctx, "LQMod_Taus_LeadingJet"));
    lq_leadingjet.reset(new LQAnalysisHists(ctx, "LQMod_LQ_LeadingJet"));

    ele_MET120.reset(new ElectronHists(ctx, "LQMod_Electrons_MET120"));
    muon_MET120.reset(new MuonHists(ctx, "LQMod_Muons_MET120"));
    jet_MET120.reset(new JetHists(ctx, "LQMod_Jets_MET120"));
    event_MET120.reset(new EventHists(ctx, "LQMod_Events_MET120"));
    tau_MET120.reset(new TauHists(ctx, "LQMod_Taus_MET120"));
    lq_MET120.reset(new LQAnalysisHists(ctx, "LQMod_LQ_MET120"));

    ele_MET200.reset(new ElectronHists(ctx, "LQMod_Electrons_MET200"));
    muon_MET200.reset(new MuonHists(ctx, "LQMod_Muons_MET200"));
    jet_MET200.reset(new JetHists(ctx, "LQMod_Jets_MET200"));
    event_MET200.reset(new EventHists(ctx, "LQMod_Events_MET200"));
    tau_MET200.reset(new TauHists(ctx, "LQMod_Taus_MET200"));
    lq_MET200.reset(new LQAnalysisHists(ctx, "LQMod_LQ_MET200"));

    ele_MET200HT1300.reset(new ElectronHists(ctx, "LQMod_Electrons_MET200HT1300"));
    muon_MET200HT1300.reset(new MuonHists(ctx, "LQMod_Muons_MET200HT1300"));
    jet_MET200HT1300.reset(new JetHists(ctx, "LQMod_Jets_MET200HT1300"));
    event_MET200HT1300.reset(new EventHists(ctx, "LQMod_Events_MET200HT1300"));
    tau_MET200HT1300.reset(new TauHists(ctx, "LQMod_Taus_MET200HT1300"));
    lq_MET200HT1300.reset(new LQAnalysisHists(ctx, "LQMod_LQ_MET200HT1300"));

    ele_MET200HT1300_TwoBTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_MET200HT1300_TwoBTagM"));
    muon_MET200HT1300_TwoBTagM.reset(new MuonHists(ctx, "LQMod_Muons_MET200HT1300_TwoBTagM"));
    jet_MET200HT1300_TwoBTagM.reset(new JetHists(ctx, "LQMod_Jets_MET200HT1300_TwoBTagM"));
    event_MET200HT1300_TwoBTagM.reset(new EventHists(ctx, "LQMod_Events_MET200HT1300_TwoBTagM"));
    tau_MET200HT1300_TwoBTagM.reset(new TauHists(ctx, "LQMod_Taus_MET200HT1300_TwoBTagM"));
    lq_MET200HT1300_TwoBTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_MET200HT1300_TwoBTagM"));

    ele_MET200HT1300_BTagT.reset(new ElectronHists(ctx, "LQMod_Electrons_MET200HT1300_BTagT"));
    muon_MET200HT1300_BTagT.reset(new MuonHists(ctx, "LQMod_Muons_MET200HT1300_BTagT"));
    jet_MET200HT1300_BTagT.reset(new JetHists(ctx, "LQMod_Jets_MET200HT1300_BTagT"));
    event_MET200HT1300_BTagT.reset(new EventHists(ctx, "LQMod_Events_MET200HT1300_BTagT"));
    tau_MET200HT1300_BTagT.reset(new TauHists(ctx, "LQMod_Taus_MET200HT1300_BTagT"));
    lq_MET200HT1300_BTagT.reset(new LQAnalysisHists(ctx, "LQMod_LQ_MET200HT1300_BTagT"));

    ele_MET200LeadingTau150.reset(new ElectronHists(ctx, "LQMod_Electrons_MET200LeadingTau150"));
    muon_MET200LeadingTau150.reset(new MuonHists(ctx, "LQMod_Muons_MET200LeadingTau150"));
    jet_MET200LeadingTau150.reset(new JetHists(ctx, "LQMod_Jets_MET200LeadingTau150"));
    event_MET200LeadingTau150.reset(new EventHists(ctx, "LQMod_Events_MET200LeadingTau150"));
    tau_MET200LeadingTau150.reset(new TauHists(ctx, "LQMod_Taus_MET200LeadingTau150"));
    lq_MET200LeadingTau150.reset(new LQAnalysisHists(ctx, "LQMod_LQ_MET200LeadingTau150"));

    ele_MET200LeadingTau150_BTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_MET200LeadingTau150_BTagM"));
    muon_MET200LeadingTau150_BTagM.reset(new MuonHists(ctx, "LQMod_Muons_MET200LeadingTau150_BTagM"));
    jet_MET200LeadingTau150_BTagM.reset(new JetHists(ctx, "LQMod_Jets_MET200LeadingTau150_BTagM"));
    event_MET200LeadingTau150_BTagM.reset(new EventHists(ctx, "LQMod_Events_MET200LeadingTau150_BTagM"));
    tau_MET200LeadingTau150_BTagM.reset(new TauHists(ctx, "LQMod_Taus_MET200LeadingTau150_BTagM"));
    lq_MET200LeadingTau150_BTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_MET200LeadingTau150_BTagM"));

    ele_MET200LeadingTau150_TwoBTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_MET200LeadingTau150_TwoBTagM"));
    muon_MET200LeadingTau150_TwoBTagM.reset(new MuonHists(ctx, "LQMod_Muons_MET200LeadingTau150_TwoBTagM"));
    jet_MET200LeadingTau150_TwoBTagM.reset(new JetHists(ctx, "LQMod_Jets_MET200LeadingTau150_TwoBTagM"));
    event_MET200LeadingTau150_TwoBTagM.reset(new EventHists(ctx, "LQMod_Events_MET200LeadingTau150_TwoBTagM"));
    tau_MET200LeadingTau150_TwoBTagM.reset(new TauHists(ctx, "LQMod_Taus_MET200LeadingTau150_TwoBTagM"));
    lq_MET200LeadingTau150_TwoBTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_MET200LeadingTau150_TwoBTagM"));

    ele_MET200LeadingTau150_BTagT.reset(new ElectronHists(ctx, "LQMod_Electrons_MET200LeadingTau150_BTagT"));
    muon_MET200LeadingTau150_BTagT.reset(new MuonHists(ctx, "LQMod_Muons_MET200LeadingTau150_BTagT"));
    jet_MET200LeadingTau150_BTagT.reset(new JetHists(ctx, "LQMod_Jets_MET200LeadingTau150_BTagT"));
    event_MET200LeadingTau150_BTagT.reset(new EventHists(ctx, "LQMod_Events_MET200LeadingTau150_BTagT"));
    tau_MET200LeadingTau150_BTagT.reset(new TauHists(ctx, "LQMod_Taus_MET200LeadingTau150_BTagT"));
    lq_MET200LeadingTau150_BTagT.reset(new LQAnalysisHists(ctx, "LQMod_LQ_MET200LeadingTau150_BTagT"));


    ele_HT1000.reset(new ElectronHists(ctx, "LQMod_Electrons_HT1000"));
    muon_HT1000.reset(new MuonHists(ctx, "LQMod_Muons_HT1000"));
    jet_HT1000.reset(new JetHists(ctx, "LQMod_Jets_HT1000"));
    event_HT1000.reset(new EventHists(ctx, "LQMod_Events_HT1000"));
    tau_HT1000.reset(new TauHists(ctx, "LQMod_Taus_HT1000"));
    lq_HT1000.reset(new LQAnalysisHists(ctx, "LQMod_LQ_HT1000"));

    ele_HT1000MET120.reset(new ElectronHists(ctx, "LQMod_Electrons_HT1000MET120"));
    muon_HT1000MET120.reset(new MuonHists(ctx, "LQMod_Muons_HT1000MET120"));
    jet_HT1000MET120.reset(new JetHists(ctx, "LQMod_Jets_HT1000MET120"));
    event_HT1000MET120.reset(new EventHists(ctx, "LQMod_Events_HT1000MET120"));
    tau_HT1000MET120.reset(new TauHists(ctx, "LQMod_Taus_HT1000MET120"));
    lq_HT1000MET120.reset(new LQAnalysisHists(ctx, "LQMod_LQ_HT1000MET120"));

    ele_HT1000MET120_BTagT.reset(new ElectronHists(ctx, "LQMod_Electrons_HT1000MET120_BTagT"));
    muon_HT1000MET120_BTagT.reset(new MuonHists(ctx, "LQMod_Muons_HT1000MET120_BTagT"));
    jet_HT1000MET120_BTagT.reset(new JetHists(ctx, "LQMod_Jets_HT1000MET120_BTagT"));
    event_HT1000MET120_BTagT.reset(new EventHists(ctx, "LQMod_Events_HT1000MET120_BTagT"));
    tau_HT1000MET120_BTagT.reset(new TauHists(ctx, "LQMod_Taus_HT1000MET120_BTagT"));
    lq_HT1000MET120_BTagT.reset(new LQAnalysisHists(ctx, "LQMod_LQ_HT1000MET120_BTagT"));

    ele_HT1000MET150.reset(new ElectronHists(ctx, "LQMod_Electrons_HT1000MET150"));
    muon_HT1000MET150.reset(new MuonHists(ctx, "LQMod_Muons_HT1000MET150"));
    jet_HT1000MET150.reset(new JetHists(ctx, "LQMod_Jets_HT1000MET150"));
    event_HT1000MET150.reset(new EventHists(ctx, "LQMod_Events_HT1000MET150"));
    tau_HT1000MET150.reset(new TauHists(ctx, "LQMod_Taus_HT1000MET150"));
    lq_HT1000MET150.reset(new LQAnalysisHists(ctx, "LQMod_LQ_HT1000MET150"));

    ele_HT1500MET150.reset(new ElectronHists(ctx, "LQMod_Electrons_HT1500MET150"));
    muon_HT1500MET150.reset(new MuonHists(ctx, "LQMod_Muons_HT1500MET150"));
    jet_HT1500MET150.reset(new JetHists(ctx, "LQMod_Jets_HT1500MET150"));
    event_HT1500MET150.reset(new EventHists(ctx, "LQMod_Events_HT1500MET150"));
    tau_HT1500MET150.reset(new TauHists(ctx, "LQMod_Taus_HT1500MET150"));
    lq_HT1500MET150.reset(new LQAnalysisHists(ctx, "LQMod_LQ_HT1500MET150"));

    ele_HT1500MET150_BTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_HT1500MET150_BTagM"));
    muon_HT1500MET150_BTagM.reset(new MuonHists(ctx, "LQMod_Muons_HT1500MET150_BTagM"));
    jet_HT1500MET150_BTagM.reset(new JetHists(ctx, "LQMod_Jets_HT1500MET150_BTagM"));
    event_HT1500MET150_BTagM.reset(new EventHists(ctx, "LQMod_Events_HT1500MET150_BTagM"));
    tau_HT1500MET150_BTagM.reset(new TauHists(ctx, "LQMod_Taus_HT1500MET150_BTagM"));
    lq_HT1500MET150_BTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_HT1500MET150_BTagM"));

    ele_HT1500MET150_TwoBTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_HT1500MET150_TwoBTagM"));
    muon_HT1500MET150_TwoBTagM.reset(new MuonHists(ctx, "LQMod_Muons_HT1500MET150_TwoBTagM"));
    jet_HT1500MET150_TwoBTagM.reset(new JetHists(ctx, "LQMod_Jets_HT1500MET150_TwoBTagM"));
    event_HT1500MET150_TwoBTagM.reset(new EventHists(ctx, "LQMod_Events_HT1500MET150_TwoBTagM"));
    tau_HT1500MET150_TwoBTagM.reset(new TauHists(ctx, "LQMod_Taus_HT1500MET150_TwoBTagM"));
    lq_HT1500MET150_TwoBTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_HT1500MET150_TwoBTagM"));

    ele_HT1600MET150_TwoBTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_HT1600MET150_TwoBTagM"));
    muon_HT1600MET150_TwoBTagM.reset(new MuonHists(ctx, "LQMod_Muons_HT1600MET150_TwoBTagM"));
    jet_HT1600MET150_TwoBTagM.reset(new JetHists(ctx, "LQMod_Jets_HT1600MET150_TwoBTagM"));
    event_HT1600MET150_TwoBTagM.reset(new EventHists(ctx, "LQMod_Events_HT1600MET150_TwoBTagM"));
    tau_HT1600MET150_TwoBTagM.reset(new TauHists(ctx, "LQMod_Taus_HT1600MET150_TwoBTagM"));
    lq_HT1600MET150_TwoBTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_HT1600MET150_TwoBTagM"));

    ele_HT1800MET150_TwoBTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_HT1800MET150_TwoBTagM"));
    muon_HT1800MET150_TwoBTagM.reset(new MuonHists(ctx, "LQMod_Muons_HT1800MET150_TwoBTagM"));
    jet_HT1800MET150_TwoBTagM.reset(new JetHists(ctx, "LQMod_Jets_HT1800MET150_TwoBTagM"));
    event_HT1800MET150_TwoBTagM.reset(new EventHists(ctx, "LQMod_Events_HT1800MET150_TwoBTagM"));
    tau_HT1800MET150_TwoBTagM.reset(new TauHists(ctx, "LQMod_Taus_HT1800MET150_TwoBTagM"));
    lq_HT1800MET150_TwoBTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_HT1800MET150_TwoBTagM"));

    ele_HT1000MET150_BTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_HT1000MET150_BTagM"));
    muon_HT1000MET150_BTagM.reset(new MuonHists(ctx, "LQMod_Muons_HT1000MET150_BTagM"));
    jet_HT1000MET150_BTagM.reset(new JetHists(ctx, "LQMod_Jets_HT1000MET150_BTagM"));
    event_HT1000MET150_BTagM.reset(new EventHists(ctx, "LQMod_Events_HT1000MET150_BTagM"));
    tau_HT1000MET150_BTagM.reset(new TauHists(ctx, "LQMod_Taus_HT1000MET150_BTagM"));
    lq_HT1000MET150_BTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_HT1000MET150_BTagM"));

    ele_HT1000MET150_TwoBTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_HT1000MET150_TwoBTagM"));
    muon_HT1000MET150_TwoBTagM.reset(new MuonHists(ctx, "LQMod_Muons_HT1000MET150_TwoBTagM"));
    jet_HT1000MET150_TwoBTagM.reset(new JetHists(ctx, "LQMod_Jets_HT1000MET150_TwoBTagM"));
    event_HT1000MET150_TwoBTagM.reset(new EventHists(ctx, "LQMod_Events_HT1000MET150_TwoBTagM"));
    tau_HT1000MET150_TwoBTagM.reset(new TauHists(ctx, "LQMod_Taus_HT1000MET150_TwoBTagM"));
    lq_HT1000MET150_TwoBTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_HT1000MET150_TwoBTagM"));

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

    if (!njet_sel->passes(event)) return false;
    ele_threejets->fill(event);
    muon_threejets->fill(event);
    jet_threejets->fill(event);
    event_threejets->fill(event);
    tau_threejets->fill(event);
    lq_threejets->fill(event);



    // leading jet pt > 100 GeV
    const auto jets = event.jets;
    if(jets->size() > 0){
      const auto & jet = (*jets)[0];
      if(jet.pt()<100) return false;
    }


    ele_leadingjet->fill(event);
    muon_leadingjet->fill(event);
    jet_leadingjet->fill(event);
    event_leadingjet->fill(event);
    tau_leadingjet->fill(event);
    lq_leadingjet->fill(event);
  


   if(met>120){
      ele_MET120->fill(event);
      muon_MET120->fill(event);
      jet_MET120->fill(event);
      event_MET120->fill(event);
      tau_MET120->fill(event);
      lq_MET120->fill(event);
    }
    
    if (st>1000){
      ele_HT1000->fill(event);
      muon_HT1000->fill(event);
      jet_HT1000->fill(event);
      event_HT1000->fill(event);
      tau_HT1000->fill(event);
      lq_HT1000->fill(event);
    }

 
    //met 200
    if(met>200){
      ele_MET200->fill(event);
      muon_MET200->fill(event);
      jet_MET200->fill(event);
      event_MET200->fill(event);
      tau_MET200->fill(event);
      lq_MET200->fill(event);
      if(st>1300){
	ele_MET200HT1300->fill(event);
	muon_MET200HT1300->fill(event);
	jet_MET200HT1300->fill(event);
	event_MET200HT1300->fill(event);
	tau_MET200HT1300->fill(event);
	lq_MET200HT1300->fill(event);
	if(TwoBTagM->passes(event)){
	  ele_MET200HT1300_TwoBTagM->fill(event);
	  muon_MET200HT1300_TwoBTagM->fill(event);
	  jet_MET200HT1300_TwoBTagM->fill(event);
	  event_MET200HT1300_TwoBTagM->fill(event);
	  tau_MET200HT1300_TwoBTagM->fill(event);
	  lq_MET200HT1300_TwoBTagM->fill(event);
	}
	if(BTagT->passes(event)){
	  ele_MET200HT1300_BTagT->fill(event);
	  muon_MET200HT1300_BTagT->fill(event);
	  jet_MET200HT1300_BTagT->fill(event);
	  event_MET200HT1300_BTagT->fill(event);
	  tau_MET200HT1300_BTagT->fill(event);
	  lq_MET200HT1300_BTagT->fill(event);
	}
      }
    }

    const auto taus = event.taus;
    if(fivejet_sel->passes(event)){
      if(met>200){
	if(taus->size() > 0){
	  const auto & tau = (*taus)[0];
	  if(tau.pt()>150){
	    ele_MET200LeadingTau150->fill(event);
	    muon_MET200LeadingTau150->fill(event);
	    jet_MET200LeadingTau150->fill(event);
	    event_MET200LeadingTau150->fill(event);
	    tau_MET200LeadingTau150->fill(event);
	    lq_MET200LeadingTau150->fill(event);
	    if(BTagM->passes(event)){
	      ele_MET200LeadingTau150_BTagM->fill(event);
	      muon_MET200LeadingTau150_BTagM->fill(event);
	      jet_MET200LeadingTau150_BTagM->fill(event);
	      event_MET200LeadingTau150_BTagM->fill(event);
	      tau_MET200LeadingTau150_BTagM->fill(event);
	      lq_MET200LeadingTau150_BTagM->fill(event);
	    }
	    if(TwoBTagM->passes(event)){
	      ele_MET200LeadingTau150_TwoBTagM->fill(event);
	      muon_MET200LeadingTau150_TwoBTagM->fill(event);
	      jet_MET200LeadingTau150_TwoBTagM->fill(event);
	      event_MET200LeadingTau150_TwoBTagM->fill(event);
	      tau_MET200LeadingTau150_TwoBTagM->fill(event);
	      lq_MET200LeadingTau150_TwoBTagM->fill(event);
	    }
	    if(BTagT->passes(event)){
	      ele_MET200LeadingTau150_BTagT->fill(event);
	      muon_MET200LeadingTau150_BTagT->fill(event);
	      jet_MET200LeadingTau150_BTagT->fill(event);
	      event_MET200LeadingTau150_BTagT->fill(event);
	      tau_MET200LeadingTau150_BTagT->fill(event);
	      lq_MET200LeadingTau150_BTagT->fill(event);
	    }
	  }
	}
      }
    }

    



    if(met<120) return false;
    if(st<1000) return false;
    
    ele_HT1000MET120->fill(event);
    muon_HT1000MET120->fill(event);
    jet_HT1000MET120->fill(event);
    event_HT1000MET120->fill(event);
    tau_HT1000MET120->fill(event);
    lq_HT1000MET120->fill(event);

    if(BTagT->passes(event)){
      ele_HT1000MET120_BTagT->fill(event);
      muon_HT1000MET120_BTagT->fill(event);
      jet_HT1000MET120_BTagT->fill(event);
      event_HT1000MET120_BTagT->fill(event);
      tau_HT1000MET120_BTagT->fill(event);
      lq_HT1000MET120_BTagT->fill(event);
    }

    // met 150
    if(met<150) return false;
    ele_HT1000MET150->fill(event);
    muon_HT1000MET150->fill(event);
    jet_HT1000MET150->fill(event);
    event_HT1000MET150->fill(event);
    tau_HT1000MET150->fill(event);
    lq_HT1000MET150->fill(event);


    if(BTagM->passes(event)){
      ele_HT1000MET150_BTagM->fill(event);
      muon_HT1000MET150_BTagM->fill(event);
      jet_HT1000MET150_BTagM->fill(event);
      event_HT1000MET150_BTagM->fill(event);
      tau_HT1000MET150_BTagM->fill(event);
      lq_HT1000MET150_BTagM->fill(event);
    }

    if(TwoBTagM->passes(event)){
      ele_HT1000MET150_TwoBTagM->fill(event);
      muon_HT1000MET150_TwoBTagM->fill(event);
      jet_HT1000MET150_TwoBTagM->fill(event);
      event_HT1000MET150_TwoBTagM->fill(event);
      tau_HT1000MET150_TwoBTagM->fill(event);
      lq_HT1000MET150_TwoBTagM->fill(event);
    }
    
    if(st<1500) return false;
    ele_HT1500MET150->fill(event);
    muon_HT1500MET150->fill(event);
    jet_HT1500MET150->fill(event);
    event_HT1500MET150->fill(event);
    tau_HT1500MET150->fill(event);
    lq_HT1500MET150->fill(event);
    if(BTagM->passes(event)){
      ele_HT1500MET150_BTagM->fill(event);
      muon_HT1500MET150_BTagM->fill(event);
      jet_HT1500MET150_BTagM->fill(event);
      event_HT1500MET150_BTagM->fill(event);
      tau_HT1500MET150_BTagM->fill(event);
      lq_HT1500MET150_BTagM->fill(event);
    }
    if(TwoBTagM->passes(event)){
      ele_HT1500MET150_TwoBTagM->fill(event);
      muon_HT1500MET150_TwoBTagM->fill(event);
      jet_HT1500MET150_TwoBTagM->fill(event);
      event_HT1500MET150_TwoBTagM->fill(event);
      tau_HT1500MET150_TwoBTagM->fill(event);
      lq_HT1500MET150_TwoBTagM->fill(event);
    }

    if(st<1600) return false;
    if(TwoBTagM->passes(event)){
      ele_HT1600MET150_TwoBTagM->fill(event);
      muon_HT1600MET150_TwoBTagM->fill(event);
      jet_HT1600MET150_TwoBTagM->fill(event);
      event_HT1600MET150_TwoBTagM->fill(event);
      tau_HT1600MET150_TwoBTagM->fill(event);
      lq_HT1600MET150_TwoBTagM->fill(event);
    }
    if(st<1800) return false;
    if(TwoBTagM->passes(event)){
      ele_HT1800MET150_TwoBTagM->fill(event);
      muon_HT1800MET150_TwoBTagM->fill(event);
      jet_HT1800MET150_TwoBTagM->fill(event);
      event_HT1800MET150_TwoBTagM->fill(event);
      tau_HT1800MET150_TwoBTagM->fill(event);
      lq_HT1800MET150_TwoBTagM->fill(event);
    }


    /*
    //printer->process(event);
    ttgenprod->process(event);
    const auto & ttbargen = event.get(h_ttbargen);
    //cout << "Decay channel is " << int(ttbargen.DecayChannel()) << endl;

    ttbardecay->Fill(ttbargen.DecayChannel());
    */

    
    // 3. decide whether or not to keep the current event in the output:
    return complete_selection;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisModule)
