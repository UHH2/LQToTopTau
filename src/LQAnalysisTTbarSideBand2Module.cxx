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
#include "UHH2/LQAnalysis/include/LQFakeTauHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/PrintingModules.h"



using namespace std;
using namespace uhh2;

/** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
 * 
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class LQAnalysisTTbarSideBand2Module: public AnalysisModule {
public:
    
    explicit LQAnalysisTTbarSideBand2Module(Context & ctx);
    virtual bool process(Event & event) override;

private:
    
    
  std::unique_ptr<JetCleaner> jetcleaner;
  std::unique_ptr<MuonIDKinematic> muonidkinematic;
  std::unique_ptr<MuonCleaner> muoncleaner;
  std::unique_ptr<MuonCleaner> muoncleaner_iso;
  std::unique_ptr<TauCleaner> taucleaner;
  std::unique_ptr<ElectronCleaner> electroncleaner;
  std::unique_ptr<ElectronCleaner> electroncleaner_iso;
  std::unique_ptr<JetLeptonCleaner> jetleptoncleaner;
  
  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  std::unique_ptr<Selection> njet_sel, threejet_sel, fourjet_sel, fivejet_sel, bsel, ntau_sel, nmuon_sel, TwoBTagL, TwoBTagM, BTagM, BTagT, TwoBTagT, samesign_sel, samesign_lead;
  std::vector<std::unique_ptr<Selection> > fullhad_sel;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  /*std::unique_ptr<Hists> h_lq_PreSel, h_tau_PreSel, h_mu_PreSel, h_ele_PreSel, h_jet_PreSel, h_event_PreSel;
  std::unique_ptr<Hists> h_lq_FourJets, h_tau_FourJets, h_mu_FourJets, h_ele_FourJets, h_jet_FourJets, h_event_FourJets;
  std::unique_ptr<Hists> h_lq_FourJets_SameSign, h_tau_FourJets_SameSign, h_mu_FourJets_SameSign, h_ele_FourJets_SameSign, h_jet_FourJets_SameSign, h_event_FourJets_SameSign;
  std::unique_ptr<Hists> h_lq_FourJets_OppositeSign, h_tau_FourJets_OppositeSign, h_mu_FourJets_OppositeSign, h_ele_FourJets_OppositeSign, h_jet_FourJets_OppositeSign, h_event_FourJets_OppositeSign;
  std::unique_ptr<Hists> h_lq_FiveJets, h_tau_FiveJets, h_mu_FiveJets, h_ele_FiveJets, h_jet_FiveJets, h_event_FiveJets;
  std::unique_ptr<Hists> h_lq_FiveJets_TwoBTagM, h_tau_FiveJets_TwoBTagM, h_mu_FiveJets_TwoBTagM, h_ele_FiveJets_TwoBTagM, h_jet_FiveJets_TwoBTagM, h_event_FiveJets_TwoBTagM;
  std::unique_ptr<Hists> h_lq_FiveJets_TwoBTagM_SameSign, h_tau_FiveJets_TwoBTagM_SameSign, h_mu_FiveJets_TwoBTagM_SameSign, h_ele_FiveJets_TwoBTagM_SameSign, h_jet_FiveJets_TwoBTagM_SameSign, h_event_FiveJets_TwoBTagM_SameSign;
  std::unique_ptr<Hists> h_lq_FiveJets_TwoBTagM_OppositeSign, h_tau_FiveJets_TwoBTagM_OppositeSign, h_mu_FiveJets_TwoBTagM_OppositeSign, h_ele_FiveJets_TwoBTagM_OppositeSign, h_jet_FiveJets_TwoBTagM_OppositeSign, h_event_FiveJets_TwoBTagM_OppositeSign;*/
  std::unique_ptr<Hists> lq_PreSelection, electron_PreSelection, muon_PreSelection, jet_PreSelection, tau_PreSelection, event_PreSelection;
  std::unique_ptr<Hists> h_lq_LeadingJet150, h_tau_LeadingJet150, h_mu_LeadingJet150, h_ele_LeadingJet150, h_jet_LeadingJet150, h_event_LeadingJet150;
  std::unique_ptr<Hists> h_lq_ThreeJets, h_tau_ThreeJets, h_mu_ThreeJets, h_ele_ThreeJets, h_jet_ThreeJets, h_event_ThreeJets;
  std::unique_ptr<Hists> h_lq_ST400, h_tau_ST400, h_mu_ST400, h_ele_ST400, h_jet_ST400, h_topjet_ST400, h_event_ST400, h_faketau_ST400;
  std::unique_ptr<Hists> h_lq_ST500, h_tau_ST500, h_mu_ST500, h_ele_ST500, h_jet_ST500, h_topjet_ST500, h_event_ST500, h_faketau_ST500;
  std::unique_ptr<Hists> h_lq_ST600, h_tau_ST600, h_mu_ST600, h_ele_ST600, h_jet_ST600, h_topjet_ST600, h_event_ST600, h_faketau_ST600;
  std::unique_ptr<Hists> h_lq_ST700, h_tau_ST700, h_mu_ST700, h_ele_ST700, h_jet_ST700, h_topjet_ST700, h_event_ST700, h_faketau_ST700;
  std::unique_ptr<Hists> h_lq_ST800, h_tau_ST800, h_mu_ST800, h_ele_ST800, h_jet_ST800, h_topjet_ST800, h_event_ST800, h_faketau_ST800;
  std::unique_ptr<Hists> h_lq_ST900, h_tau_ST900, h_mu_ST900, h_ele_ST900, h_jet_ST900, h_topjet_ST900, h_event_ST900, h_faketau_ST900;
  std::unique_ptr<Hists> h_lq_ST1000, h_tau_ST1000, h_mu_ST1000, h_ele_ST1000, h_jet_ST1000, h_topjet_ST1000, h_event_ST1000, h_faketau_ST1000;
  std::unique_ptr<Hists> h_lq_ST1100, h_tau_ST1100, h_mu_ST1100, h_ele_ST1100, h_jet_ST1100, h_topjet_ST1100, h_event_ST1100, h_faketau_ST1100;
  std::unique_ptr<Hists> h_lq_ST1200, h_tau_ST1200, h_mu_ST1200, h_ele_ST1200, h_jet_ST1200, h_topjet_ST1200, h_event_ST1200, h_faketau_ST1200;
  std::unique_ptr<Hists> h_lq_ST1300, h_tau_ST1300, h_mu_ST1300, h_ele_ST1300, h_jet_ST1300, h_topjet_ST1300, h_event_ST1300;
  std::unique_ptr<Hists> h_lq_ST1400, h_tau_ST1400, h_mu_ST1400, h_ele_ST1400, h_jet_ST1400, h_topjet_ST1400, h_event_ST1400;
  std::unique_ptr<Hists> h_lq_ST1500, h_tau_ST1500, h_mu_ST1500, h_ele_ST1500, h_jet_ST1500, h_topjet_ST1500, h_event_ST1500;
  std::unique_ptr<Hists> h_lq_ST1600, h_tau_ST1600, h_mu_ST1600, h_ele_ST1600, h_jet_ST1600, h_topjet_ST1600, h_event_ST1600;
  std::unique_ptr<Hists> h_lq_ST1700, h_tau_ST1700, h_mu_ST1700, h_ele_ST1700, h_jet_ST1700, h_topjet_ST1700, h_event_ST1700;
  std::unique_ptr<Hists> h_lq_ST1800, h_tau_ST1800, h_mu_ST1800, h_ele_ST1800, h_jet_ST1800, h_topjet_ST1800, h_event_ST1800;
  std::unique_ptr<Hists> h_lq_FactorTwo, h_tau_FactorTwo, h_mu_FactorTwo, h_ele_FactorTwo, h_jet_FactorTwo, h_event_FactorTwo;
  std::unique_ptr<Hists> h_ele_full, h_tau_full, h_event_full, h_lq_full, h_jet_full, h_muon_full;
  std::unique_ptr<Hists> ele_HT1000, muon_HT1000, tau_HT1000, event_HT1000, jet_HT1000, lq_HT1000;
  std::unique_ptr<Hists> ele_HT1100MET100_BTagM, muon_HT1100MET100_BTagM, tau_HT1100MET100_BTagM, event_HT1100MET100_BTagM, jet_HT1100MET100_BTagM, lq_HT1100MET100_BTagM;
  std::unique_ptr<Hists> ele_HT1100MET100, muon_HT1100MET100, tau_HT1100MET100, event_HT1100MET100, jet_HT1100MET100, lq_HT1100MET100;
  std::unique_ptr<Hists> ele_HT1000_BTagM, muon_HT1000_BTagM, tau_HT1000_BTagM, event_HT1000_BTagM, jet_HT1000_BTagM, lq_HT1000_BTagM;
  std::unique_ptr<Hists> ele_HT1000_BTagT, muon_HT1000_BTagT, tau_HT1000_BTagT, event_HT1000_BTagT, jet_HT1000_BTagT, lq_HT1000_BTagT;
  std::unique_ptr<Hists> ele_HT1000_TwoBTagL, muon_HT1000_TwoBTagL, tau_HT1000_TwoBTagL, event_HT1000_TwoBTagL, jet_HT1000_TwoBTagL, lq_HT1000_TwoBTagL;
  std::unique_ptr<Hists> ele_PreSel_SameSign, muon_PreSel_SameSign, tau_PreSel_SameSign, event_PreSel_SameSign, jet_PreSel_SameSign, lq_PreSel_SameSign;
  std::unique_ptr<Hists> ele_PreSel_OppositeSign, muon_PreSel_OppositeSign, tau_PreSel_OppositeSign, event_PreSel_OppositeSign, jet_PreSel_OppositeSign, lq_PreSel_OppositeSign;
  std::unique_ptr<Hists> ele_TwoMuon, muon_TwoMuon, tau_TwoMuon, event_TwoMuon, jet_TwoMuon, lq_TwoMuon;
  std::unique_ptr<Hists> ele_TwoMuon_ZCut, muon_TwoMuon_ZCut, tau_TwoMuon_ZCut, event_TwoMuon_ZCut, jet_TwoMuon_ZCut, lq_TwoMuon_ZCut;











  std::unique_ptr<AnalysisModule> muons_before, muons_after, jets_before, jets_after;


  JetId BTagLoose, BTagMedium, BTagTight;
  MuonId MuIso;
  ElectronId EleIso;

};


LQAnalysisTTbarSideBand2Module::LQAnalysisTTbarSideBand2Module(Context & ctx){
    // In the constructor, the typical tasks are to create
    // other modules like cleaners (1), selections (2) and Hist classes (3).
    // But you can do more and e.g. access the configuration, as shown below.
    
    cout << "Hello World from LQAnalysisTTbarSideBand2Module!" << endl;
    
    // If needed, access the configuration of the module here, e.g.:
    string testvalue = ctx.get("TestKey", "<not set>");
    cout << "TestKey in the configuration was: " << testvalue << endl;
    
    // If running in SFrame, the keys "dataset_version", "dataset_type" and "dataset_lumi"
    // are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }

    jets_before.reset(new JetPrinter("jets before", 30.0));
    jets_after.reset(new JetPrinter("jets after", 30.0));
    muons_before.reset(new MuonPrinter());
    muons_after.reset(new MuonPrinter());

    /*
    //test
    jet_lepton_cleaner.reset(new JetLeptonCleaner(JERFiles::PHYS14_L123_MC));
    jets_before.reset(new JetPrinter("jets before", 30.0));
    jets_after.reset(new JetPrinter("jets after", 30.0));
    muons_before.reset(new MuonPrinter());
    electrons_before.reset(new ElectronPrinter("before"));
    electrons_after.reset(new ElectronPrinter("after"));
    muon_cleaner.reset(new MuonCleaner(AndId<Muon>(PtEtaCut(30., 2.4), MuonIDTight())));
    ele_cleaner.reset(new ElectronCleaner(AndId<Electron>(PtEtaCut(30.0, 2.4), &ElectronID_CSA14_50ns_medium)));
    */

    // 1. setup other modules.
    MuIso = MuonIso(0.12);
    EleIso = ElectronIso(0.12);
    jetcleaner.reset(new JetCleaner(30.0, 2.5));
    muonidkinematic.reset(new MuonIDKinematic(30.0,3.0));
    muoncleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), MuonIDKinematic(30.0, 2.4))));
    muoncleaner_iso.reset(new MuonCleaner(AndId<Muon>(MuIso, MuonIDKinematic(30.0, 2.4))));
    electroncleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_PHYS14_25ns_medium, PtEtaCut(30.0, 2.5))));
    electroncleaner_iso.reset(new ElectronCleaner(AndId<Electron>(EleIso, PtEtaCut(30.0, 2.5))));
    taucleaner.reset(new TauCleaner(AndId<Tau>(TauIDMediumInverted(), PtEtaCut(30.0, 2.1))));
    BTagLoose = CSVBTag(CSVBTag::WP_LOOSE);
    BTagMedium = CSVBTag(CSVBTag::WP_MEDIUM);
    BTagTight = CSVBTag(CSVBTag::WP_TIGHT);
    jetleptoncleaner.reset(new JetLeptonCleaner(JERFiles::PHYS14_L123_MC));


    // 2. set up selections:
    njet_sel.reset(new NJetCut(2,-1,50,3.0));
    threejet_sel.reset(new NJetCut(4,-1,30,5.0));
    fourjet_sel.reset(new NJetCut(4,-1,30,5.0));
    fivejet_sel.reset(new NJetCut(5,-1,30,5.0));
    ntau_sel.reset(new NTauSelection(1,-1));
    nmuon_sel.reset(new NMuonSelection(1,-1));
    bsel.reset(new NBTagSelection(1,-1));
    TwoBTagL.reset(new NJetSelection(2,999,BTagLoose));
    TwoBTagM.reset(new NJetSelection(2,999,BTagMedium));
    BTagM.reset(new NJetSelection(1,999,BTagMedium));
    BTagT.reset(new NJetSelection(1,999,BTagTight));
    TwoBTagT.reset(new NJetSelection(2,999,BTagTight));
    samesign_sel.reset(new SameSignCut());
    samesign_lead.reset(new SameSignCutLeadingLep());


    int n_cuts = 4;
    fullhad_sel.resize(n_cuts);
    fullhad_sel[0].reset(new NJetSelection(2,-1));
    fullhad_sel[1].reset(new NTauSelection(1));
    //fullhad_sel[2].reset(new NJetSelection(1,999,BTagMedium));
    fullhad_sel[2].reset(new NJetCut(2,-1,50,3.0));
    fullhad_sel[3].reset(new NMuonSelection(1));
    //fullhad_sel[4].reset(new METCut(100,-1));



    // 3. Set up Hists classes:
    /*    h_lq_PreSel.reset(new LQAnalysisHists(ctx, "LQMod_LQ_PreSel"));
    h_tau_PreSel.reset(new TauHists(ctx, "LQMod_Taus_PreSel"));
    h_mu_PreSel.reset(new MuonHists(ctx, "LQMod_Muons_PreSel"));
    h_ele_PreSel.reset(new ElectronHists(ctx, "LQMod_Electrons_PreSel"));
    h_jet_PreSel.reset(new JetHists(ctx, "LQMod_Jets_PreSel"));
    h_event_PreSel.reset(new EventHists(ctx, "LQMod_Events_PreSel"));

    h_lq_FourJets.reset(new LQAnalysisHists(ctx, "LQMod_LQ_FourJets"));
    h_tau_FourJets.reset(new TauHists(ctx, "LQMod_Taus_FourJets"));
    h_mu_FourJets.reset(new MuonHists(ctx, "LQMod_Muons_FourJets"));
    h_ele_FourJets.reset(new ElectronHists(ctx, "LQMod_Electrons_FourJets"));
    h_jet_FourJets.reset(new JetHists(ctx, "LQMod_Jets_FourJets"));
    h_event_FourJets.reset(new EventHists(ctx, "LQMod_Events_FourJets"));

    h_lq_FourJets_SameSign.reset(new LQAnalysisHists(ctx, "LQMod_LQ_FourJets_SameSign"));
    h_tau_FourJets_SameSign.reset(new TauHists(ctx, "LQMod_Taus_FourJets_SameSign"));
    h_mu_FourJets_SameSign.reset(new MuonHists(ctx, "LQMod_Muons_FourJets_SameSign"));
    h_ele_FourJets_SameSign.reset(new ElectronHists(ctx, "LQMod_Electrons_FourJets_SameSign"));
    h_jet_FourJets_SameSign.reset(new JetHists(ctx, "LQMod_Jets_FourJets_SameSign"));
    h_event_FourJets_SameSign.reset(new EventHists(ctx, "LQMod_Events_FourJets_SameSign"));

    h_lq_FourJets_OppositeSign.reset(new LQAnalysisHists(ctx, "LQMod_LQ_FourJets_OppositeSign"));
    h_tau_FourJets_OppositeSign.reset(new TauHists(ctx, "LQMod_Taus_FourJets_OppositeSign"));
    h_mu_FourJets_OppositeSign.reset(new MuonHists(ctx, "LQMod_Muons_FourJets_OppositeSign"));
    h_ele_FourJets_OppositeSign.reset(new ElectronHists(ctx, "LQMod_Electrons_FourJets_OppositeSign"));
    h_jet_FourJets_OppositeSign.reset(new JetHists(ctx, "LQMod_Jets_FourJets_OppositeSign"));
    h_event_FourJets_OppositeSign.reset(new EventHists(ctx, "LQMod_Events_FourJets_OppositeSign"));

    h_lq_FiveJets.reset(new LQAnalysisHists(ctx, "LQMod_LQ_FiveJets"));
    h_tau_FiveJets.reset(new TauHists(ctx, "LQMod_Taus_FiveJets"));
    h_mu_FiveJets.reset(new MuonHists(ctx, "LQMod_Muons_FiveJets"));
    h_ele_FiveJets.reset(new ElectronHists(ctx, "LQMod_Electrons_FiveJets"));
    h_jet_FiveJets.reset(new JetHists(ctx, "LQMod_Jets_FiveJets"));
    h_event_FiveJets.reset(new EventHists(ctx, "LQMod_Events_FiveJets"));

    h_lq_FiveJets_TwoBTagM.reset(new LQAnalysisHists(ctx, "LQMod_LQ_FiveJets_TwoBTagM"));
    h_tau_FiveJets_TwoBTagM.reset(new TauHists(ctx, "LQMod_Taus_FiveJets_TwoBTagM"));
    h_mu_FiveJets_TwoBTagM.reset(new MuonHists(ctx, "LQMod_Muons_FiveJets_TwoBTagM"));
    h_ele_FiveJets_TwoBTagM.reset(new ElectronHists(ctx, "LQMod_Electrons_FiveJets_TwoBTagM"));
    h_jet_FiveJets_TwoBTagM.reset(new JetHists(ctx, "LQMod_Jets_FiveJets_TwoBTagM"));
    h_event_FiveJets_TwoBTagM.reset(new EventHists(ctx, "LQMod_Events_FiveJets_TwoBTagM"));
  
    h_lq_FiveJets_TwoBTagM_SameSign.reset(new LQAnalysisHists(ctx, "LQMod_LQ_FiveJets_TwoBTagM_SameSign"));
    h_tau_FiveJets_TwoBTagM_SameSign.reset(new TauHists(ctx, "LQMod_Taus_FiveJets_TwoBTagM_SameSign"));
    h_mu_FiveJets_TwoBTagM_SameSign.reset(new MuonHists(ctx, "LQMod_Muons_FiveJets_TwoBTagM_SameSign"));
    h_ele_FiveJets_TwoBTagM_SameSign.reset(new ElectronHists(ctx, "LQMod_Electrons_FiveJets_TwoBTagM_SameSign"));
    h_jet_FiveJets_TwoBTagM_SameSign.reset(new JetHists(ctx, "LQMod_Jets_FiveJets_TwoBTagM_SameSign"));
    h_event_FiveJets_TwoBTagM_SameSign.reset(new EventHists(ctx, "LQMod_Events_FiveJets_TwoBTagM_SameSign"));

    h_lq_FiveJets_TwoBTagM_OppositeSign.reset(new LQAnalysisHists(ctx, "LQMod_LQ_FiveJets_TwoBTagM_OppositeSign"));
    h_tau_FiveJets_TwoBTagM_OppositeSign.reset(new TauHists(ctx, "LQMod_Taus_FiveJets_TwoBTagM_OppositeSign"));
    h_mu_FiveJets_TwoBTagM_OppositeSign.reset(new MuonHists(ctx, "LQMod_Muons_FiveJets_TwoBTagM_OppositeSign"));
    h_ele_FiveJets_TwoBTagM_OppositeSign.reset(new ElectronHists(ctx, "LQMod_Electrons_FiveJets_TwoBTagM_OppositeSign"));
    h_jet_FiveJets_TwoBTagM_OppositeSign.reset(new JetHists(ctx, "LQMod_Jets_FiveJets_TwoBTagM_OppositeSign"));
    h_event_FiveJets_TwoBTagM_OppositeSign.reset(new EventHists(ctx, "LQMod_Events_FiveJets_TwoBTagM_OppositeSign"));  */
 
    electron_PreSelection.reset(new ElectronHists(ctx, "LQMod_Electrons_PreSel"));
    muon_PreSelection.reset(new MuonHists(ctx, "LQMod_Muons_PreSel"));
    tau_PreSelection.reset(new TauHists(ctx, "LQMod_Taus_PreSel"));
    jet_PreSelection.reset(new JetHists(ctx, "LQMod_Jets_PreSel"));
    event_PreSelection.reset(new EventHists(ctx, "LQMod_Events_PreSel"));
    lq_PreSelection.reset(new LQAnalysisHists(ctx, "LQMod_LQ_PreSel"));

    ele_TwoMuon.reset(new ElectronHists(ctx, "LQMod_Electrons_TwoMuon"));
    muon_TwoMuon.reset(new MuonHists(ctx, "LQMod_Muons_TwoMuon"));
    tau_TwoMuon.reset(new TauHists(ctx, "LQMod_Taus_TwoMuon"));
    jet_TwoMuon.reset(new JetHists(ctx, "LQMod_Jets_TwoMuon"));
    event_TwoMuon.reset(new EventHists(ctx, "LQMod_Events_TwoMuon"));
    lq_TwoMuon.reset(new LQAnalysisHists(ctx, "LQMod_LQ_TwoMuon"));

    ele_TwoMuon_ZCut.reset(new ElectronHists(ctx, "LQMod_Electrons_TwoMuon_ZCut"));
    muon_TwoMuon_ZCut.reset(new MuonHists(ctx, "LQMod_Muons_TwoMuon_ZCut"));
    tau_TwoMuon_ZCut.reset(new TauHists(ctx, "LQMod_Taus_TwoMuon_ZCut"));
    jet_TwoMuon_ZCut.reset(new JetHists(ctx, "LQMod_Jets_TwoMuon_ZCut"));
    event_TwoMuon_ZCut.reset(new EventHists(ctx, "LQMod_Events_TwoMuon_ZCut"));
    lq_TwoMuon_ZCut.reset(new LQAnalysisHists(ctx, "LQMod_LQ_TwoMuon_ZCut"));

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
    h_topjet_ST400.reset(new TopJetHists(ctx, "LQMod_TopJets_ST400"));
    h_event_ST400.reset(new EventHists(ctx, "LQMod_Events_ST400"));
    h_tau_ST400.reset(new TauHists(ctx, "LQMod_Taus_ST400"));
    h_lq_ST400.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST400"));
    h_faketau_ST400.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_ST400"));

    h_ele_ST500.reset(new ElectronHists(ctx, "LQMod_Electrons_ST500"));
    h_mu_ST500.reset(new MuonHists(ctx, "LQMod_Muons_ST500"));
    h_jet_ST500.reset(new JetHists(ctx, "LQMod_Jets_ST500"));
    h_topjet_ST500.reset(new TopJetHists(ctx, "LQMod_TopJets_ST500"));
    h_event_ST500.reset(new EventHists(ctx, "LQMod_Events_ST500"));
    h_tau_ST500.reset(new TauHists(ctx, "LQMod_Taus_ST500"));
    h_lq_ST500.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST500"));
    h_faketau_ST500.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_ST500"));

    h_ele_ST600.reset(new ElectronHists(ctx, "LQMod_Electrons_ST600"));
    h_mu_ST600.reset(new MuonHists(ctx, "LQMod_Muons_ST600"));
    h_jet_ST600.reset(new JetHists(ctx, "LQMod_Jets_ST600"));
    h_topjet_ST600.reset(new TopJetHists(ctx, "LQMod_TopJets_ST600"));
    h_event_ST600.reset(new EventHists(ctx, "LQMod_Events_ST600"));
    h_tau_ST600.reset(new TauHists(ctx, "LQMod_Taus_ST600"));
    h_lq_ST600.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST600"));
    h_faketau_ST600.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_ST600"));

    h_ele_ST700.reset(new ElectronHists(ctx, "LQMod_Electrons_ST700"));
    h_mu_ST700.reset(new MuonHists(ctx, "LQMod_Muons_ST700"));
    h_jet_ST700.reset(new JetHists(ctx, "LQMod_Jets_ST700"));
    h_topjet_ST700.reset(new TopJetHists(ctx, "LQMod_TopJets_ST700"));
    h_event_ST700.reset(new EventHists(ctx, "LQMod_Events_ST700"));
    h_tau_ST700.reset(new TauHists(ctx, "LQMod_Taus_ST700"));
    h_lq_ST700.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST700"));
    h_faketau_ST700.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_ST700"));

    h_ele_ST800.reset(new ElectronHists(ctx, "LQMod_Electrons_ST800"));
    h_mu_ST800.reset(new MuonHists(ctx, "LQMod_Muons_ST800"));
    h_jet_ST800.reset(new JetHists(ctx, "LQMod_Jets_ST800"));
    h_topjet_ST800.reset(new TopJetHists(ctx, "LQMod_TopJets_ST800"));
    h_event_ST800.reset(new EventHists(ctx, "LQMod_Events_ST800"));
    h_tau_ST800.reset(new TauHists(ctx, "LQMod_Taus_ST800"));
    h_lq_ST800.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST800"));
    h_faketau_ST800.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_ST800"));

    h_ele_ST900.reset(new ElectronHists(ctx, "LQMod_Electrons_ST900"));
    h_mu_ST900.reset(new MuonHists(ctx, "LQMod_Muons_ST900"));
    h_jet_ST900.reset(new JetHists(ctx, "LQMod_Jets_ST900"));
    h_topjet_ST900.reset(new TopJetHists(ctx, "LQMod_TopJets_ST900"));
    h_event_ST900.reset(new EventHists(ctx, "LQMod_Events_ST900"));
    h_tau_ST900.reset(new TauHists(ctx, "LQMod_Taus_ST900"));
    h_lq_ST900.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST900"));
    h_faketau_ST900.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_ST900"));

    h_lq_ST1000.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1000"));
    h_tau_ST1000.reset(new TauHists(ctx, "LQMod_Taus_ST1000"));
    h_mu_ST1000.reset(new MuonHists(ctx, "LQMod_Muons_ST1000"));
    h_ele_ST1000.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1000"));
    h_jet_ST1000.reset(new JetHists(ctx, "LQMod_Jets_ST1000"));
    h_topjet_ST1000.reset(new TopJetHists(ctx, "LQMod_TopJets_ST1000"));
    h_event_ST1000.reset(new EventHists(ctx, "LQMod_Events_ST1000"));
    h_faketau_ST1000.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_ST1000"));

    h_ele_ST1100.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1100"));
    h_mu_ST1100.reset(new MuonHists(ctx, "LQMod_Muons_ST1100"));
    h_jet_ST1100.reset(new JetHists(ctx, "LQMod_Jets_ST1100"));
    h_topjet_ST1100.reset(new TopJetHists(ctx, "LQMod_TopJets_ST1100"));
    h_event_ST1100.reset(new EventHists(ctx, "LQMod_Events_ST1100"));
    h_tau_ST1100.reset(new TauHists(ctx, "LQMod_Taus_ST1100"));
    h_lq_ST1100.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1100"));
    h_faketau_ST1100.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_ST1100"));

    h_ele_ST1200.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1200"));
    h_mu_ST1200.reset(new MuonHists(ctx, "LQMod_Muons_ST1200"));
    h_jet_ST1200.reset(new JetHists(ctx, "LQMod_Jets_ST1200"));
    h_topjet_ST1200.reset(new TopJetHists(ctx, "LQMod_TopJets_ST1200"));
    h_event_ST1200.reset(new EventHists(ctx, "LQMod_Events_ST1200"));
    h_tau_ST1200.reset(new TauHists(ctx, "LQMod_Taus_ST1200"));
    h_lq_ST1200.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1200"));
    h_faketau_ST1200.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_ST1200"));

    h_ele_ST1300.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1300"));
    h_mu_ST1300.reset(new MuonHists(ctx, "LQMod_Muons_ST1300"));
    h_jet_ST1300.reset(new JetHists(ctx, "LQMod_Jets_ST1300"));
    h_topjet_ST1300.reset(new TopJetHists(ctx, "LQMod_TopJets_ST1300"));
    h_event_ST1300.reset(new EventHists(ctx, "LQMod_Events_ST1300"));
    h_tau_ST1300.reset(new TauHists(ctx, "LQMod_Taus_ST1300"));
    h_lq_ST1300.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1300"));

    h_ele_ST1400.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1400"));
    h_mu_ST1400.reset(new MuonHists(ctx, "LQMod_Muons_ST1400"));
    h_jet_ST1400.reset(new JetHists(ctx, "LQMod_Jets_ST1400"));
    h_topjet_ST1400.reset(new TopJetHists(ctx, "LQMod_TopJets_ST1400"));
    h_event_ST1400.reset(new EventHists(ctx, "LQMod_Events_ST1400"));
    h_tau_ST1400.reset(new TauHists(ctx, "LQMod_Taus_ST1400"));
    h_lq_ST1400.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1400"));

    h_ele_ST1500.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1500"));
    h_mu_ST1500.reset(new MuonHists(ctx, "LQMod_Muons_ST1500"));
    h_jet_ST1500.reset(new JetHists(ctx, "LQMod_Jets_ST1500"));
    h_topjet_ST1500.reset(new TopJetHists(ctx, "LQMod_TopJets_ST1500"));
    h_event_ST1500.reset(new EventHists(ctx, "LQMod_Events_ST1500"));
    h_tau_ST1500.reset(new TauHists(ctx, "LQMod_Taus_ST1500"));
    h_lq_ST1500.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1500"));

    h_ele_ST1600.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1600"));
    h_mu_ST1600.reset(new MuonHists(ctx, "LQMod_Muons_ST1600"));
    h_jet_ST1600.reset(new JetHists(ctx, "LQMod_Jets_ST1600"));
    h_topjet_ST1600.reset(new TopJetHists(ctx, "LQMod_TopJets_ST1600"));
    h_event_ST1600.reset(new EventHists(ctx, "LQMod_Events_ST1600"));
    h_tau_ST1600.reset(new TauHists(ctx, "LQMod_Taus_ST1600"));
    h_lq_ST1600.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1600"));

    h_ele_ST1700.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1700"));
    h_mu_ST1700.reset(new MuonHists(ctx, "LQMod_Muons_ST1700"));
    h_jet_ST1700.reset(new JetHists(ctx, "LQMod_Jets_ST1700"));
    h_topjet_ST1700.reset(new TopJetHists(ctx, "LQMod_TopJets_ST1700"));
    h_event_ST1700.reset(new EventHists(ctx, "LQMod_Events_ST1700"));
    h_tau_ST1700.reset(new TauHists(ctx, "LQMod_Taus_ST1700"));
    h_lq_ST1700.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1700"));

    h_ele_ST1800.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1800"));
    h_mu_ST1800.reset(new MuonHists(ctx, "LQMod_Muons_ST1800"));
    h_jet_ST1800.reset(new JetHists(ctx, "LQMod_Jets_ST1800"));
    h_topjet_ST1800.reset(new TopJetHists(ctx, "LQMod_TopJets_ST1800"));
    h_event_ST1800.reset(new EventHists(ctx, "LQMod_Events_ST1800"));
    h_tau_ST1800.reset(new TauHists(ctx, "LQMod_Taus_ST1800"));
    h_lq_ST1800.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1800"));

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

    ele_PreSel_OppositeSign.reset(new ElectronHists(ctx, "ele_PreSel_OppositeSign"));
    muon_PreSel_OppositeSign.reset(new MuonHists(ctx, "muon_PreSel_OppositeSign"));
    jet_PreSel_OppositeSign.reset(new JetHists(ctx, "jet_PreSel_OppositeSign"));
    event_PreSel_OppositeSign.reset(new EventHists(ctx, "event_PreSel_OppositeSign"));
    tau_PreSel_OppositeSign.reset(new TauHists(ctx, "tau_PreSel_OppositeSign"));
    lq_PreSel_OppositeSign.reset(new LQAnalysisHists(ctx, "lq_PreSel_OppositeSign"));


}


bool LQAnalysisTTbarSideBand2Module::process(Event & event) {
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
    
  /*
  taucleaner->process(event);
  muoncleaner->process(event);
  electroncleaner->process(event);

  jetcleaner->process(event);
  */


  /*
  //print all trigger names    
  for (unsigned int i=0; i<event.get_current_triggernames().size();i++)
  cout<< event.get_current_triggernames()[i]<<"\n";
  cout<<"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
  */

    //jetlepcleantest->process(event);
  /*
    h_lq_PreSel->fill(event);
    h_tau_PreSel->fill(event);
    h_mu_PreSel->fill(event);
    h_ele_PreSel->fill(event);
    h_jet_PreSel->fill(event);
    h_event_PreSel->fill(event);
  */

  electron_PreSelection->fill(event);
  muon_PreSelection->fill(event);
  tau_PreSelection->fill(event);
  jet_PreSelection->fill(event);
  event_PreSelection->fill(event);
  lq_PreSelection->fill(event);





    
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
    

  const auto jets = event.jets;
  if(jets->size() > 0){
    const auto & jet = (*jets)[0];
    if(jet.pt()<150) return false;
  }


  for(const auto & tau : *event.taus){
    if(tau.pt()<30.0) return false;
  }



    //if(!fivejet_sel->passes(event)) return false;
    /*h_lq_FiveJets->fill(event);
    h_tau_FiveJets->fill(event);
    h_mu_FiveJets->fill(event);
    h_ele_FiveJets->fill(event);
    h_jet_FiveJets->fill(event);
    h_event_FiveJets->fill(event);*/

    /*
    for(const auto & tau : *event.taus){
      //if(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/tau.pt()>1.0 && tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/tau.pt()<1.1) return false;
      //if(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/tau.pt()>0.35 && tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()/tau.pt()<0.7) return false; // looks quite good if we use only 4 pt tau 1 bins...

      if(tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()>70 && tau.byCombinedIsolationDeltaBetaCorrRaw3Hits()<140) return false;

    }
    */



    const auto taus = event.taus;
    const auto & tau = (*event.taus)[0];

    //if(tau.pt()<100) return false;

    //if(!mbtau_sel->passes(event)) return false;
    
    bool OS_sel(false);
    if(OS_sel){
      if(samesign_lead->passes(event)) return false;
      if(!fourjet_sel->passes(event)) return false;
      if(met<100) return false;
    }
    else{
      if(!samesign_lead->passes(event)) return false;
      if(!BTagM->passes(event)) return false;
    }
    



  /*
    for(const auto & muon : *event.muons)
    {
    for(const auto & tau : *event.taus)
    {
    if (muon.charge() != tau.charge()){
    cout << "fail" << endl; 
    }
    else{
    cout << "success" << endl;
    }
    }     
    }
  */

  //if(!TwoBTagT->passes(event)) return false;
  /*h_lq_FiveJets_TwoBTagM->fill(event);
    h_tau_FiveJets_TwoBTagM->fill(event);
    h_mu_FiveJets_TwoBTagM->fill(event);
    h_ele_FiveJets_TwoBTagM->fill(event);
    h_jet_FiveJets_TwoBTagM->fill(event);
    h_event_FiveJets_TwoBTagM->fill(event);*/

    if(st<400) return false;
    h_ele_ST400->fill(event);
    h_mu_ST400->fill(event);
    h_jet_ST400->fill(event);
    h_topjet_ST400->fill(event);
    h_event_ST400->fill(event);
    h_tau_ST400->fill(event);
    h_lq_ST400->fill(event);
    h_faketau_ST400->fill(event);



    if(st<600) return false;
    h_ele_ST600->fill(event);
    h_mu_ST600->fill(event);
    h_jet_ST600->fill(event);
    h_topjet_ST600->fill(event);
    h_event_ST600->fill(event);
    h_tau_ST600->fill(event);
    h_lq_ST600->fill(event);
    h_faketau_ST600->fill(event);

    if(st<800) return false;
    h_ele_ST800->fill(event);
    h_mu_ST800->fill(event);
    h_jet_ST800->fill(event);
    h_topjet_ST800->fill(event);
    h_event_ST800->fill(event);
    h_tau_ST800->fill(event);
    h_lq_ST800->fill(event);
    h_faketau_ST800->fill(event);

    if(st<1000) return false;
    h_ele_ST1000->fill(event);
    h_mu_ST1000->fill(event);
    h_jet_ST1000->fill(event);
    h_topjet_ST1000->fill(event);
    h_event_ST1000->fill(event);
    h_tau_ST1000->fill(event);
    h_lq_ST1000->fill(event);
    h_faketau_ST1000->fill(event);

    if(st<1200) return false;
    h_ele_ST1200->fill(event);
    h_mu_ST1200->fill(event);
    h_jet_ST1200->fill(event);
    h_topjet_ST1200->fill(event);
    h_event_ST1200->fill(event);
    h_tau_ST1200->fill(event);
    h_lq_ST1200->fill(event);
    h_faketau_ST1200->fill(event);


    if(samesign_sel->passes(event)){
      /*h_lq_FiveJets_TwoBTagM_SameSign->fill(event);
      h_tau_FiveJets_TwoBTagM_SameSign->fill(event);
      h_mu_FiveJets_TwoBTagM_SameSign->fill(event);
      h_ele_FiveJets_TwoBTagM_SameSign->fill(event);
      h_jet_FiveJets_TwoBTagM_SameSign->fill(event);
      h_event_FiveJets_TwoBTagM_SameSign->fill(event);*/
    }

    if(!samesign_sel->passes(event)){
      /*h_lq_FiveJets_TwoBTagM_OppositeSign->fill(event);
      h_tau_FiveJets_TwoBTagM_OppositeSign->fill(event);
      h_mu_FiveJets_TwoBTagM_OppositeSign->fill(event);
      h_ele_FiveJets_TwoBTagM_OppositeSign->fill(event);
      h_jet_FiveJets_TwoBTagM_OppositeSign->fill(event);
      h_event_FiveJets_TwoBTagM_OppositeSign->fill(event);*/
    }

    
    // 3. decide whether or not to keep the current event in the output:
    return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisTTbarSideBand2Module)
