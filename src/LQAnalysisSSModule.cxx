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
class LQAnalysisSSModule: public AnalysisModule {
public:
    
    explicit LQAnalysisSSModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
    
    
  std::unique_ptr<JetCleaner> jetcleaner;
  std::unique_ptr<MuonIDKinematic> muonidkinematic;
  std::unique_ptr<MuonCleaner> muoncleaner;
  std::unique_ptr<TauCleaner> taucleaner;
  std::unique_ptr<ElectronCleaner> electroncleaner;
  
  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  std::unique_ptr<Selection> njet_sel, TwoBTagL, BTagM, BTagT, ntau_sel, ele_sel, muon_sel, onemuon_sel, twomuons_sel, isomuon_sel, leadingjet_sel, leadingjet300_sel, secondjet100_sel, samesign_sel, met_sel, mbtau_sel,Mmumu_cut;
  std::vector<std::unique_ptr<Selection> > fullhad_sel;

  std::vector<std::unique_ptr<AnalysisModule>> pre_modules;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> lq_PreSelection, electron_PreSelection, muon_PreSelection, jet_PreSelection, tau_PreSelection, event_PreSelection;
  std::unique_ptr<Hists> h_lq_LeadingJet150, h_tau_LeadingJet150, h_mu_LeadingJet150, h_ele_LeadingJet150, h_jet_LeadingJet150, h_event_LeadingJet150;
  std::unique_ptr<Hists> h_lq_ThreeJets, h_tau_ThreeJets, h_mu_ThreeJets, h_ele_ThreeJets, h_jet_ThreeJets, h_event_ThreeJets;
  std::unique_ptr<Hists> h_lq_ST410, h_tau_ST410, h_mu_ST410, h_ele_ST410, h_jet_ST410, h_event_ST410;
  std::unique_ptr<Hists> h_lq_ST470, h_tau_ST470, h_mu_ST470, h_ele_ST470, h_jet_ST470, h_event_ST470;
  std::unique_ptr<Hists> h_lq_ST470_TwoBTagL, h_tau_ST470_TwoBTagL, h_mu_ST470_TwoBTagL, h_ele_ST470_TwoBTagL, h_jet_ST470_TwoBTagL, h_event_ST470_TwoBTagL;
  std::unique_ptr<Hists> h_lq_ST500, h_tau_ST500, h_mu_ST500, h_ele_ST500, h_jet_ST500, h_event_ST500;
  std::unique_ptr<Hists> h_lq_ST600, h_tau_ST600, h_mu_ST600, h_ele_ST600, h_jet_ST600, h_event_ST600;
  std::unique_ptr<Hists> h_lq_ST770, h_tau_ST770, h_mu_ST770, h_ele_ST770, h_jet_ST770, h_event_ST770;
  std::unique_ptr<Hists> h_lq_ST800, h_tau_ST800, h_mu_ST800, h_ele_ST800, h_jet_ST800, h_event_ST800;
  std::unique_ptr<Hists> h_lq_ST900, h_tau_ST900, h_mu_ST900, h_ele_ST900, h_jet_ST900, h_event_ST900;
  std::unique_ptr<Hists> h_lq_ST1100, h_tau_ST1100, h_mu_ST1100, h_ele_ST1100, h_jet_ST1100, h_event_ST1100;
  std::unique_ptr<Hists> h_lq_ST1200, h_tau_ST1200, h_mu_ST1200, h_ele_ST1200, h_jet_ST1200, h_event_ST1200;
  std::unique_ptr<Hists> h_lq_ST1300, h_tau_ST1300, h_mu_ST1300, h_ele_ST1300, h_jet_ST1300, h_event_ST1300;
  std::unique_ptr<Hists> h_lq_ST1400, h_tau_ST1400, h_mu_ST1400, h_ele_ST1400, h_jet_ST1400, h_event_ST1400;
  std::unique_ptr<Hists> h_lq_ST1500, h_tau_ST1500, h_mu_ST1500, h_ele_ST1500, h_jet_ST1500, h_event_ST1500;

  JetId BTagLoose, BTagMedium, BTagTight;
  MuonId MuIso;


};


LQAnalysisSSModule::LQAnalysisSSModule(Context & ctx){
    // In the constructor, the typical tasks are to create
    // other modules like cleaners (1), selections (2) and Hist classes (3).
    // But you can do more and e.g. access the configuration, as shown below.
    
    cout << "Hello World from LQAnalysisSSModule!" << endl;
    
    // If needed, access the configuration of the module here, e.g.:
    string testvalue = ctx.get("TestKey", "<not set>");
    cout << "TestKey in the configuration was: " << testvalue << endl;
    
    // If running in SFrame, the keys "dataset_version", "dataset_type" and "dataset_lumi"
    // are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    /*for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
	}*/
    
    // 1. setup other modules.
    jetcleaner.reset(new JetCleaner(ctx, 30.0, 2.4));
    muonidkinematic.reset(new MuonIDKinematic(30.0,3.0));
    muoncleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(), MuonIDKinematic(30.0, 2.1))));
    electroncleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_PHYS14_25ns_medium, PtEtaCut(20.0, 2.5))));
    taucleaner.reset(new TauCleaner(AndId<Tau>(TauIDTight(), PtEtaCut(20.0, 2.1))));
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
    met_sel.reset(new METCut(150,-1));
    mbtau_sel.reset(new MbtauSelection(150,-1));
    Mmumu_cut.reset(new InvMass2MuVeto(81.,101.));


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

    h_ele_ST410.reset(new ElectronHists(ctx, "LQMod_Electrons_ST410"));
    h_mu_ST410.reset(new MuonHists(ctx, "LQMod_Muons_ST410"));
    h_jet_ST410.reset(new JetHists(ctx, "LQMod_Jets_ST410"));
    h_event_ST410.reset(new EventHists(ctx, "LQMod_Events_ST410"));
    h_tau_ST410.reset(new TauHists(ctx, "LQMod_Taus_ST410"));
    h_lq_ST410.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST410"));

    h_ele_ST470.reset(new ElectronHists(ctx, "LQMod_Electrons_ST470"));
    h_mu_ST470.reset(new MuonHists(ctx, "LQMod_Muons_ST470"));
    h_jet_ST470.reset(new JetHists(ctx, "LQMod_Jets_ST470"));
    h_event_ST470.reset(new EventHists(ctx, "LQMod_Events_ST470"));
    h_tau_ST470.reset(new TauHists(ctx, "LQMod_Taus_ST470"));
    h_lq_ST470.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST470"));

    h_ele_ST470_TwoBTagL.reset(new ElectronHists(ctx, "LQMod_Electrons_ST470_TwoBTagL"));
    h_mu_ST470_TwoBTagL.reset(new MuonHists(ctx, "LQMod_Muons_ST470_TwoBTagL"));
    h_jet_ST470_TwoBTagL.reset(new JetHists(ctx, "LQMod_Jets_ST470_TwoBTagL"));
    h_event_ST470_TwoBTagL.reset(new EventHists(ctx, "LQMod_Events_ST470_TwoBTagL"));
    h_tau_ST470_TwoBTagL.reset(new TauHists(ctx, "LQMod_Taus_ST470_TwoBTagL"));
    h_lq_ST470_TwoBTagL.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST470_TwoBTagL"));

    h_ele_ST500.reset(new ElectronHists(ctx, "LQMod_Electrons_ST500"));
    h_mu_ST500.reset(new MuonHists(ctx, "LQMod_Muons_ST500"));
    h_jet_ST500.reset(new JetHists(ctx, "LQMod_Jets_ST500"));
    h_event_ST500.reset(new EventHists(ctx, "LQMod_Events_ST500"));
    h_tau_ST500.reset(new TauHists(ctx, "LQMod_Taus_ST500"));
    h_lq_ST500.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST500"));

    h_ele_ST600.reset(new ElectronHists(ctx, "LQMod_Electrons_ST600"));
    h_mu_ST600.reset(new MuonHists(ctx, "LQMod_Muons_ST600"));
    h_jet_ST600.reset(new JetHists(ctx, "LQMod_Jets_ST600"));
    h_event_ST600.reset(new EventHists(ctx, "LQMod_Events_ST600"));
    h_tau_ST600.reset(new TauHists(ctx, "LQMod_Taus_ST600"));
    h_lq_ST600.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST600"));

    h_ele_ST770.reset(new ElectronHists(ctx, "LQMod_Electrons_ST770"));
    h_mu_ST770.reset(new MuonHists(ctx, "LQMod_Muons_ST770"));
    h_jet_ST770.reset(new JetHists(ctx, "LQMod_Jets_ST770"));
    h_event_ST770.reset(new EventHists(ctx, "LQMod_Events_ST770"));
    h_tau_ST770.reset(new TauHists(ctx, "LQMod_Taus_ST770"));
    h_lq_ST770.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST770"));

    h_ele_ST800.reset(new ElectronHists(ctx, "LQMod_Electrons_ST800"));
    h_mu_ST800.reset(new MuonHists(ctx, "LQMod_Muons_ST800"));
    h_jet_ST800.reset(new JetHists(ctx, "LQMod_Jets_ST800"));
    h_event_ST800.reset(new EventHists(ctx, "LQMod_Events_ST800"));
    h_tau_ST800.reset(new TauHists(ctx, "LQMod_Taus_ST800"));
    h_lq_ST800.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST800"));

    h_ele_ST900.reset(new ElectronHists(ctx, "LQMod_Electrons_ST900"));
    h_mu_ST900.reset(new MuonHists(ctx, "LQMod_Muons_ST900"));
    h_jet_ST900.reset(new JetHists(ctx, "LQMod_Jets_ST900"));
    h_event_ST900.reset(new EventHists(ctx, "LQMod_Events_ST900"));
    h_tau_ST900.reset(new TauHists(ctx, "LQMod_Taus_ST900"));
    h_lq_ST900.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST900"));

    h_ele_ST1100.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1100"));
    h_mu_ST1100.reset(new MuonHists(ctx, "LQMod_Muons_ST1100"));
    h_jet_ST1100.reset(new JetHists(ctx, "LQMod_Jets_ST1100"));
    h_event_ST1100.reset(new EventHists(ctx, "LQMod_Events_ST1100"));
    h_tau_ST1100.reset(new TauHists(ctx, "LQMod_Taus_ST1100"));
    h_lq_ST1100.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1100"));

    h_ele_ST1200.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1200"));
    h_mu_ST1200.reset(new MuonHists(ctx, "LQMod_Muons_ST1200"));
    h_jet_ST1200.reset(new JetHists(ctx, "LQMod_Jets_ST1200"));
    h_event_ST1200.reset(new EventHists(ctx, "LQMod_Events_ST1200"));
    h_tau_ST1200.reset(new TauHists(ctx, "LQMod_Taus_ST1200"));
    h_lq_ST1200.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1200"));

    h_ele_ST1300.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1300"));
    h_mu_ST1300.reset(new MuonHists(ctx, "LQMod_Muons_ST1300"));
    h_jet_ST1300.reset(new JetHists(ctx, "LQMod_Jets_ST1300"));
    h_event_ST1300.reset(new EventHists(ctx, "LQMod_Events_ST1300"));
    h_tau_ST1300.reset(new TauHists(ctx, "LQMod_Taus_ST1300"));
    h_lq_ST1300.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1300"));

    h_ele_ST1400.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1400"));
    h_mu_ST1400.reset(new MuonHists(ctx, "LQMod_Muons_ST1400"));
    h_jet_ST1400.reset(new JetHists(ctx, "LQMod_Jets_ST1400"));
    h_event_ST1400.reset(new EventHists(ctx, "LQMod_Events_ST1400"));
    h_tau_ST1400.reset(new TauHists(ctx, "LQMod_Taus_ST1400"));
    h_lq_ST1400.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1400"));

    h_ele_ST1500.reset(new ElectronHists(ctx, "LQMod_Electrons_ST1500"));
    h_mu_ST1500.reset(new MuonHists(ctx, "LQMod_Muons_ST1500"));
    h_jet_ST1500.reset(new JetHists(ctx, "LQMod_Jets_ST1500"));
    h_event_ST1500.reset(new EventHists(ctx, "LQMod_Events_ST1500"));
    h_tau_ST1500.reset(new TauHists(ctx, "LQMod_Taus_ST1500"));
    h_lq_ST1500.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ST1500"));

 
}


bool LQAnalysisSSModule::process(Event & event) {
    // This is the main procedure, called for each event. Typically,
    // do some pre-processing by calling the modules' process method
    // of the modules constructed in the constructor (1).
    // Then, test whether the event passes some selection and -- if yes --
    // use it to fill the histograms (2).
    // Finally, decide whether or not to keep the event in the output (3);
    // this is controlled by the return value of this method: If it
    // returns true, the event is kept; if it returns false, the event
    // is thrown away.
    
    //cout << "LQAnalysisSSModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
    
  //if(event.weight<0) return false;

    // 1. run all modules; here: only jet cleaning.
  /*
    taucleaner->process(event);
    muoncleaner->process(event);
    electroncleaner->process(event);
    jetcleaner->process(event);
  */
    for (auto & mod : pre_modules) {
      mod->process(event);
    }
    
    if(!muon_sel->passes(event)) return false;
    if(!ntau_sel->passes(event)) return false;
    //if(!met_sel->passes(event)) return false;
    //if(!mbtau_sel->passes(event)) return false;
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

    /*
    for(const auto & muon : *event.muons){
      if( (sqrt(2*muon.pt()*event.met->pt()* (1-cos(event.met->phi()-muon.phi())) ))<40) return false;
    }
    */

    if(!Mmumu_cut->passes(event)) return false;

    electron_PreSelection->fill(event);
    muon_PreSelection->fill(event);
    tau_PreSelection->fill(event);
    jet_PreSelection->fill(event);
    event_PreSelection->fill(event);
    lq_PreSelection->fill(event);




 

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


    for(const auto & tau : *event.taus){
      if(tau.pt()<35) return false;
    }
    if(!ntau_sel->passes(event)) return false;
    if(st<410) return false;
    h_lq_ST410->fill(event);
    h_tau_ST410->fill(event);
    h_mu_ST410->fill(event);
    h_ele_ST410->fill(event);
    h_jet_ST410->fill(event);
    h_event_ST410->fill(event);


    for(const auto & tau : *event.taus){
      if(tau.pt()<50) return false;
    }
    if(!ntau_sel->passes(event)) return false;
    if(st<470) return false;
    h_lq_ST470->fill(event);
    h_tau_ST470->fill(event);
    h_mu_ST470->fill(event);
    h_ele_ST470->fill(event);
    h_jet_ST470->fill(event);
    h_event_ST470->fill(event);

    if(TwoBTagL->passes(event)){
      h_lq_ST470_TwoBTagL->fill(event);
      h_tau_ST470_TwoBTagL->fill(event);
      h_mu_ST470_TwoBTagL->fill(event);
      h_ele_ST470_TwoBTagL->fill(event);
      h_jet_ST470_TwoBTagL->fill(event);
      h_event_ST470_TwoBTagL->fill(event);
    }



    //if(samesign_sel->passes(event)) return false;
    //if(met<100) return false;
    //if (!BTagM->passes(event)) return false;

    if(st<500) return false;
    h_lq_ST500->fill(event);
    h_tau_ST500->fill(event);
    h_mu_ST500->fill(event);
    h_ele_ST500->fill(event);
    h_jet_ST500->fill(event);
    h_event_ST500->fill(event);

    if(st<600) return false;
    h_lq_ST600->fill(event);
    h_tau_ST600->fill(event);
    h_mu_ST600->fill(event);
    h_ele_ST600->fill(event);
    h_jet_ST600->fill(event);
    h_event_ST600->fill(event);


    for(const auto & tau : *event.taus){
      if(tau.pt()<65) return false;
    }
    if(!ntau_sel->passes(event)) return false;
    if(st<770) return false;
    h_lq_ST770->fill(event);
    h_tau_ST770->fill(event);
    h_mu_ST770->fill(event);
    h_ele_ST770->fill(event);
    h_jet_ST770->fill(event);
    h_event_ST770->fill(event);

    
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

    if(st<800) return false;
    h_lq_ST800->fill(event);
    h_tau_ST800->fill(event);
    h_mu_ST800->fill(event);
    h_ele_ST800->fill(event);
    h_jet_ST800->fill(event);
    h_event_ST800->fill(event);

    if(st<900) return false;
    h_lq_ST900->fill(event);
    h_tau_ST900->fill(event);
    h_mu_ST900->fill(event);
    h_ele_ST900->fill(event);
    h_jet_ST900->fill(event);
    h_event_ST900->fill(event);

    if(st<1100) return false;
    h_lq_ST1100->fill(event);
    h_tau_ST1100->fill(event);
    h_mu_ST1100->fill(event);
    h_ele_ST1100->fill(event);
    h_jet_ST1100->fill(event);
    h_event_ST1100->fill(event);
    if(st>1200){
      h_lq_ST1200->fill(event);
      h_tau_ST1200->fill(event);
      h_mu_ST1200->fill(event);
      h_ele_ST1200->fill(event);
      h_jet_ST1200->fill(event);
      h_event_ST1200->fill(event);
    }
    if(st>1300){
      h_lq_ST1300->fill(event);
      h_tau_ST1300->fill(event);
      h_mu_ST1300->fill(event);
      h_ele_ST1300->fill(event);
      h_jet_ST1300->fill(event);
      h_event_ST1300->fill(event);
    }
    if(st>1400){
      h_lq_ST1400->fill(event);
      h_tau_ST1400->fill(event);
      h_mu_ST1400->fill(event);
      h_ele_ST1400->fill(event);
      h_jet_ST1400->fill(event);
      h_event_ST1400->fill(event);
    }
    if(st>1500){
      h_lq_ST1500->fill(event);
      h_tau_ST1500->fill(event);
      h_mu_ST1500->fill(event);
      h_ele_ST1500->fill(event);
      h_jet_ST1500->fill(event);
      h_event_ST1500->fill(event);
    }


    if(met<100) return false;
    //if(!leadingjet300_sel->passes(event)) return false;
    //if(!secondjet100_sel->passes(event)) return false;



    
    // 3. decide whether or not to keep the current event in the output:
    return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisSSModule)
