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
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/LQAnalysis/include/LQAnalysisSelections.h"
#include "UHH2/LQAnalysis/include/LQAnalysisHists.h"
#include "UHH2/LQAnalysis/include/LQAnalysisPDFHists.h"
#include "UHH2/LQAnalysis/include/LQFakeTauHists.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/TauUncerts.h"
#include "UHH2/common/include/TopPtReweight.h"

#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQAnalysis/include/TTbarFullhadRecoHypothesisDiscriminators.h"
#include "UHH2/LQAnalysis/include/TTbarFullhadRecoHypothesis.h"
#include "UHH2/LQAnalysis/include/TTbarFullhadReco.h"

#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/ReconstructionHypothesis.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"
#include "UHH2/common/include/HypothesisHists.h"

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
  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  std::unique_ptr<Selection> njet_sel, twojet_sel, fourjet_sel, BTagM, ntau_sel, ExactlyOneTau, twotau_sel, muon_sel, oppositesign_sel, samesign_sel, samesign_lead, elesamesign_lead, MTcut, metcut, metcut_presel, getfaketaus, getrealtaus, ptleadingjetcut, possible_matches, OneIsoOneNonIso;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> lq_PreSelection, electron_PreSelection, muon_PreSelection, jet_PreSelection, tau_PreSelection, event_PreSelection, faketau_PreSelection;
  std::unique_ptr<Hists> lq_TwoTaus, electron_TwoTaus, muon_TwoTaus, jet_TwoTaus, tau_TwoTaus, event_TwoTaus, faketau_TwoTaus;
  std::unique_ptr<Hists> lq_TwoTaus_final, electron_TwoTaus_final, muon_TwoTaus_final, jet_TwoTaus_final, tau_TwoTaus_final, event_TwoTaus_final, faketau_TwoTaus_final;
  std::unique_ptr<Hists> h_lq_LeadingJet150, h_tau_LeadingJet150, h_mu_LeadingJet150, h_ele_LeadingJet150, h_jet_LeadingJet150, h_event_LeadingJet150, h_faketau_LeadingJet150, h_lumi_LeadingJet150;
  std::unique_ptr<Hists> h_lq_LeadingJet150_matchable, h_tau_LeadingJet150_matchable, h_mu_LeadingJet150_matchable, h_ele_LeadingJet150_matchable, h_jet_LeadingJet150_matchable, h_event_LeadingJet150_matchable, h_faketau_LeadingJet150_matchable;
  std::unique_ptr<Hists> h_lq_LeadingJet150_CorrectMatch, h_tau_LeadingJet150_CorrectMatch, h_mu_LeadingJet150_CorrectMatch, h_ele_LeadingJet150_CorrectMatch, h_jet_LeadingJet150_CorrectMatch, h_event_LeadingJet150_CorrectMatch, h_faketau_LeadingJet150_CorrectMatch;
  std::unique_ptr<Hists> h_lq_ThreeJets, h_tau_ThreeJets, h_mu_ThreeJets, h_ele_ThreeJets, h_jet_ThreeJets, h_event_ThreeJets, h_faketau_ThreeJets;
  std::unique_ptr<Hists> h_lq_OScut, h_tau_OScut, h_mu_OScut, h_ele_OScut, h_jet_OScut, h_event_OScut, h_faketau_OScut;
  std::unique_ptr<Hists> h_lq_OSfourjets, h_tau_OSfourjets, h_mu_OSfourjets, h_ele_OSfourjets, h_jet_OSfourjets, h_event_OSfourjets, h_faketau_OSfourjets;
  std::unique_ptr<Hists> h_lq_OSmet, h_tau_OSmet, h_mu_OSmet, h_ele_OSmet, h_jet_OSmet, h_event_OSmet, h_faketau_OSmet;
  std::unique_ptr<Hists> h_lq_OShtlep, h_tau_OShtlep, h_mu_OShtlep, h_ele_OShtlep, h_jet_OShtlep, h_event_OShtlep, h_faketau_OShtlep;
  std::unique_ptr<Hists> h_lq_OSpttau, h_tau_OSpttau, h_mu_OSpttau, h_ele_OSpttau, h_jet_OSpttau, h_event_OSpttau, h_faketau_OSpttau;
  std::unique_ptr<Hists> h_lq_SScut, h_tau_SScut, h_mu_SScut, h_ele_SScut, h_jet_SScut, h_event_SScut, h_faketau_SScut;
  std::unique_ptr<Hists> h_lq_SShtlep, h_tau_SShtlep, h_mu_SShtlep, h_ele_SShtlep, h_jet_SShtlep, h_event_SShtlep, h_faketau_SShtlep;
  std::unique_ptr<Hists> h_lq_SSpttau, h_tau_SSpttau, h_mu_SSpttau, h_ele_SSpttau, h_jet_SSpttau, h_event_SSpttau, h_faketau_SSpttau;
  std::unique_ptr<Hists> h_OS_lq_ST400, h_OS_tau_ST400, h_OS_mu_ST400, h_OS_ele_ST400, h_OS_jet_ST400, h_OS_topjet_ST400, h_OS_event_ST400, h_OS_faketau_ST400, h_OS_btageff_ST400;
  std::unique_ptr<Hists> h_OS_lq_lowmasses, h_OS_tau_lowmasses, h_OS_mu_lowmasses, h_OS_ele_lowmasses, h_OS_jet_lowmasses, h_OS_topjet_lowmasses, h_OS_event_lowmasses, h_OS_faketau_lowmasses, h_OS_pdf_lowmasses;
  std::unique_ptr<Hists> h_OS_lq_ST1200, h_OS_tau_ST1200, h_OS_mu_ST1200, h_OS_ele_ST1200, h_OS_jet_ST1200, h_OS_topjet_ST1200, h_OS_event_ST1200, h_OS_faketau_ST1200, h_OS_pdf_ST1200;
  std::unique_ptr<Hists> h_OS_lq_lowmasses_matchable, h_OS_tau_lowmasses_matchable, h_OS_mu_lowmasses_matchable, h_OS_ele_lowmasses_matchable, h_OS_jet_lowmasses_matchable, h_OS_topjet_lowmasses_matchable, h_OS_event_lowmasses_matchable, h_OS_faketau_lowmasses_matchable;
  std::unique_ptr<Hists> h_OS_lq_lowmasses_CorrectMatch, h_OS_tau_lowmasses_CorrectMatch, h_OS_mu_lowmasses_CorrectMatch, h_OS_ele_lowmasses_CorrectMatch, h_OS_jet_lowmasses_CorrectMatch, h_OS_topjet_lowmasses_CorrectMatch, h_OS_event_lowmasses_CorrectMatch, h_OS_faketau_lowmasses_CorrectMatch;
  std::unique_ptr<Hists> h_OS_lq_ST1200_matchable, h_OS_tau_ST1200_matchable, h_OS_mu_ST1200_matchable, h_OS_ele_ST1200_matchable, h_OS_jet_ST1200_matchable, h_OS_topjet_ST1200_matchable, h_OS_event_ST1200_matchable, h_OS_faketau_ST1200_matchable;
  std::unique_ptr<Hists> h_OS_lq_ST1200_CorrectMatch, h_OS_tau_ST1200_CorrectMatch, h_OS_mu_ST1200_CorrectMatch, h_OS_ele_ST1200_CorrectMatch, h_OS_jet_ST1200_CorrectMatch, h_OS_topjet_ST1200_CorrectMatch, h_OS_event_ST1200_CorrectMatch, h_OS_faketau_ST1200_CorrectMatch;
  std::unique_ptr<Hists> h_SS_lq_ST400, h_SS_tau_ST400, h_SS_mu_ST400, h_SS_ele_ST400, h_SS_jet_ST400, h_SS_topjet_ST400, h_SS_event_ST400, h_SS_faketau_ST400, h_SS_btageff_ST400;
  std::unique_ptr<Hists> h_SS_lq_lowmasses, h_SS_tau_lowmasses, h_SS_mu_lowmasses, h_SS_ele_lowmasses, h_SS_jet_lowmasses, h_SS_topjet_lowmasses, h_SS_event_lowmasses, h_SS_faketau_lowmasses, h_SS_pdf_lowmasses;
  std::unique_ptr<Hists> h_SS_lq_ST1200, h_SS_tau_ST1200, h_SS_mu_ST1200, h_SS_ele_ST1200, h_SS_jet_ST1200, h_SS_topjet_ST1200, h_SS_event_ST1200, h_SS_faketau_ST1200, h_SS_pdf_ST1200;
  std::unique_ptr<Hists> h_SS_lq_lowmasses_matchable, h_SS_tau_lowmasses_matchable, h_SS_mu_lowmasses_matchable, h_SS_ele_lowmasses_matchable, h_SS_jet_lowmasses_matchable, h_SS_topjet_lowmasses_matchable, h_SS_event_lowmasses_matchable, h_SS_faketau_lowmasses_matchable;
  std::unique_ptr<Hists> h_SS_lq_lowmasses_CorrectMatch, h_SS_tau_lowmasses_CorrectMatch, h_SS_mu_lowmasses_CorrectMatch, h_SS_ele_lowmasses_CorrectMatch, h_SS_jet_lowmasses_CorrectMatch, h_SS_topjet_lowmasses_CorrectMatch, h_SS_event_lowmasses_CorrectMatch, h_SS_faketau_lowmasses_CorrectMatch;
  std::unique_ptr<Hists> h_SS_lq_ST1200_matchable, h_SS_tau_ST1200_matchable, h_SS_mu_ST1200_matchable, h_SS_ele_ST1200_matchable, h_SS_jet_ST1200_matchable, h_SS_topjet_ST1200_matchable, h_SS_event_ST1200_matchable, h_SS_faketau_ST1200_matchable;
  std::unique_ptr<Hists> h_SS_lq_ST1200_CorrectMatch, h_SS_tau_ST1200_CorrectMatch, h_SS_mu_ST1200_CorrectMatch, h_SS_ele_ST1200_CorrectMatch, h_SS_jet_ST1200_CorrectMatch, h_SS_topjet_ST1200_CorrectMatch, h_SS_event_ST1200_CorrectMatch, h_SS_faketau_ST1200_CorrectMatch;

  std::unique_ptr<Hists> h_TopPtReweight;


  //std::unique_ptr<AnalysisModule> jetlepcleantest;
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<JetCleaner> jetcleaner;

  bool is_data, faketaus, realtaus, do_scale_variation, do_taueff_variation, do_taucharge_variation, do_tauenergy_variation, do_taures_variation, do_pdf_variations;
  MuonId MuId;
  ElectronId EleId;
  TauId TauonId;
  JetId BTagMedium;

  std::vector<std::unique_ptr<AnalysisModule>> pre_modules, recomodules;

  std::unique_ptr<AnalysisModule> ttgenprod, SF_muonID, SF_muonTrigger, SF_muonIso, SF_OS_btag, SF_SS_btag, syst_module, taueff_module, taucharge_module, tauenergy_module, taures_module;
  Event::Handle<TTbarGen> h_ttbargen;

  std::vector<std::unique_ptr<AnalysisModule>> recomodules_new;
  Event::Handle<TTbarGen> h_ttbargen_new;
  Event::Handle<std::vector<TTbarFullhadRecoHypothesis>> h_hadr_hyps_new;

  CSVBTag::wp wp_btag_medium;
  std::string channel_, Sys_PU, region_;
  std::unique_ptr<TopPtReweight> toppt_reweight;

};


LQAnalysisMuModule::LQAnalysisMuModule(Context & ctx){   
  // If running in SFrame, the keys "dataset_version", "dataset_type" and "dataset_lumi"
  // are set to the according values in the xml file. For CMSSW, these are
  // not set automatically, but can be set in the python config file.
  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }
    
  // 1. setup other modules.
  jetcleaner.reset(new JetCleaner(ctx, 30.0, 2.4));
  BTagMedium = CSVBTag(CSVBTag::WP_MEDIUM);
  wp_btag_medium = CSVBTag::WP_MEDIUM;
  toppt_reweight.reset(new TopPtReweight(ctx,0.0615,-0.0005,"ttbargen","weight_ttbar",false,1.02033407984698177));

  EleId = AndId<Electron>(ElectronID_Spring15_25ns_medium, PtEtaCut(30.0, 2.5));
  MuId = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.1),MuonIso(0.15));
  TauonId = AndId<Tau>(TauIDMedium(), PtEtaCut(20.0, 2.1));

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

  //systematics modules
  syst_module.reset(new MCScaleVariation(ctx));
  taueff_module.reset(new TauEffVariation(ctx));
  taucharge_module.reset(new TauChargeVariation(ctx));
  tauenergy_module.reset(new TauEnergySmearing(ctx));
  taures_module.reset(new TauEnergyResolutionShifter(ctx));

  channel_ = ctx.get("channel");
  region_ = ctx.get("region");
  is_data = ctx.get("dataset_type") == "DATA";
  faketaus = ctx.get("dataset_version") == "TTbar_onlyfakes";
  realtaus = ctx.get("dataset_version") == "TTbar_onlyrealtaus";
  auto muid_var=ctx.get("muonvariation");
  auto mutrig_var=ctx.get("muontriggervariation");
  auto Sys_MuonIso = ctx.get("Systematic_MuonIso");
  auto os_btag_var=ctx.get("OS_MCBtagEfficiencies");
  auto ss_btag_var=ctx.get("SS_MCBtagEfficiencies");

  do_scale_variation = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");
  do_taueff_variation = (ctx.get("TauIdVariation") == "nominal" || ctx.get("TauIdVariation") == "up" || ctx.get("TauIdVariation") == "down");
  do_taucharge_variation = (ctx.get("TauChargeVariation") == "up" || ctx.get("TauChargeVariation") == "down");
  do_tauenergy_variation = (ctx.get("TauEnergyVariation") == "up" || ctx.get("TauEnergyVariation") == "down");
  do_taures_variation = (ctx.get("TauEnergyResolutionVariation") == "up" || ctx.get("TauEnergyResolutionVariation") == "down");
  do_pdf_variations = ctx.get("b_PDFUncertainties") == "true";
  /*
  SF_muonID.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/mstoev/CMSSW_7_6_3/src/UHH2/common/data/MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1", 1, "tightID", muid_var));
  SF_muonTrigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/mstoev/CMSSW_7_6_3/src/UHH2/common/data/SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins", 0.5, "trigger", mutrig_var));
  SF_muonIso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/mstoev/CMSSW_7_6_3/src/UHH2/common/data/MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1", 1, "iso", Sys_MuonIso));
  
  */

  SF_OS_btag.reset(new MCBTagScaleFactor(ctx,wp_btag_medium,"jets",os_btag_var ,"mujets","incl","MCBtagEfficiencies_OS"));
  SF_SS_btag.reset(new MCBTagScaleFactor(ctx,wp_btag_medium,"jets",ss_btag_var ,"mujets","incl","MCBtagEfficiencies_SS"));


  // ttbar GEN
  ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
  // ttbar RECO hypotheses new
  h_hadr_hyps_new = ctx.get_handle<std::vector<TTbarFullhadRecoHypothesis>>("HighMassHadronicTTbarFullhadReco");
  recomodules_new.emplace_back(new HighMassHadronicTTbarReco(ctx));
  recomodules_new.emplace_back(new TTbarFullhadRecoChi2Discriminator(ctx,"HighMassHadronicTTbarFullhadReco"));
  //if(!is_data) recomodules_new.emplace_back(new LQCorrectMatchDiscriminator(ctx,"HighMassHadronicTTbarFullhadReco"));
  if(!is_data) recomodules_new.emplace_back(new LQCorrectMatchDiscriminator(ctx,"HighMassHadronicTTbarFullhadReco"));


  pre_modules.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));

  // 2. set up selections:
  twojet_sel.reset(new NJetSelection(2,-1));
  njet_sel.reset(new NJetSelection(3,-1));
  fourjet_sel.reset(new NJetSelection(4,-1));
  ntau_sel.reset(new NTauSelection(1,-1));
  ExactlyOneTau.reset(new NTauSelection(1,1));
  twotau_sel.reset(new NTauSelection(2,-1));
  BTagM.reset(new NJetSelection(1,999,BTagMedium));
  muon_sel.reset(new NMuonSelection(1,-1));
  metcut_presel.reset(new METCut(50,-1));
  metcut.reset(new METCut(100,-1));
  MTcut.reset(new MTCut(100,-1));
  ptleadingjetcut.reset(new PtLeadingJetSelection(150,-1));
  getfaketaus.reset(new GetFakeTaus());
  getrealtaus.reset(new GetRealTaus());
  OneIsoOneNonIso.reset(new IsoNonIso());
  oppositesign_sel.reset(new OppositeSignCut());
  samesign_sel.reset(new SameSignCut());
  if(channel_ == "muon"){
    samesign_lead.reset(new SameSignCutLeadingLep());
  }
  if(channel_ == "electron"){
    samesign_lead.reset(new EleTauSameSignCut());
  }
  elesamesign_lead.reset(new EleTauSameSignCut());
  possible_matches.reset(new LQMatch());


  // 3. Set up Hists classes:
  electron_PreSelection.reset(new ElectronHists(ctx, "LQMod_Electrons_PreSel"));
  muon_PreSelection.reset(new MuonHists(ctx, "LQMod_Muons_PreSel"));
  tau_PreSelection.reset(new TauHists(ctx, "LQMod_Taus_PreSel"));
  jet_PreSelection.reset(new JetHists(ctx, "LQMod_Jets_PreSel"));
  event_PreSelection.reset(new EventHists(ctx, "LQMod_Events_PreSel"));
  faketau_PreSelection.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_PreSel"));
  lq_PreSelection.reset(new LQAnalysisHists(ctx, "LQMod_LQ_PreSel"));

  electron_TwoTaus.reset(new ElectronHists(ctx, "LQMod_Electrons_TwoTaus"));
  muon_TwoTaus.reset(new MuonHists(ctx, "LQMod_Muons_TwoTaus"));
  tau_TwoTaus.reset(new TauHists(ctx, "LQMod_Taus_TwoTaus"));
  jet_TwoTaus.reset(new JetHists(ctx, "LQMod_Jets_TwoTaus"));
  event_TwoTaus.reset(new EventHists(ctx, "LQMod_Events_TwoTaus"));
  faketau_TwoTaus.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_TwoTaus"));
  lq_TwoTaus.reset(new LQAnalysisHists(ctx, "LQMod_LQ_TwoTaus"));
  
  electron_TwoTaus_final.reset(new ElectronHists(ctx, "LQMod_Electrons_TwoTaus_final"));
  muon_TwoTaus_final.reset(new MuonHists(ctx, "LQMod_Muons_TwoTaus_final"));
  tau_TwoTaus_final.reset(new TauHists(ctx, "LQMod_Taus_TwoTaus_final"));
  jet_TwoTaus_final.reset(new JetHists(ctx, "LQMod_Jets_TwoTaus_final"));
  event_TwoTaus_final.reset(new EventHists(ctx, "LQMod_Events_TwoTaus_final"));
  faketau_TwoTaus_final.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_TwoTaus_final"));
  lq_TwoTaus_final.reset(new LQAnalysisHists(ctx, "LQMod_LQ_TwoTaus_final"));

  h_ele_ThreeJets.reset(new ElectronHists(ctx, "LQMod_Electrons_ThreeJets"));
  h_mu_ThreeJets.reset(new MuonHists(ctx, "LQMod_Muons_ThreeJets"));
  h_jet_ThreeJets.reset(new JetHists(ctx, "LQMod_Jets_ThreeJets"));
  h_event_ThreeJets.reset(new EventHists(ctx, "LQMod_Events_ThreeJets"));
  h_tau_ThreeJets.reset(new TauHists(ctx, "LQMod_Taus_ThreeJets"));
  h_lq_ThreeJets.reset(new LQAnalysisHists(ctx, "LQMod_LQ_ThreeJets"));
  h_faketau_ThreeJets.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_ThreeJets"));

  h_lq_LeadingJet150.reset(new LQAnalysisHists(ctx, "LQMod_LQ_LeadingJet150"));
  h_tau_LeadingJet150.reset(new TauHists(ctx, "LQMod_Taus_LeadingJet150"));
  h_mu_LeadingJet150.reset(new MuonHists(ctx, "LQMod_Muons_LeadingJet150"));
  h_ele_LeadingJet150.reset(new ElectronHists(ctx, "LQMod_Electrons_LeadingJet150"));
  h_jet_LeadingJet150.reset(new JetHists(ctx, "LQMod_Jets_LeadingJet150"));
  h_event_LeadingJet150.reset(new EventHists(ctx, "LQMod_Events_LeadingJet150"));
  h_faketau_LeadingJet150.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_LeadingJet150"));
  h_lumi_LeadingJet150.reset(new LuminosityHists(ctx, "LQMod_Lumi_LeadingJet150"));

  h_lq_LeadingJet150_matchable.reset(new LQAnalysisHists(ctx, "LQMod_LQ_LeadingJet150_matchable"));
  h_tau_LeadingJet150_matchable.reset(new TauHists(ctx, "LQMod_Taus_LeadingJet150_matchable"));
  h_mu_LeadingJet150_matchable.reset(new MuonHists(ctx, "LQMod_Muons_LeadingJet150_matchable"));
  h_ele_LeadingJet150_matchable.reset(new ElectronHists(ctx, "LQMod_Electrons_LeadingJet150_matchable"));
  h_jet_LeadingJet150_matchable.reset(new JetHists(ctx, "LQMod_Jets_LeadingJet150_matchable"));
  h_event_LeadingJet150_matchable.reset(new EventHists(ctx, "LQMod_Events_LeadingJet150_matchable"));
  h_faketau_LeadingJet150_matchable.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_LeadingJet150_matchable"));

  h_lq_LeadingJet150_CorrectMatch.reset(new LQAnalysisHists(ctx, "LQMod_LQ_LeadingJet150_CorrectMatch"));
  h_tau_LeadingJet150_CorrectMatch.reset(new TauHists(ctx, "LQMod_Taus_LeadingJet150_CorrectMatch"));
  h_mu_LeadingJet150_CorrectMatch.reset(new MuonHists(ctx, "LQMod_Muons_LeadingJet150_CorrectMatch"));
  h_ele_LeadingJet150_CorrectMatch.reset(new ElectronHists(ctx, "LQMod_Electrons_LeadingJet150_CorrectMatch"));
  h_jet_LeadingJet150_CorrectMatch.reset(new JetHists(ctx, "LQMod_Jets_LeadingJet150_CorrectMatch"));
  h_event_LeadingJet150_CorrectMatch.reset(new EventHists(ctx, "LQMod_Events_LeadingJet150_CorrectMatch"));
  h_faketau_LeadingJet150_CorrectMatch.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_LeadingJet150_CorrectMatch"));

  h_lq_OScut.reset(new LQAnalysisHists(ctx, "LQMod_LQ_OScut"));
  h_tau_OScut.reset(new TauHists(ctx, "LQMod_Taus_OScut"));
  h_mu_OScut.reset(new MuonHists(ctx, "LQMod_Muons_OScut"));
  h_ele_OScut.reset(new ElectronHists(ctx, "LQMod_Electrons_OScut"));
  h_jet_OScut.reset(new JetHists(ctx, "LQMod_Jets_OScut"));
  h_event_OScut.reset(new EventHists(ctx, "LQMod_Events_OScut"));
  h_faketau_OScut.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_OScut"));

  h_lq_OSfourjets.reset(new LQAnalysisHists(ctx, "LQMod_LQ_OSfourjets"));
  h_tau_OSfourjets.reset(new TauHists(ctx, "LQMod_Taus_OSfourjets"));
  h_mu_OSfourjets.reset(new MuonHists(ctx, "LQMod_Muons_OSfourjets"));
  h_ele_OSfourjets.reset(new ElectronHists(ctx, "LQMod_Electrons_OSfourjets"));
  h_jet_OSfourjets.reset(new JetHists(ctx, "LQMod_Jets_OSfourjets"));
  h_event_OSfourjets.reset(new EventHists(ctx, "LQMod_Events_OSfourjets"));
  h_faketau_OSfourjets.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_OSfourjets"));

  h_lq_OSmet.reset(new LQAnalysisHists(ctx, "LQMod_LQ_OSmet"));
  h_tau_OSmet.reset(new TauHists(ctx, "LQMod_Taus_OSmet"));
  h_mu_OSmet.reset(new MuonHists(ctx, "LQMod_Muons_OSmet"));
  h_ele_OSmet.reset(new ElectronHists(ctx, "LQMod_Electrons_OSmet"));
  h_jet_OSmet.reset(new JetHists(ctx, "LQMod_Jets_OSmet"));
  h_event_OSmet.reset(new EventHists(ctx, "LQMod_Events_OSmet"));
  h_faketau_OSmet.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_OSmet"));

  h_lq_OShtlep.reset(new LQAnalysisHists(ctx, "LQMod_LQ_OShtlep"));
  h_tau_OShtlep.reset(new TauHists(ctx, "LQMod_Taus_OShtlep"));
  h_mu_OShtlep.reset(new MuonHists(ctx, "LQMod_Muons_OShtlep"));
  h_ele_OShtlep.reset(new ElectronHists(ctx, "LQMod_Electrons_OShtlep"));
  h_jet_OShtlep.reset(new JetHists(ctx, "LQMod_Jets_OShtlep"));
  h_event_OShtlep.reset(new EventHists(ctx, "LQMod_Events_OShtlep"));
  h_faketau_OShtlep.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_OShtlep"));

  h_lq_OSpttau.reset(new LQAnalysisHists(ctx, "LQMod_LQ_OSpttau"));
  h_tau_OSpttau.reset(new TauHists(ctx, "LQMod_Taus_OSpttau"));
  h_mu_OSpttau.reset(new MuonHists(ctx, "LQMod_Muons_OSpttau"));
  h_ele_OSpttau.reset(new ElectronHists(ctx, "LQMod_Electrons_OSpttau"));
  h_jet_OSpttau.reset(new JetHists(ctx, "LQMod_Jets_OSpttau"));
  h_event_OSpttau.reset(new EventHists(ctx, "LQMod_Events_OSpttau"));
  h_faketau_OSpttau.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_OSpttau"));

  h_lq_SScut.reset(new LQAnalysisHists(ctx, "LQMod_LQ_SScut"));
  h_tau_SScut.reset(new TauHists(ctx, "LQMod_Taus_SScut"));
  h_mu_SScut.reset(new MuonHists(ctx, "LQMod_Muons_SScut"));
  h_ele_SScut.reset(new ElectronHists(ctx, "LQMod_Electrons_SScut"));
  h_jet_SScut.reset(new JetHists(ctx, "LQMod_Jets_SScut"));
  h_event_SScut.reset(new EventHists(ctx, "LQMod_Events_SScut"));
  h_faketau_SScut.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_SScut"));

  h_lq_SShtlep.reset(new LQAnalysisHists(ctx, "LQMod_LQ_SShtlep"));
  h_tau_SShtlep.reset(new TauHists(ctx, "LQMod_Taus_SShtlep"));
  h_mu_SShtlep.reset(new MuonHists(ctx, "LQMod_Muons_SShtlep"));
  h_ele_SShtlep.reset(new ElectronHists(ctx, "LQMod_Electrons_SShtlep"));
  h_jet_SShtlep.reset(new JetHists(ctx, "LQMod_Jets_SShtlep"));
  h_event_SShtlep.reset(new EventHists(ctx, "LQMod_Events_SShtlep"));
  h_faketau_SShtlep.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_SShtlep"));

  h_lq_SSpttau.reset(new LQAnalysisHists(ctx, "LQMod_LQ_SSpttau"));
  h_tau_SSpttau.reset(new TauHists(ctx, "LQMod_Taus_SSpttau"));
  h_mu_SSpttau.reset(new MuonHists(ctx, "LQMod_Muons_SSpttau"));
  h_ele_SSpttau.reset(new ElectronHists(ctx, "LQMod_Electrons_SSpttau"));
  h_jet_SSpttau.reset(new JetHists(ctx, "LQMod_Jets_SSpttau"));
  h_event_SSpttau.reset(new EventHists(ctx, "LQMod_Events_SSpttau"));
  h_faketau_SSpttau.reset(new LQFakeTauHists(ctx, "LQMod_FakeTau_SSpttau"));

  h_OS_ele_ST400.reset(new ElectronHists(ctx, "LQMod_OS_Electrons_ST400"));
  h_OS_mu_ST400.reset(new MuonHists(ctx, "LQMod_OS_Muons_ST400"));
  h_OS_jet_ST400.reset(new JetHists(ctx, "LQMod_OS_Jets_ST400"));
  h_OS_topjet_ST400.reset(new TopJetHists(ctx, "LQMod_OS_TopJets_ST400"));
  h_OS_event_ST400.reset(new EventHists(ctx, "LQMod_OS_Events_ST400"));
  h_OS_tau_ST400.reset(new TauHists(ctx, "LQMod_OS_Taus_ST400"));
  h_OS_lq_ST400.reset(new LQAnalysisHists(ctx, "LQMod_OS_LQ_ST400"));
  h_OS_faketau_ST400.reset(new LQFakeTauHists(ctx, "LQMod_OS_FakeTau_ST400"));
  h_OS_btageff_ST400.reset(new BTagMCEfficiencyHists(ctx, "LQMod_OS_BTagEff_ST400",wp_btag_medium));

  h_OS_ele_lowmasses.reset(new ElectronHists(ctx, "LQMod_OS_Electrons_LowMasses"));
  h_OS_mu_lowmasses.reset(new MuonHists(ctx, "LQMod_OS_Muons_LowMasses"));
  h_OS_jet_lowmasses.reset(new JetHists(ctx, "LQMod_OS_Jets_LowMasses"));
  h_OS_topjet_lowmasses.reset(new TopJetHists(ctx, "LQMod_OS_TopJets_LowMasses"));
  h_OS_event_lowmasses.reset(new EventHists(ctx, "LQMod_OS_Events_LowMasses"));
  h_OS_tau_lowmasses.reset(new TauHists(ctx, "LQMod_OS_Taus_LowMasses"));
  h_OS_lq_lowmasses.reset(new LQAnalysisHists(ctx, "LQMod_OS_LQ_LowMasses"));
  h_OS_faketau_lowmasses.reset(new LQFakeTauHists(ctx, "LQMod_OS_FakeTau_LowMasses"));
  h_OS_pdf_lowmasses.reset(new LQAnalysisPDFHists(ctx, "LQMod_OS_PDF_LowMasses", do_pdf_variations));

  h_OS_ele_ST1200.reset(new ElectronHists(ctx, "LQMod_OS_Electrons_ST1200"));
  h_OS_mu_ST1200.reset(new MuonHists(ctx, "LQMod_OS_Muons_ST1200"));
  h_OS_jet_ST1200.reset(new JetHists(ctx, "LQMod_OS_Jets_ST1200"));
  h_OS_topjet_ST1200.reset(new TopJetHists(ctx, "LQMod_OS_TopJets_ST1200"));
  h_OS_event_ST1200.reset(new EventHists(ctx, "LQMod_OS_Events_ST1200"));
  h_OS_tau_ST1200.reset(new TauHists(ctx, "LQMod_OS_Taus_ST1200"));
  h_OS_lq_ST1200.reset(new LQAnalysisHists(ctx, "LQMod_OS_LQ_ST1200"));
  h_OS_faketau_ST1200.reset(new LQFakeTauHists(ctx, "LQMod_OS_FakeTau_ST1200"));
  h_OS_pdf_ST1200.reset(new LQAnalysisPDFHists(ctx, "LQMod_OS_PDF_ST1200", do_pdf_variations));

  h_OS_ele_lowmasses_matchable.reset(new ElectronHists(ctx, "LQMod_OS_Electrons_LowMasses_matchable"));
  h_OS_mu_lowmasses_matchable.reset(new MuonHists(ctx, "LQMod_OS_Muons_LowMasses_matchable"));
  h_OS_jet_lowmasses_matchable.reset(new JetHists(ctx, "LQMod_OS_Jets_LowMasses_matchable"));
  h_OS_topjet_lowmasses_matchable.reset(new TopJetHists(ctx, "LQMod_OS_TopJets_LowMasses_matchable"));
  h_OS_event_lowmasses_matchable.reset(new EventHists(ctx, "LQMod_OS_Events_LowMasses_matchable"));
  h_OS_tau_lowmasses_matchable.reset(new TauHists(ctx, "LQMod_OS_Taus_LowMasses_matchable"));
  h_OS_lq_lowmasses_matchable.reset(new LQAnalysisHists(ctx, "LQMod_OS_LQ_LowMasses_matchable"));
  h_OS_faketau_lowmasses_matchable.reset(new LQFakeTauHists(ctx, "LQMod_OS_FakeTau_LowMasses_matchable"));

  h_OS_ele_ST1200_matchable.reset(new ElectronHists(ctx, "LQMod_OS_Electrons_ST1200_matchable"));
  h_OS_mu_ST1200_matchable.reset(new MuonHists(ctx, "LQMod_OS_Muons_ST1200_matchable"));
  h_OS_jet_ST1200_matchable.reset(new JetHists(ctx, "LQMod_OS_Jets_ST1200_matchable"));
  h_OS_topjet_ST1200_matchable.reset(new TopJetHists(ctx, "LQMod_OS_TopJets_ST1200_matchable"));
  h_OS_event_ST1200_matchable.reset(new EventHists(ctx, "LQMod_OS_Events_ST1200_matchable"));
  h_OS_tau_ST1200_matchable.reset(new TauHists(ctx, "LQMod_OS_Taus_ST1200_matchable"));
  h_OS_lq_ST1200_matchable.reset(new LQAnalysisHists(ctx, "LQMod_OS_LQ_ST1200_matchable"));
  h_OS_faketau_ST1200_matchable.reset(new LQFakeTauHists(ctx, "LQMod_OS_FakeTau_ST1200_matchable"));

  h_OS_ele_lowmasses_CorrectMatch.reset(new ElectronHists(ctx, "LQMod_OS_Electrons_LowMasses_CorrectMatch"));
  h_OS_mu_lowmasses_CorrectMatch.reset(new MuonHists(ctx, "LQMod_OS_Muons_LowMasses_CorrectMatch"));
  h_OS_jet_lowmasses_CorrectMatch.reset(new JetHists(ctx, "LQMod_OS_Jets_LowMasses_CorrectMatch"));
  h_OS_topjet_lowmasses_CorrectMatch.reset(new TopJetHists(ctx, "LQMod_OS_TopJets_LowMasses_CorrectMatch"));
  h_OS_event_lowmasses_CorrectMatch.reset(new EventHists(ctx, "LQMod_OS_Events_LowMasses_CorrectMatch"));
  h_OS_tau_lowmasses_CorrectMatch.reset(new TauHists(ctx, "LQMod_OS_Taus_LowMasses_CorrectMatch"));
  h_OS_lq_lowmasses_CorrectMatch.reset(new LQAnalysisHists(ctx, "LQMod_OS_LQ_LowMasses_CorrectMatch"));
  h_OS_faketau_lowmasses_CorrectMatch.reset(new LQFakeTauHists(ctx, "LQMod_OS_FakeTau_LowMasses_CorrectMatch"));

  h_OS_ele_ST1200_CorrectMatch.reset(new ElectronHists(ctx, "LQMod_OS_Electrons_ST1200_CorrectMatch"));
  h_OS_mu_ST1200_CorrectMatch.reset(new MuonHists(ctx, "LQMod_OS_Muons_ST1200_CorrectMatch"));
  h_OS_jet_ST1200_CorrectMatch.reset(new JetHists(ctx, "LQMod_OS_Jets_ST1200_CorrectMatch"));
  h_OS_topjet_ST1200_CorrectMatch.reset(new TopJetHists(ctx, "LQMod_OS_TopJets_ST1200_CorrectMatch"));
  h_OS_event_ST1200_CorrectMatch.reset(new EventHists(ctx, "LQMod_OS_Events_ST1200_CorrectMatch"));
  h_OS_tau_ST1200_CorrectMatch.reset(new TauHists(ctx, "LQMod_OS_Taus_ST1200_CorrectMatch"));
  h_OS_lq_ST1200_CorrectMatch.reset(new LQAnalysisHists(ctx, "LQMod_OS_LQ_ST1200_CorrectMatch"));
  h_OS_faketau_ST1200_CorrectMatch.reset(new LQFakeTauHists(ctx, "LQMod_OS_FakeTau_ST1200_CorrectMatch"));

  h_SS_ele_ST400.reset(new ElectronHists(ctx, "LQMod_SS_Electrons_ST400"));
  h_SS_mu_ST400.reset(new MuonHists(ctx, "LQMod_SS_Muons_ST400"));
  h_SS_jet_ST400.reset(new JetHists(ctx, "LQMod_SS_Jets_ST400"));
  h_SS_topjet_ST400.reset(new TopJetHists(ctx, "LQMod_SS_TopJets_ST400"));
  h_SS_event_ST400.reset(new EventHists(ctx, "LQMod_SS_Events_ST400"));
  h_SS_tau_ST400.reset(new TauHists(ctx, "LQMod_SS_Taus_ST400"));
  h_SS_lq_ST400.reset(new LQAnalysisHists(ctx, "LQMod_SS_LQ_ST400"));
  h_SS_faketau_ST400.reset(new LQFakeTauHists(ctx, "LQMod_SS_FakeTau_ST400"));
  h_SS_btageff_ST400.reset(new BTagMCEfficiencyHists(ctx, "LQMod_SS_BTagEff_ST400",wp_btag_medium));

  h_SS_ele_lowmasses.reset(new ElectronHists(ctx, "LQMod_SS_Electrons_LowMasses"));
  h_SS_mu_lowmasses.reset(new MuonHists(ctx, "LQMod_SS_Muons_LowMasses"));
  h_SS_jet_lowmasses.reset(new JetHists(ctx, "LQMod_SS_Jets_LowMasses"));
  h_SS_topjet_lowmasses.reset(new TopJetHists(ctx, "LQMod_SS_TopJets_LowMasses"));
  h_SS_event_lowmasses.reset(new EventHists(ctx, "LQMod_SS_Events_LowMasses"));
  h_SS_tau_lowmasses.reset(new TauHists(ctx, "LQMod_SS_Taus_LowMasses"));
  h_SS_lq_lowmasses.reset(new LQAnalysisHists(ctx, "LQMod_SS_LQ_LowMasses"));
  h_SS_faketau_lowmasses.reset(new LQFakeTauHists(ctx, "LQMod_SS_FakeTau_LowMasses"));
  h_SS_pdf_lowmasses.reset(new LQAnalysisPDFHists(ctx, "LQMod_SS_PDF_LowMasses", do_pdf_variations));

  h_SS_ele_ST1200.reset(new ElectronHists(ctx, "LQMod_SS_Electrons_ST1200"));
  h_SS_mu_ST1200.reset(new MuonHists(ctx, "LQMod_SS_Muons_ST1200"));
  h_SS_jet_ST1200.reset(new JetHists(ctx, "LQMod_SS_Jets_ST1200"));
  h_SS_topjet_ST1200.reset(new TopJetHists(ctx, "LQMod_SS_TopJets_ST1200"));
  h_SS_event_ST1200.reset(new EventHists(ctx, "LQMod_SS_Events_ST1200"));
  h_SS_tau_ST1200.reset(new TauHists(ctx, "LQMod_SS_Taus_ST1200"));
  h_SS_lq_ST1200.reset(new LQAnalysisHists(ctx, "LQMod_SS_LQ_ST1200"));
  h_SS_faketau_ST1200.reset(new LQFakeTauHists(ctx, "LQMod_SS_FakeTau_ST1200"));
  h_SS_pdf_ST1200.reset(new LQAnalysisPDFHists(ctx, "LQMod_SS_PDF_ST1200", do_pdf_variations));
  
  h_SS_ele_lowmasses_matchable.reset(new ElectronHists(ctx, "LQMod_SS_Electrons_LowMasses_matchable"));
  h_SS_mu_lowmasses_matchable.reset(new MuonHists(ctx, "LQMod_SS_Muons_LowMasses_matchable"));
  h_SS_jet_lowmasses_matchable.reset(new JetHists(ctx, "LQMod_SS_Jets_LowMasses_matchable"));
  h_SS_topjet_lowmasses_matchable.reset(new TopJetHists(ctx, "LQMod_SS_TopJets_LowMasses_matchable"));
  h_SS_event_lowmasses_matchable.reset(new EventHists(ctx, "LQMod_SS_Events_LowMasses_matchable"));
  h_SS_tau_lowmasses_matchable.reset(new TauHists(ctx, "LQMod_SS_Taus_LowMasses_matchable"));
  h_SS_lq_lowmasses_matchable.reset(new LQAnalysisHists(ctx, "LQMod_SS_LQ_LowMasses_matchable"));
  h_SS_faketau_lowmasses_matchable.reset(new LQFakeTauHists(ctx, "LQMod_SS_FakeTau_LowMasses_matchable"));

  h_SS_ele_ST1200_matchable.reset(new ElectronHists(ctx, "LQMod_SS_Electrons_ST1200_matchable"));
  h_SS_mu_ST1200_matchable.reset(new MuonHists(ctx, "LQMod_SS_Muons_ST1200_matchable"));
  h_SS_jet_ST1200_matchable.reset(new JetHists(ctx, "LQMod_SS_Jets_ST1200_matchable"));
  h_SS_topjet_ST1200_matchable.reset(new TopJetHists(ctx, "LQMod_SS_TopJets_ST1200_matchable"));
  h_SS_event_ST1200_matchable.reset(new EventHists(ctx, "LQMod_SS_Events_ST1200_matchable"));
  h_SS_tau_ST1200_matchable.reset(new TauHists(ctx, "LQMod_SS_Taus_ST1200_matchable"));
  h_SS_lq_ST1200_matchable.reset(new LQAnalysisHists(ctx, "LQMod_SS_LQ_ST1200_matchable"));
  h_SS_faketau_ST1200_matchable.reset(new LQFakeTauHists(ctx, "LQMod_SS_FakeTau_ST1200_matchable"));

  h_SS_ele_lowmasses_CorrectMatch.reset(new ElectronHists(ctx, "LQMod_SS_Electrons_LowMasses_CorrectMatch"));
  h_SS_mu_lowmasses_CorrectMatch.reset(new MuonHists(ctx, "LQMod_SS_Muons_LowMasses_CorrectMatch"));
  h_SS_jet_lowmasses_CorrectMatch.reset(new JetHists(ctx, "LQMod_SS_Jets_LowMasses_CorrectMatch"));
  h_SS_topjet_lowmasses_CorrectMatch.reset(new TopJetHists(ctx, "LQMod_SS_TopJets_LowMasses_CorrectMatch"));
  h_SS_event_lowmasses_CorrectMatch.reset(new EventHists(ctx, "LQMod_SS_Events_LowMasses_CorrectMatch"));
  h_SS_tau_lowmasses_CorrectMatch.reset(new TauHists(ctx, "LQMod_SS_Taus_LowMasses_CorrectMatch"));
  h_SS_lq_lowmasses_CorrectMatch.reset(new LQAnalysisHists(ctx, "LQMod_SS_LQ_LowMasses_CorrectMatch"));
  h_SS_faketau_lowmasses_CorrectMatch.reset(new LQFakeTauHists(ctx, "LQMod_SS_FakeTau_LowMasses_CorrectMatch"));

  h_SS_ele_ST1200_CorrectMatch.reset(new ElectronHists(ctx, "LQMod_SS_Electrons_ST1200_CorrectMatch"));
  h_SS_mu_ST1200_CorrectMatch.reset(new MuonHists(ctx, "LQMod_SS_Muons_ST1200_CorrectMatch"));
  h_SS_jet_ST1200_CorrectMatch.reset(new JetHists(ctx, "LQMod_SS_Jets_ST1200_CorrectMatch"));
  h_SS_topjet_ST1200_CorrectMatch.reset(new TopJetHists(ctx, "LQMod_SS_TopJets_ST1200_CorrectMatch"));
  h_SS_event_ST1200_CorrectMatch.reset(new EventHists(ctx, "LQMod_SS_Events_ST1200_CorrectMatch"));
  h_SS_tau_ST1200_CorrectMatch.reset(new TauHists(ctx, "LQMod_SS_Taus_ST1200_CorrectMatch"));
  h_SS_lq_ST1200_CorrectMatch.reset(new LQAnalysisHists(ctx, "LQMod_SS_LQ_ST1200_CorrectMatch"));
  h_SS_faketau_ST1200_CorrectMatch.reset(new LQFakeTauHists(ctx, "LQMod_SS_FakeTau_ST1200_CorrectMatch"));

  h_TopPtReweight.reset(new TopPtReweightHist(ctx, "LQMod_TopPtReweighting", "weight_ttbar"));

}


bool LQAnalysisMuModule::process(Event & event) {
  //cout << "LQAnalysisMuModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;

  // process muon scale factors
  /*
  if(channel_ == "muon"){
    SF_muonTrigger->process(event);
    SF_muonID->process(event);
    SF_muonIso->process(event);
  }
  */

  bool pass_common = common->process(event);
  if(!pass_common) return false;
  jetcleaner->process(event);
  // HT calculator
  for (auto & mod : pre_modules) {
    mod->process(event);
  }

  if(!is_data){
    ttgenprod->process(event);
  }
  // top reconstruction
  for (auto & th : recomodules_new) {
    th->process(event);
  }
  
  toppt_reweight->process(event);
  h_TopPtReweight->fill(event);

  //if(!muon_sel->passes(event)) return false;
  if(!ntau_sel->passes(event)) return false;

  // process some systematic uncertainties
  if(do_scale_variation) syst_module->process(event);      
  if(do_taueff_variation) taueff_module->process(event);  
  if(do_taucharge_variation) taucharge_module->process(event);
  if(do_tauenergy_variation) tauenergy_module->process(event);
  if(do_taures_variation) taures_module->process(event);


  // cut events with real taus
  if(!is_data && faketaus){
    if(!getfaketaus->passes(event)) return false;
  }
  // cut events with fake taus
  if(!is_data && realtaus){
    if(!getrealtaus->passes(event)) return false;
  }
  

  if (!twojet_sel->passes(event)) return false;
  const auto jets = event.jets;
  if(jets->size() > 1){
    const auto & jet1 = (*jets)[0];
    const auto & jet2 = (*jets)[1];
    if(jet1.pt()<50) return false;
    if(jet2.pt()<50) return false;
  }

  if(!metcut_presel->passes(event)) return false;
  /*
  double chi2 = event.get(h_hadr_hyps_new)[0].chi2();
  if(chi2>4) return false;
  */


  //if(!possible_matches->passes(event)) return false;


  electron_PreSelection->fill(event);
  muon_PreSelection->fill(event);
  tau_PreSelection->fill(event);
  jet_PreSelection->fill(event);
  event_PreSelection->fill(event);
  faketau_PreSelection->fill(event);
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


  if (!njet_sel->passes(event)) return false;
  h_lq_ThreeJets->fill(event);
  h_tau_ThreeJets->fill(event);
  h_mu_ThreeJets->fill(event);
  h_ele_ThreeJets->fill(event);
  h_jet_ThreeJets->fill(event);
  h_event_ThreeJets->fill(event);
  h_faketau_ThreeJets->fill(event);


  //const auto taus = event.taus;  
  const auto & tau = (*event.taus)[0];

  if(twotau_sel->passes(event)){
    const auto & tau2 = (*event.taus)[1];
    if(region_!="CRtwo" && tau.charge()!=tau2.charge()){
      electron_TwoTaus->fill(event);
      muon_TwoTaus->fill(event);
      tau_TwoTaus->fill(event);
      jet_TwoTaus->fill(event);
      event_TwoTaus->fill(event);
      faketau_TwoTaus->fill(event);
      lq_TwoTaus->fill(event);
      if(tau.pt()>65 && tau2.pt()>35){
	electron_TwoTaus_final->fill(event);
	muon_TwoTaus_final->fill(event);
	tau_TwoTaus_final->fill(event);
	jet_TwoTaus_final->fill(event);
	event_TwoTaus_final->fill(event);
	faketau_TwoTaus_final->fill(event);
	lq_TwoTaus_final->fill(event);
      }
    }
    if(region_=="CRtwo" && tau.pt()<45 && !OneIsoOneNonIso->passes(event)){
      electron_TwoTaus->fill(event);
      muon_TwoTaus->fill(event);
      tau_TwoTaus->fill(event);
      jet_TwoTaus->fill(event);
      event_TwoTaus->fill(event);
      faketau_TwoTaus->fill(event);
      lq_TwoTaus->fill(event);
      if(BTagM->passes(event) && !samesign_lead->passes(event) && MTcut->passes(event)){
	electron_TwoTaus_final->fill(event);
	muon_TwoTaus_final->fill(event);
	tau_TwoTaus_final->fill(event);
	jet_TwoTaus_final->fill(event);
	event_TwoTaus_final->fill(event);
	faketau_TwoTaus_final->fill(event);
	lq_TwoTaus_final->fill(event);
      }
    }
  }

  //if(!ExactlyOneTau->passes(event)) return false;

  if(!ptleadingjetcut->passes(event)) return false;
  h_lq_LeadingJet150->fill(event);
  h_tau_LeadingJet150->fill(event);
  h_mu_LeadingJet150->fill(event);
  h_ele_LeadingJet150->fill(event);
  h_jet_LeadingJet150->fill(event);
  h_event_LeadingJet150->fill(event);
  h_faketau_LeadingJet150->fill(event);
  h_lumi_LeadingJet150->fill(event);

  /*
    for(const auto & jet : *event.jets){
    if(  fabs(event.get(h_hadr_hyps_new)[0].tophad1_v4().phi()-jet.phi()) >2.4  ){
    cout << "jet pt: " << jet.pt() << endl;
    cout << "mass:   " << (event.get(h_hadr_hyps_new)[0].tophad1_v4()+jet.v4()).M() << endl;
    }
    }
  */
  //const auto jets = event.jets;
  /*
    for(unsigned int i=0; i<jets->size(); i++){
    for(unsigned int j=0; j<jets->size(); j++){
    Jet jet1 = jets->at(i);
    Jet jet2 = jets->at(j);
    if(i!=j && fabs(jet1.phi()-jet2.phi()) > 2.5 ){
    cout << (jet1.v4()+jet2.v4()).M() << endl;
    }
    }}
    cout << endl;
  */


  //cout << "chi2 for t: " <<  event.get(h_hadr_hyps_new)[0].tchi2() << endl;
  //cout << "chi2 for W: " <<  event.get(h_hadr_hyps_new)[0].wchi2() << endl << endl;
  //cout << "wmass: " << event.get(h_hadr_hyps_new)[0].whad1_v4().M() << endl;

  if(!is_data){
    if(possible_matches->passes(event)){
      h_tau_LeadingJet150_matchable->fill(event);
      h_mu_LeadingJet150_matchable->fill(event);
      h_ele_LeadingJet150_matchable->fill(event);
      h_jet_LeadingJet150_matchable->fill(event);
      h_event_LeadingJet150_matchable->fill(event);
      h_faketau_LeadingJet150_matchable->fill(event);
      h_lq_LeadingJet150_matchable->fill(event);

      if (event.get(h_hadr_hyps_new)[0].discriminator("CorrectMatch") <99){
	//cout << "--> correctly matched +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl << endl;
	h_lq_LeadingJet150_CorrectMatch->fill(event);
	h_tau_LeadingJet150_CorrectMatch->fill(event);
	h_mu_LeadingJet150_CorrectMatch->fill(event);
	h_ele_LeadingJet150_CorrectMatch->fill(event);
	h_jet_LeadingJet150_CorrectMatch->fill(event);
	h_event_LeadingJet150_CorrectMatch->fill(event);
	h_faketau_LeadingJet150_CorrectMatch->fill(event);
      }
      //else{ cout << "not matched +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;}
    }
  }



  if(!samesign_lead->passes(event)){
    h_lq_OScut->fill(event);
    h_tau_OScut->fill(event);
    h_mu_OScut->fill(event);
    h_ele_OScut->fill(event);
    h_jet_OScut->fill(event);
    h_event_OScut->fill(event);
    h_faketau_OScut->fill(event);
    if(fourjet_sel->passes(event)){
      h_lq_OSfourjets->fill(event);
      h_tau_OSfourjets->fill(event);
      h_mu_OSfourjets->fill(event);
      h_ele_OSfourjets->fill(event);
      h_jet_OSfourjets->fill(event);
      h_event_OSfourjets->fill(event);
      h_faketau_OSfourjets->fill(event);
      if(metcut->passes(event)){
	h_lq_OSmet->fill(event);
	h_tau_OSmet->fill(event);
	h_mu_OSmet->fill(event);
	h_ele_OSmet->fill(event);
	h_jet_OSmet->fill(event);
	h_event_OSmet->fill(event);
	h_faketau_OSmet->fill(event);
	if(ht_lep>=0){ // was 200 before
	  h_lq_OShtlep->fill(event);
	  h_tau_OShtlep->fill(event);
	  h_mu_OShtlep->fill(event);
	  h_ele_OShtlep->fill(event);
	  h_jet_OShtlep->fill(event);
	  h_event_OShtlep->fill(event);
	  h_faketau_OShtlep->fill(event);
	  if(tau.pt()>=100){
	    h_lq_OSpttau->fill(event);
	    h_tau_OSpttau->fill(event);
	    h_mu_OSpttau->fill(event);
	    h_ele_OSpttau->fill(event);
	    h_jet_OSpttau->fill(event);
	    h_event_OSpttau->fill(event);
	    h_faketau_OSpttau->fill(event);
	    if(BTagM->passes(event)){
	      SF_OS_btag->process(event);
	      if(st>=400){
		h_OS_lq_ST400->fill(event);
		h_OS_tau_ST400->fill(event);
		h_OS_mu_ST400->fill(event);
		h_OS_ele_ST400->fill(event);
		h_OS_jet_ST400->fill(event);
		h_OS_topjet_ST400->fill(event);
		h_OS_event_ST400->fill(event);
		h_OS_faketau_ST400->fill(event);
		//h_OS_btageff_ST400->fill(event);
		if(st<1200){
		  h_OS_lq_lowmasses->fill(event);
		  h_OS_tau_lowmasses->fill(event);
		  h_OS_mu_lowmasses->fill(event);
		  h_OS_ele_lowmasses->fill(event);
		  h_OS_jet_lowmasses->fill(event);
		  h_OS_topjet_lowmasses->fill(event);
		  h_OS_event_lowmasses->fill(event);
		  h_OS_faketau_lowmasses->fill(event);
		  h_OS_pdf_lowmasses->fill(event);
		  if(!is_data){
		    if(possible_matches->passes(event)){
		      h_OS_lq_lowmasses_matchable->fill(event);
		      h_OS_tau_lowmasses_matchable->fill(event);
		      h_OS_mu_lowmasses_matchable->fill(event);
		      h_OS_ele_lowmasses_matchable->fill(event);
		      h_OS_jet_lowmasses_matchable->fill(event);
		      h_OS_topjet_lowmasses_matchable->fill(event);
		      h_OS_event_lowmasses_matchable->fill(event);
		      h_OS_faketau_lowmasses_matchable->fill(event);
		      if (event.get(h_hadr_hyps_new)[0].discriminator("CorrectMatch") <99){
			h_OS_lq_lowmasses_CorrectMatch->fill(event);
			h_OS_tau_lowmasses_CorrectMatch->fill(event);
			h_OS_mu_lowmasses_CorrectMatch->fill(event);
			h_OS_ele_lowmasses_CorrectMatch->fill(event);
			h_OS_jet_lowmasses_CorrectMatch->fill(event);
			h_OS_topjet_lowmasses_CorrectMatch->fill(event);
			h_OS_event_lowmasses_CorrectMatch->fill(event);
			h_OS_faketau_lowmasses_CorrectMatch->fill(event);
		      }
		    }
		  }
		}
		if(st>=1200){
		  h_OS_lq_ST1200->fill(event);
		  h_OS_tau_ST1200->fill(event);
		  h_OS_mu_ST1200->fill(event);
		  h_OS_ele_ST1200->fill(event);
		  h_OS_jet_ST1200->fill(event);
		  h_OS_topjet_ST1200->fill(event);
		  h_OS_event_ST1200->fill(event);
		  h_OS_faketau_ST1200->fill(event);
		  h_OS_pdf_ST1200->fill(event);
		  if(!is_data){
		    if(possible_matches->passes(event)){
		      h_OS_lq_ST1200_matchable->fill(event);
		      h_OS_tau_ST1200_matchable->fill(event);
		      h_OS_mu_ST1200_matchable->fill(event);
		      h_OS_ele_ST1200_matchable->fill(event);
		      h_OS_jet_ST1200_matchable->fill(event);
		      h_OS_topjet_ST1200_matchable->fill(event);
		      h_OS_event_ST1200_matchable->fill(event);
		      h_OS_faketau_ST1200_matchable->fill(event);
		      if (event.get(h_hadr_hyps_new)[0].discriminator("CorrectMatch") <99){
			h_OS_lq_ST1200_CorrectMatch->fill(event);
			h_OS_tau_ST1200_CorrectMatch->fill(event);
			h_OS_mu_ST1200_CorrectMatch->fill(event);
			h_OS_ele_ST1200_CorrectMatch->fill(event);
			h_OS_jet_ST1200_CorrectMatch->fill(event);
			h_OS_topjet_ST1200_CorrectMatch->fill(event);
			h_OS_event_ST1200_CorrectMatch->fill(event);
			h_OS_faketau_ST1200_CorrectMatch->fill(event);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }


  if(samesign_lead->passes(event)){
    h_lq_SScut->fill(event);
    h_tau_SScut->fill(event);
    h_mu_SScut->fill(event);
    h_ele_SScut->fill(event);
    h_jet_SScut->fill(event);
    h_event_SScut->fill(event);
    h_faketau_SScut->fill(event);
    //if(fourjet_sel->passes(event)){
    //if(metcut->passes(event)){
    if(ht_lep>=0){ //180
      h_lq_SShtlep->fill(event);
      h_tau_SShtlep->fill(event);
      h_mu_SShtlep->fill(event);
      h_ele_SShtlep->fill(event);
      h_jet_SShtlep->fill(event);
      h_event_SShtlep->fill(event);
      h_faketau_SShtlep->fill(event);
      if(tau.pt()>=100){ //60
	h_lq_SSpttau->fill(event);
	h_tau_SSpttau->fill(event);
	h_mu_SSpttau->fill(event);
	h_ele_SSpttau->fill(event);
	h_jet_SSpttau->fill(event);
	h_event_SSpttau->fill(event);
	h_faketau_SSpttau->fill(event);
	if(BTagM->passes(event)){
	  SF_SS_btag->process(event);
	  if(st>400){
	    h_SS_lq_ST400->fill(event);
	    h_SS_tau_ST400->fill(event);
	    h_SS_mu_ST400->fill(event);
	    h_SS_ele_ST400->fill(event);
	    h_SS_jet_ST400->fill(event);
	    h_SS_topjet_ST400->fill(event);
	    h_SS_event_ST400->fill(event);
	    h_SS_faketau_ST400->fill(event);
	    //h_SS_btageff_ST400->fill(event);
	    if(st<1200){
	      h_SS_lq_lowmasses->fill(event);
	      h_SS_tau_lowmasses->fill(event);
	      h_SS_mu_lowmasses->fill(event);
	      h_SS_ele_lowmasses->fill(event);
	      h_SS_jet_lowmasses->fill(event);
	      h_SS_topjet_lowmasses->fill(event);
	      h_SS_event_lowmasses->fill(event);
	      h_SS_faketau_lowmasses->fill(event);
	      h_SS_pdf_lowmasses->fill(event);
	      if(!is_data){
		if(possible_matches->passes(event)){
		  h_SS_lq_lowmasses_matchable->fill(event);
		  h_SS_tau_lowmasses_matchable->fill(event);
		  h_SS_mu_lowmasses_matchable->fill(event);
		  h_SS_ele_lowmasses_matchable->fill(event);
		  h_SS_jet_lowmasses_matchable->fill(event);
		  h_SS_topjet_lowmasses_matchable->fill(event);
		  h_SS_event_lowmasses_matchable->fill(event);
		  h_SS_faketau_lowmasses_matchable->fill(event);
		  if (event.get(h_hadr_hyps_new)[0].discriminator("CorrectMatch") <99){
		    h_SS_lq_lowmasses_CorrectMatch->fill(event);
		    h_SS_tau_lowmasses_CorrectMatch->fill(event);
		    h_SS_mu_lowmasses_CorrectMatch->fill(event);
		    h_SS_ele_lowmasses_CorrectMatch->fill(event);
		    h_SS_jet_lowmasses_CorrectMatch->fill(event);
		    h_SS_topjet_lowmasses_CorrectMatch->fill(event);
		    h_SS_event_lowmasses_CorrectMatch->fill(event);
		    h_SS_faketau_lowmasses_CorrectMatch->fill(event);
		  }
		}
	      }
	    }
	    if(st>=1200){
	      h_SS_lq_ST1200->fill(event);
	      h_SS_tau_ST1200->fill(event);
	      h_SS_mu_ST1200->fill(event);
	      h_SS_ele_ST1200->fill(event);
	      h_SS_jet_ST1200->fill(event);
	      h_SS_topjet_ST1200->fill(event);
	      h_SS_event_ST1200->fill(event);
	      h_SS_faketau_ST1200->fill(event);
	      h_SS_pdf_ST1200->fill(event);
	      if(!is_data){
		if(possible_matches->passes(event)){
		  h_SS_lq_ST1200_matchable->fill(event);
		  h_SS_tau_ST1200_matchable->fill(event);
		  h_SS_mu_ST1200_matchable->fill(event);
		  h_SS_ele_ST1200_matchable->fill(event);
		  h_SS_jet_ST1200_matchable->fill(event);
		  h_SS_topjet_ST1200_matchable->fill(event);
		  h_SS_event_ST1200_matchable->fill(event);
		  h_SS_faketau_ST1200_matchable->fill(event);
		  if (event.get(h_hadr_hyps_new)[0].discriminator("CorrectMatch") <99){
		    h_SS_lq_ST1200_CorrectMatch->fill(event);
		    h_SS_tau_ST1200_CorrectMatch->fill(event);
		    h_SS_mu_ST1200_CorrectMatch->fill(event);
		    h_SS_ele_ST1200_CorrectMatch->fill(event);
		    h_SS_jet_ST1200_CorrectMatch->fill(event);
		    h_SS_topjet_ST1200_CorrectMatch->fill(event);
		    h_SS_event_ST1200_CorrectMatch->fill(event);
		    h_SS_faketau_ST1200_CorrectMatch->fill(event);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  //}//
  //}//

  /*
  if(!is_data){
    for(const auto & muon : *event.muons){
      double dR = 1000;
      for(auto genp : *event.genparticles){
	if(abs(genp.pdgId())==13){
	  double tmp = deltaR(muon,genp);
	  if(tmp<dR){
	    dR = tmp;
	  }
	}
      }
      if(dR<0.1){
      }
      else{
	return false;
      }
    }
  }

  
  if(!is_data){
    for(const auto & tau : *event.taus){
      for(auto genp : *event.genparticles){
	//if(abs(genp.pdgId())!=15) continue;
	double dR = deltaR(tau,genp);
	if(abs(genp.pdgId())==15 && dR<0.4){
	  const GenParticle* mother1 = genp.mother(event.genparticles, 1);
	  //auto daughter1 = genp.daughter(event.genparticles, 1);
	  cout << "tau mother:  " << mother1->pdgId() << endl;
	  cout << "gen tau charge: " << genp.charge() << endl;
	  if(genp.charge()==tau.charge()){
	    cout << "true" << endl;
	  }
	  else{
	    cout << "false" << endl;
	  }
	  //cout << "tau daugher: " << daughter1->pdgId() << endl;
	}
      }
    }
  }

  if(!is_data){
    for(const auto & muon : *event.muons){
      for(auto genp : *event.genparticles){
	double dR = deltaR(muon,genp);
	if(abs(genp.pdgId())==13 && dR<0.1){
	  const GenParticle* mother1 = genp.mother(event.genparticles, 1);
	  cout << "muo mother:  " << mother1->pdgId() << endl;
	  cout << "gen muo charge: " << genp.charge() << endl;
	}
      }
    }
  }

  const auto & muon = (*event.muons)[0];
  //const auto & tau = (*event.taus)[0];
  cout << "tau charge: " << tau.charge() << endl;
  cout << "muo charge: " << muon.charge() << endl;
  */
    
  // 3. decide whether or not to keep the current event in the output:
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the LQAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(LQAnalysisMuModule)
