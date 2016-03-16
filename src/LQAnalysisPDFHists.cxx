#include "UHH2/LQAnalysis/include/LQAnalysisPDFHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include <math.h>
#include <sstream>

#include "TH1F.h"
#include "TH2D.h"
#include <iostream>

using namespace std;
using namespace uhh2;



LQAnalysisPDFHists::LQAnalysisPDFHists(Context & ctx, const string & dirname, bool use_pdf_weights_): Hists(ctx, dirname), use_pdf_weights(use_pdf_weights_){  
  
  is_mc = ctx.get("dataset_type") == "MC";
  double pttoprebins[5]={0,100,200,300,1200};
  double pttopbins[4]={0,100,200,1200};

  for(int i=0; i<100; i++){
    stringstream ss_name;
    ss_name << "N_PDF_"  << i+1 ;
    stringstream ss_name2;
    ss_name2 << "Pt_tophad_PDF_"  << i+1 ;
    stringstream ss_name3;
    ss_name3 << "Pt2_tophad_PDF_"  << i+1 ;
    stringstream ss_title;
    ss_title << "Counting exp. for PDF No. "  << i+1 << " out of 100" ;
    stringstream ss_title2;
    ss_title2 << "Pt tophad for PDF No. "  << i+1 << " out of 100" ;
    stringstream ss_title3;
    ss_title3 << "Pt2 tophad for PDF No. "  << i+1 << " out of 100" ;

    string s_name = ss_name.str();
    string s_title = ss_title.str();
    string s_name2 = ss_name2.str();
    string s_title2 = ss_title2.str();
    string s_name3 = ss_name3.str();
    string s_title3 = ss_title3.str();

    const char* char_name = s_name.c_str();
    const char* char_title = s_title.c_str();
    const char* char_name2 = s_name2.c_str();
    const char* char_title2 = s_title2.c_str();
    const char* char_name3 = s_name3.c_str();
    const char* char_title3 = s_title3.c_str();

    histo_names[i] = s_name;
    histo_names2[i] = s_name2;
    histo_names3[i] = s_name3;

    book<TH1F>(char_name, char_title, 1, 0.5, 1.5);
    book<TH1F>(char_name2, char_title2, 4, pttoprebins);
    book<TH1F>(char_name3, char_title3, 3, pttopbins);

  }

  h_hadr_hyps = ctx.get_handle<std::vector<TTbarFullhadRecoHypothesis>>("HighMassHadronicTTbarFullhadReco");

}

void LQAnalysisPDFHists::fill(const Event & event){
  double weight = event.weight;

  std::vector<TTbarFullhadRecoHypothesis> hadr_hyps = event.get(h_hadr_hyps);
  const TTbarFullhadRecoHypothesis* hadr_hyp = get_best_hypothesis( hadr_hyps, "Chi2Hadronic" );

  double ptTopHad = hadr_hyp->tophad1_v4().pt();

  if(is_mc){
    if(event.genInfo->systweights().size()){
      for(int i=0; i<100; i++){
	if(use_pdf_weights){
	  double pdf_weight = event.genInfo->systweights().at(i+9);
	  double fillweight = weight * pdf_weight/event.genInfo->originalXWGTUP();
	  const char* name = histo_names[i].c_str();
	  hist(name)->Fill(1,fillweight);
	}
      }
    }
    if(event.genInfo->systweights().size()){
      for(int i=0; i<100; i++){
	if(use_pdf_weights){
	  double pdf_weight = event.genInfo->systweights().at(i+9);
	  double fillweight = weight * pdf_weight/event.genInfo->originalXWGTUP();
	  const char* name2 = histo_names2[i].c_str();
	  const char* name3 = histo_names3[i].c_str();
	  hist(name2)->Fill(ptTopHad,fillweight);
	  hist(name3)->Fill(ptTopHad,fillweight);
	}
      }
    }

  } 

}

LQAnalysisPDFHists::~LQAnalysisPDFHists(){}














