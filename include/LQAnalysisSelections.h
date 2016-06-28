#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"

namespace uhh2 {
      
  class InvMass2MuVeto: public uhh2::Selection {
  public:
    explicit InvMass2MuVeto(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };

  class OppositeSignCut: public Selection{
  public:
    OppositeSignCut();
    ~OppositeSignCut(){};
    virtual bool passes(const Event & event);
  private:
  };

  class SameSignCut: public Selection{
  public:
    SameSignCut();
    ~SameSignCut(){};
    virtual bool passes(const Event & event);
  private:
  };

  class SameSignCutLeadingLep: public Selection{
  public:
    SameSignCutLeadingLep();
    ~SameSignCutLeadingLep(){};
    virtual bool passes(const Event & event);
  private:
  };

  class EleTauSameSignCut: public Selection{
  public:
    EleTauSameSignCut();
    ~EleTauSameSignCut(){};
    virtual bool passes(const Event & event);
  private:
  };

  class JetTauCleaning: public Selection{
  public:
    JetTauCleaning();
    ~JetTauCleaning(){};
    virtual bool passes(const Event & event);
  private:
  };

  class GetFakeTaus: public Selection{
  public:
    GetFakeTaus();
    ~GetFakeTaus(){};
    virtual bool passes(const Event & event);
  private:
  };

  class GetRealTaus: public Selection{
  public:
    GetRealTaus();
    ~GetRealTaus(){};
    virtual bool passes(const Event & event);
  private:
  };

  class MbtauSelection: public Selection{
  public:
    explicit MbtauSelection(double minMbtau=0., double maxMbtau=-1);
    virtual bool passes(const Event &);
  private:
    double minMbtau_, maxMbtau_;
  };

  class PtLeadingJetSelection: public uhh2::Selection {
  public:
    explicit PtLeadingJetSelection(double pt_min = 30., double m_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double pt_min, pt_max;
  };

  class NJetCut: public Selection {
  public:
    /// In case nmax=-1, no cut on the maximum is applied.
    explicit NJetCut(int nmin, int nmax = -1, double ptmin=0., double etamax = -1);
    virtual bool passes(const Event & event) override;
  private:
    int nmin, nmax;
    double ptmin, etamax;
  };


  class METCut: public Selection {
  public:
    explicit METCut(double min_met=0., double max_met=-1);
    virtual bool passes(const Event &);
  private:
    double min_met, max_met;
  };

  class HtSelection: public uhh2::Selection {
  public:
    explicit HtSelection(double ht_min=0., double ht_max=-1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ht_min, ht_max;
  };



  namespace btagging {
    
    enum class csv_wp {
      loose, medium, tight
	};

    /// convert a CSV working point to a numerical threshold of the discriminator.
    float csv_threshold(const csv_wp & wp);

  }

  /// Select events with certain minimum / maximum number of b-tagged jets using the CSV tagger
  class NBTagSelection: public Selection {
  public:
    /// In case nmax=-1, no cut on the maximum is applied.
    explicit NBTagSelection(int nmin, int nmax = -1, btagging::csv_wp wp = btagging::csv_wp::medium);
    virtual bool passes(const Event & event) override;
    
  private:
    int nmin, nmax;
    float min_csv;
  };

 
  class TestCut{
  public:
  TestCut(double min_met_, double max_met_): min_met(min_met_), max_met(max_met_){}
    bool operator()(const uhh2::Event & event) const{
      return event.met->pt() > min_met && event.met->pt() < max_met;
    }
  private:
    float min_met, max_met;
  };


}

class ElectronIso {
 public:
  ElectronIso(double iso = 0.15);
  bool operator()(const Electron & electron, const uhh2::Event & event) const;
 private:
  double iso;
};




