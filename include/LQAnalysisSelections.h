#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"

namespace uhh2 {
    
/** NOTE: These classes are here as a (small) example only. If you need them, do NOT copy+paste these; the
 *  same (or even improved) functionality is also available through classes already in UHH2/common!
 */

/// Select events with certain minimum / maximum number of jets
/*
class NJetSelection: public Selection {
public:
    /// In case nmax=-1, no cut on the maximum is applied.
    explicit NJetSelection(int nmin, int nmax = -1);

    virtual bool passes(const Event & event) override;
    
private:
    int nmin, nmax;
};
*/

/** \brief Various definitions of b-tagging, in particular working points
 * 
 * This is useful for various selection modules, and thus defined outside of a particular Selection class.
 */
  
  class InvMass2MuVeto: public uhh2::Selection {
  public:
    explicit InvMass2MuVeto(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };


  class SameSignCut: public Selection{
  public:
    SameSignCut();
    ~SameSignCut(){};
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

  class MbtauSelection: public Selection{
  public:
    explicit MbtauSelection(double minMbtau=0., double maxMbtau=-1);
    virtual bool passes(const Event &);
  private:
    double minMbtau_, maxMbtau_;
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


/*
namespace lqanalysis_sel {
  using namespace uhh2;
  class NBTags: public Selection {
  public:
    /// In case nmax=-1, no cut on the maximum is applied.
    explicit NBTags(Context &ctx, int nmin_, int nmax_ = -1):
    hndl(ctx.get_handle<int>("n_btags")), nmin(nmin_), nmax(nmax_) {}
    virtual bool passes(const Event & event) override {
      int n = event.get(hndl);
      return n >= nmin && (nmax < 0 || n <= nmax);
    }
  private:
    Event::Handle<int> hndl;
    int nmin, nmax;
  }; 
}
*/


class ElectronIso {
public:
ElectronIso(double iso = 0.15);
bool operator()(const Electron & electron, const uhh2::Event & event) const;
private:
double iso;
};


