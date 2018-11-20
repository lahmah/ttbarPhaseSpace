// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_TTBAR : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_TTBAR);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      //declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");
      declare(FinalState(), "FS");

      // Book histograms
      _h_pT_t1 = bookHisto1D("pT_t1", 50, 0.0, 600.0);
      _h_pT_t2 = bookHisto1D("pT_t2", 50, 0.0, 600.0); 
      _h_E_t1 = bookHisto1D("E_t1", 50, 0.0, 600.0);
      _h_E_t2 = bookHisto1D("E_t2", 50, 0.0, 600.0);
      _h_eta_t1 = bookHisto1D("eta_t1", 50, -6.0, 6.0);
      _h_eta_t2 = bookHisto1D("eta_t2", 50, -6.0, 6.0);
      _h_angle_tt = bookHisto1D("angle_tt", 50, 0, M_PI);
      _h_t1_mass = bookHisto1D("t1_mass", 50, 0.0, 200.0); 
      _h_t2_mass = bookHisto1D("t2_mass", 50, 0.0, 200.0);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      FourMomentum mom_b;
      FourMomentum mom_W;

      const FinalState& fs = apply<FinalState>(event, "FS");
      for (const Particle& p : fs.particles()) {
        if (p.pid() == 6) {
          _h_pT_t1->fill(p.pT()/GeV, weight);
          _h_E_t1->fill(p.E()/GeV, weight);
          _h_eta_t1->fill(p.eta(), weight);
          _h_t1_mass->fill(p.mass()/GeV, weight);
	      }
        else if (p.pid() == -6) {
          _h_pT_t2->fill(p.pT()/GeV, weight);
          _h_E_t2->fill(p.E()/GeV, weight);
          _h_eta_t2->fill(p.eta(), weight);
          _h_t2_mass->fill(p.mass()/GeV, weight);
	      }
      }

      const Particles ps = fs.particles();
      // loop over all combinations
      for (unsigned int i=0; i<ps.size(); i++) {
        for (unsigned int j=0; j<ps.size(); j++) {
          if (ps[i].pid()==6 && ps[j].pid()==-6) {
            _h_angle_tt->fill(ps[i].angle(ps[j]), weight);
	        } 
	  
	      }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize
      normalize(_h_pT_t1);
      normalize(_h_E_t1);
      normalize(_h_eta_t1);
      normalize(_h_angle_tt);
      normalize(_h_t1_mass);
      normalize(_h_pT_t2);
      normalize(_h_E_t2);
      normalize(_h_eta_t2);
      normalize(_h_t2_mass);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_pT_t1, _h_E_t1, _h_eta_t1, _h_angle_tt, _h_t1_mass, _h_pT_t2, _h_E_t2, _h_eta_t2, _h_t2_mass; 
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TTBAR);


}
