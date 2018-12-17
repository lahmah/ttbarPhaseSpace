// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_TTBAR_SINGLEDECAY : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_TTBAR_SINGLEDECAY);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      //declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");
      declare(FinalState(), "FS");

      // Book histograms
      _h_pT_t = bookHisto1D("pT_t", 50, 0.0, 600.0);
      _h_pT_b = bookHisto1D("pT_b", 50, 0.0, 600.0);
      _h_pT_e = bookHisto1D("pT_e", 50, 0.0, 600.0);
      _h_pT_nu = bookHisto1D("pT_nu", 50, 0.0, 600.0);
      _h_E_t = bookHisto1D("E_t", 50, 0.0, 600.0);
      _h_E_b = bookHisto1D("E_b", 50, 0.0, 600.0);
      _h_E_e = bookHisto1D("E_e", 50, 0.0, 600.0);
      _h_E_nu = bookHisto1D("E_nu", 50, 0.0, 600.0);
      _h_eta_t = bookHisto1D("eta_t", 50, -6.0, 6.0);
      _h_eta_b = bookHisto1D("eta_b", 50, -6.0, 6.0);
      _h_eta_e = bookHisto1D("eta_e", 50, -6.0, 6.0);
      _h_eta_nu = bookHisto1D("eta_nu", 50, -6.0, 6.0);
      _h_angle_tb = bookHisto1D("angle_tb", 50, 0, M_PI);
      _h_angle_te = bookHisto1D("angle_te", 50, 0, M_PI);
      _h_angle_tnu = bookHisto1D("angle_tnu", 50, 0, M_PI);
      _h_angle_be = bookHisto1D("angle_be", 50, 0, M_PI);
      _h_angle_bnu = bookHisto1D("angle_bnu", 50, 0, M_PI);
      _h_angle_enu = bookHisto1D("angle_enu", 50, 0, M_PI);
      _h_M_enu = bookHisto1D("M_enu", 50, 0.0, 1000.0);
      _h_M_trec = bookHisto1D("M_trec", 100,0,1000);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      FourMomentum mom_t(0., 0., 0., 0.);
      FourMomentum mom_b, mom_e, mom_nu;
      bool found_top = 0;

      const FinalState& fs = apply<FinalState>(event, "FS");
      for (const Particle& p : fs.particles()) {
        if (p.pid() == 6) {
          found_top = 1;
          _h_pT_t->fill(p.pT()/GeV, weight);
          _h_E_t->fill(p.E()/GeV, weight);
          _h_eta_t->fill(p.eta(), weight);
	}
	else if (p.pid() == -5) {
          mom_b = p.momentum();
          _h_pT_b->fill(p.pT()/GeV, weight);
          _h_E_b->fill(p.E()/GeV, weight);
          _h_eta_b->fill(p.eta(), weight);
	}
	else if (p.pid() == 11) {
          mom_e = p.momentum();
          _h_pT_e->fill(p.pT()/GeV, weight);
          _h_E_e->fill(p.E()/GeV, weight);
          _h_eta_e->fill(p.eta(), weight);
	}
	else if (p.pid() == -12) {
          mom_nu = p.momentum();
          _h_pT_nu->fill(p.pT()/GeV, weight);
          _h_E_nu->fill(p.E()/GeV, weight);
          _h_eta_nu->fill(p.eta(), weight);
	}
        // Hack for Sherpa
	else if (p.pid() == 5) {
	  mom_t += p.momentum();
	}
	else if (p.pid() == -11) {
	  mom_t += p.momentum();
	}
	else if (p.pid() == 12) {
	  mom_t += p.momentum();
	}
      }

      if (!found_top) {
        _h_pT_t->fill(mom_t.pT()/GeV, weight);
        _h_E_t->fill(mom_t.E()/GeV, weight);
        _h_eta_t->fill(mom_t.eta(), weight);
        _h_angle_tb->fill(mom_t.angle(mom_b), weight);
        _h_angle_te->fill(mom_t.angle(mom_e), weight);
        _h_angle_tnu->fill(mom_t.angle(mom_nu), weight);
       
      }
      
      const Particles ps = fs.particles();
      // loop over all combinations
      for (unsigned int i=0; i<ps.size(); i++) {
        for (unsigned int j=0; j<ps.size(); j++) {
          if (ps[i].pid()==6 && ps[j].pid()==-5) {
            _h_angle_tb->fill(ps[i].angle(ps[j]), weight);
	  }
	  else if (ps[i].pid()==6 && ps[j].pid()==11) {
            _h_angle_te->fill(ps[i].angle(ps[j]), weight);
	  }
	  else if (ps[i].pid()==6 && ps[j].pid()==-12) {
            _h_angle_tnu->fill(ps[i].angle(ps[j]), weight);
	  }
	  else if (ps[i].pid()==-5 && ps[j].pid()==11) {
            _h_angle_be->fill(ps[i].angle(ps[j]), weight);
	  }
	  else if (ps[i].pid()==-5 && ps[j].pid()==-12) {
            _h_angle_bnu->fill(ps[i].angle(ps[j]), weight);
	  }
	  else if (ps[i].pid()==11 && ps[j].pid()==-12) {
            _h_angle_enu->fill(ps[i].angle(ps[j]), weight);
	    FourMomentum sum_enu = ps[i].momentum() + ps[j].momentum();
            _h_M_enu->fill(sum_enu.mass()/GeV, weight);
            for (unsigned int k=0; k<ps.size(); k++){
              if ( ps[k].pid()==-5){
                FourMomentum sum_trec = sum_enu + ps[k].momentum();
                 _h_M_trec->fill(sum_trec.mass(),weight);

              }
            

            }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize
      normalize(_h_pT_t);
      normalize(_h_pT_b);
      normalize(_h_pT_e);
      normalize(_h_pT_nu);
      normalize(_h_E_t);
      normalize(_h_E_b);
      normalize(_h_E_e);
      normalize(_h_E_nu);
      normalize(_h_eta_t);
      normalize(_h_eta_b);
      normalize(_h_eta_e);
      normalize(_h_eta_nu);
      normalize(_h_angle_tb);
      normalize(_h_angle_te);
      normalize(_h_angle_tnu);
      normalize(_h_angle_be);
      normalize(_h_angle_bnu);
      normalize(_h_angle_enu);
      normalize(_h_M_enu);
      normalize(_h_M_trec);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_pT_t, _h_pT_b, _h_pT_e, _h_pT_nu, _h_E_t, _h_E_b, _h_E_e, _h_E_nu, _h_eta_t, _h_eta_b, _h_eta_e, _h_eta_nu, _h_angle_tb, _h_angle_te, _h_angle_tnu, _h_angle_be, _h_angle_bnu, _h_angle_enu, _h_M_enu, _h_M_trec;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TTBAR_SINGLEDECAY);


}
