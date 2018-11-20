// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_TTBAR_DOUBLEDECAY : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_TTBAR_DOUBLEDECAY);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      //declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");
      declare(FinalState(), "FS");

      // Book histograms
      _h_pT_b = bookHisto1D("pT_b", 50, 0.0, 600.0);
      _h_pT_bb = bookHisto1D("pT_bb", 50, 0.0, 600.0);
      _h_pT_e = bookHisto1D("pT_e", 50, 0.0, 600.0);
      _h_pT_nue = bookHisto1D("pT_nue", 50, 0.0, 600.0);
      _h_pT_mu = bookHisto1D("pT_mu", 50, 0.0, 600.0);
      _h_pT_numu = bookHisto1D("pT_numu", 50, 0.0, 600.0);
      _h_E_b = bookHisto1D("E_b", 50, 0.0, 600.0);
      _h_E_bb = bookHisto1D("E_bb", 50, 0.0, 600.0);
      _h_E_e = bookHisto1D("E_e", 50, 0.0, 600.0);
      _h_E_nue = bookHisto1D("E_nue", 50, 0.0, 600.0);
      _h_E_mu = bookHisto1D("E_mu", 50, 0.0, 600.0);
      _h_E_numu = bookHisto1D("E_numu", 50, 0.0, 600.0);
      _h_eta_b = bookHisto1D("eta_b", 50, -6.0, 6.0);
      _h_eta_bb = bookHisto1D("eta_bb", 50, -6.0, 6.0);
      _h_eta_e = bookHisto1D("eta_e", 50, -6.0, 6.0);
      _h_eta_nue = bookHisto1D("eta_nue", 50, -6.0, 6.0);
      _h_eta_mu = bookHisto1D("eta_mu", 50, -6.0, 6.0);
      _h_eta_numu = bookHisto1D("eta_numu", 50, -6.0, 6.0);
      _h_angle_bb = bookHisto1D("angle_bb", 50, 0, M_PI);
      _h_angle_munu = bookHisto1D("angle_munu", 50, 0, M_PI);
      _h_angle_enu = bookHisto1D("angle_enu", 50, 0, M_PI);
      _h_angle_ebb = bookHisto1D("angle_ebb", 50, 0, M_PI);
      _h_angle_nubb = bookHisto1D("angle_nubb", 50, 0, M_PI);
      _h_angle_Wmbb = bookHisto1D("angle_Wmbb", 50, 0, M_PI);
      _h_angle_mub = bookHisto1D("angle_mub", 50, 0, M_PI);
      _h_angle_nub = bookHisto1D("angle_nub", 50, 0, M_PI);
      _h_angle_Wpb = bookHisto1D("angle_Wpb", 50, 0, M_PI);
      _h_Wp_mass = bookHisto1D("Wp_mass", 50, 0.0, 1000.0);
      _h_Wm_mass = bookHisto1D("Wm_mass", 50, 0.0, 1000.0);
      _h_t_mass = bookHisto1D("t_mass", 50, 0.0, 1000.0);
      _h_tb_mass = bookHisto1D("tb_mass", 50, 0.0, 1000.0);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const FinalState& fs = apply<FinalState>(event, "FS");

      FourMomentum mom_b;
      FourMomentum mom_bb;
      FourMomentum mom_Wm;
      FourMomentum mom_Wp;

      for (const Particle& p : fs.particles()) {
        if (p.pid() == 5) {
          _h_pT_b->fill(p.pT()/GeV, weight);
          _h_E_b->fill(p.E()/GeV, weight);
          _h_eta_b->fill(p.eta(), weight);
	  mom_b = p.momentum();
	}
	else if (p.pid() == -5) {
          _h_pT_bb->fill(p.pT()/GeV, weight);
          _h_E_bb->fill(p.E()/GeV, weight);
          _h_eta_bb->fill(p.eta(), weight);
	  mom_bb = p.momentum();
	}
	else if (p.pid() == 11) {
          _h_pT_e->fill(p.pT()/GeV, weight);
          _h_E_e->fill(p.E()/GeV, weight);
          _h_eta_e->fill(p.eta(), weight);
	}
	else if (p.pid() == -12) {
          _h_pT_nue->fill(p.pT()/GeV, weight);
          _h_E_nue->fill(p.E()/GeV, weight);
          _h_eta_nue->fill(p.eta(), weight);
	}
	else if (p.pid() == -13) {
          _h_pT_mu->fill(p.pT()/GeV, weight);
          _h_E_mu->fill(p.E()/GeV, weight);
          _h_eta_mu->fill(p.eta(), weight);
	}
	else if (p.pid() == 14) {
          _h_pT_numu->fill(p.pT()/GeV, weight);
          _h_E_numu->fill(p.E()/GeV, weight);
          _h_eta_numu->fill(p.eta(), weight);
	}
      }

      const Particles ps = fs.particles();
      // loop over all combinations
      for (unsigned int i=0; i<ps.size(); i++) {
        for (unsigned int j=0; j<ps.size(); j++) {
          if (ps[i].pid()==5 && ps[j].pid()==-5) {
            _h_angle_bb->fill(ps[i].angle(ps[j]), weight);
	  }
	  else if (ps[i].pid()==11 && ps[j].pid()==-12) {
            _h_angle_enu->fill(ps[i].angle(ps[j]), weight);
	    mom_Wm = ps[i].momentum() + ps[j].momentum();
            _h_Wm_mass->fill(mom_Wm.mass()/GeV, weight);
	    FourMomentum mom_tb = mom_Wm + mom_bb;
            _h_tb_mass->fill(mom_tb.mass()/GeV, weight);
	    _h_angle_Wmbb->fill(mom_Wm.angle(mom_bb), weight);
	  }
	  else if (ps[i].pid()==-13 && ps[j].pid()==14) {
            _h_angle_munu->fill(ps[i].angle(ps[j]), weight);
	    mom_Wp = ps[i].momentum() + ps[j].momentum();
            _h_Wp_mass->fill(mom_Wp.mass()/GeV, weight);
	    FourMomentum mom_t = mom_Wp + mom_b;
            _h_t_mass->fill(mom_t.mass()/GeV, weight);
	    _h_angle_Wpb->fill(mom_Wp.angle(mom_b), weight);
	  }
	  else if (ps[i].pid()==11 && ps[j].pid()==-5) {
            _h_angle_ebb->fill(ps[i].angle(ps[j]), weight);
	  }
	  else if (ps[i].pid()==-13 && ps[j].pid()==5) {
            _h_angle_mub->fill(ps[i].angle(ps[j]), weight);
	  }
	  else if (ps[i].pid()==-12 && ps[j].pid()==-5) {
            _h_angle_nubb->fill(ps[i].angle(ps[j]), weight);
	  }
	  else if (ps[i].pid()==14 && ps[j].pid()==5) {
            _h_angle_nub->fill(ps[i].angle(ps[j]), weight);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize
      normalize(_h_pT_b);
      normalize(_h_pT_bb);
      normalize(_h_pT_e);
      normalize(_h_pT_nue);
      normalize(_h_pT_mu);
      normalize(_h_pT_numu);
      normalize(_h_E_b);
      normalize(_h_E_bb);
      normalize(_h_E_e);
      normalize(_h_E_nue);
      normalize(_h_E_mu);
      normalize(_h_E_numu);
      normalize(_h_eta_b);
      normalize(_h_eta_bb);
      normalize(_h_eta_e);
      normalize(_h_eta_nue);
      normalize(_h_eta_mu);
      normalize(_h_eta_numu);
      normalize(_h_angle_bb);
      normalize(_h_angle_munu);
      normalize(_h_angle_enu);
      normalize(_h_angle_ebb);
      normalize(_h_angle_nubb);
      normalize(_h_angle_Wmbb);
      normalize(_h_angle_mub);
      normalize(_h_angle_nub);
      normalize(_h_angle_Wpb);
      normalize(_h_Wp_mass);
      normalize(_h_Wm_mass);
      normalize(_h_t_mass);
      normalize(_h_tb_mass);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_pT_b, _h_pT_bb, _h_pT_e, _h_pT_nue, _h_pT_mu, _h_pT_numu, _h_E_b, _h_E_bb, _h_E_e, _h_E_nue, _h_E_mu, _h_E_numu, _h_eta_b, _h_eta_bb, _h_eta_e, _h_eta_nue, _h_eta_mu, _h_eta_numu, _h_angle_bb, _h_angle_munu, _h_angle_enu, _h_angle_ebb, _h_angle_nubb, _h_angle_Wmbb, _h_angle_mub, _h_angle_nub, _h_angle_Wpb, _h_Wp_mass, _h_Wm_mass, _h_t_mass, _h_tb_mass;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TTBAR_DOUBLEDECAY);


}
