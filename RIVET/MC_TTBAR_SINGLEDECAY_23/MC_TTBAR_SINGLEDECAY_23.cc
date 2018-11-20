// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_TTBAR_SINGLEDECAY_23 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_TTBAR_SINGLEDECAY_23);


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
      _h_pT_W = bookHisto1D("pT_W", 50, 0.0, 600.0);
      _h_E_t = bookHisto1D("E_t", 50, 0.0, 600.0);
      _h_E_b = bookHisto1D("E_b", 50, 0.0, 600.0);
      _h_E_W = bookHisto1D("E_W", 50, 0.0, 600.0);
      _h_eta_t = bookHisto1D("eta_t", 50, -6.0, 6.0);
      _h_eta_b = bookHisto1D("eta_b", 50, -6.0, 6.0);
      _h_eta_W = bookHisto1D("eta_W", 50, -6.0, 6.0);
      _h_angle_tb = bookHisto1D("angle_tb", 50, 0, M_PI);
      _h_angle_tW = bookHisto1D("angle_tW", 50, 0, M_PI);
      _h_angle_bW = bookHisto1D("angle_bW", 50, 0, M_PI);
      _h_W_mass = bookHisto1D("W_mass", 50, 0.0, 200.0);
      _h_t_mass = bookHisto1D("Top_mass", 50, 0.0, 200.0);
      _h_tr_mass = bookHisto1D("recT_mass", 50, 0.0, 200.0);
      _h_phi_t = bookHisto1D("phi_t",50,0, M_PI);
      _h_phi_w = bookHisto1D("phi_w",50,0, M_PI);
      _h_phi_b = bookHisto1D("phi_b",50,0, M_PI);
      _h_b_mass = bookHisto1D("b_mass", 50, -1, 10);
      _h_angel_labt = bookHisto1D("angel_labt",50,0, M_PI);
      _h_angel_labw = bookHisto1D("angel_labw",50,0, M_PI);
      _h_angel_labb = bookHisto1D("angel_labb",50,0, M_PI);


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      FourMomentum mom_b;
      FourMomentum mom_W;
      FourMomentum mom_lab;
      mom_lab.setPx(0);
      mom_lab.setPy(0);
      mom_lab.setPz(1); 
      mom_lab.setE(1);

      const FinalState& fs = apply<FinalState>(event, "FS");
      for (const Particle& p : fs.particles()) {
        if (p.pid() == 6) {
          _h_pT_t->fill(p.pT()/GeV, weight);
          _h_E_t->fill(p.E()/GeV, weight);
          _h_eta_t->fill(p.eta(), weight);
          _h_t_mass->fill(p.mass()/GeV, weight);
          _h_phi_t->fill(p.phi(), weight);
          _h_angel_labt->fill(p.angle(mom_lab), weight);

	}
	else if (p.pid() == -5) {
          mom_b = p.momentum();
          _h_pT_b->fill(p.pT()/GeV, weight);
          _h_E_b->fill(p.E()/GeV, weight);
          _h_eta_b->fill(p.eta(), weight);
          _h_phi_b->fill(p.phi(), weight);
          _h_b_mass->fill(p.mass()/GeV, weight);
          _h_angel_labb->fill(p.angle(mom_lab), weight);



	}
	else if (p.pid() == -24) {
          mom_W = p.momentum();
          _h_pT_W->fill(p.pT()/GeV, weight);
          _h_E_W->fill(p.E()/GeV, weight);
          _h_eta_W->fill(p.eta(), weight);
          _h_W_mass->fill(p.mass()/GeV, weight);
          _h_phi_w->fill(p.phi(), weight);
          _h_angel_labw->fill(p.angle(mom_lab), weight);

	}
      }

      const Particles ps = fs.particles();
      // loop over all combinations
      for (unsigned int i=0; i<ps.size(); i++) {
        for (unsigned int j=0; j<ps.size(); j++) {
          if (ps[i].pid()==6 && ps[j].pid()==-5) {
            _h_angle_tb->fill(ps[i].angle(ps[j]), weight);
	  }
	  else if (ps[i].pid()==6 && ps[j].pid()==-24) {
            _h_angle_tW->fill(ps[i].angle(ps[j]), weight);
	  }
	  else if (ps[i].pid()==-5 && ps[j].pid()==-24) {
            _h_angle_bW->fill(ps[i].angle(ps[j]), weight);
            FourMomentum mom_t = ps[i].momentum() + ps[j].momentum();
            _h_tr_mass->fill(mom_t.mass()/GeV,weight);
	  }
	  // Hack for Sherpa because the top decays
	  else if (ps[i].pid()==5 && ps[j].pid()==24) {
            FourMomentum mom_t = ps[i].momentum() + ps[j].momentum();
	    _h_E_t->fill(mom_t.E()/GeV, weight);
            _h_pT_t->fill(mom_t.pT()/GeV, weight);
            _h_eta_t->fill(mom_t.eta(), weight);
            _h_angle_tb->fill(mom_t.angle(mom_b), weight);
            _h_angle_tW->fill(mom_t.angle(mom_W), weight);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize
      normalize(_h_pT_t);
      normalize(_h_pT_b);
      normalize(_h_pT_W);
      normalize(_h_E_t);
      normalize(_h_E_b);
      normalize(_h_E_W);
      normalize(_h_eta_t);
      normalize(_h_eta_b);
      normalize(_h_eta_W);
      normalize(_h_angle_tb);
      normalize(_h_angle_tW);
      normalize(_h_angle_bW);
      normalize(_h_W_mass);
      normalize(_h_t_mass);
      normalize(_h_tr_mass);
      normalize(_h_phi_t);
      normalize(_h_phi_b);
      normalize(_h_phi_w);
      normalize(_h_b_mass);
      normalize(_h_angel_labt);
      normalize(_h_angel_labw);
      normalize(_h_angel_labb);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_pT_t, _h_pT_b, _h_pT_W, _h_E_t, _h_E_b, _h_E_W, _h_eta_t, _h_eta_b, _h_eta_W, _h_angle_tb, _h_angle_tW, _h_angle_bW, _h_W_mass, _h_t_mass, _h_tr_mass, _h_phi_t, _h_phi_w, _h_phi_b, _h_b_mass, _h_angel_labt, _h_angel_labw,_h_angel_labb;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TTBAR_SINGLEDECAY_23);


}

