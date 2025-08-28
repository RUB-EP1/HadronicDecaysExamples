// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Amplitude analysis of the Lambda_c -> p K- pi+ decay
  class LHCB_2023_I2683025 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2023_I2683025);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid== 4122);
      declare(ufs, "UFS");
      DecayedParticles Lambda_c(ufs);
      Lambda_c.addStable(PID::PROTON);
      Lambda_c.addStable(PID::KMINUS);
      Lambda_c.addStable(PID::PIPLUS);
      declare(Lambda_c,"Lambda_c");
      // histograms
      book(_h_pK,  1,1,1);
      book(_h_Kpi, 2,1,1);
      book(_dalitz, 3,1,1);
      book(_dalitz_test, "dalitz",100,2.03,4.62,100,0.39,1.82);

    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles Lambda_c = apply<DecayedParticles>(event, "Lambda_c");
      // loop over particles
      for(unsigned int ix=0; ix<Lambda_c.decaying().size(); ++ix) {
	        const Particle & p  = Lambda_c.decayProducts()[ix].at(PID::PROTON)[0];
	        const Particle & K  = Lambda_c.decayProducts()[ix].at(PID::KMINUS)[0];
	        const Particle & pi = Lambda_c.decayProducts()[ix].at(PID::PIPLUS)[0];
	        double m_pK  = (p.momentum() + K.momentum()).mass2();
	        double m_Kpi = (K.momentum() + pi.momentum()).mass2();
	        _h_pK ->fill(m_pK);
	        _h_Kpi->fill(m_Kpi);
	        _dalitz->fill(m_pK, m_Kpi);
	        _dalitz_test->fill(m_pK, m_Kpi);
      }
    }




    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pK);
      normalize(_h_Kpi);
      normalize(_dalitz);
      normalize(_dalitz_test);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pK,_h_Kpi;
    Histo2DPtr _dalitz, _dalitz_test;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2023_I2683025);

}
