// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"


namespace Rivet {


  /// First signal analysis pp->γjj->γbbjj
  class dijets : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(dijets);


    /// @name Analysis methods
    //@{

    /// Initialise projections before the run
    void init() {
      declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 30*GeV), "FS");
      FinalState fs;
      FastJets fj(fs, FastJets::ANTIKT, 0.4);
      declare(fj, "JET");


      // photon
      PromptFinalState photonfs(Cuts::abspid == PID::PHOTON && Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);
      declare(photonfs, "photons");



      ///Histos after cuts
      _h_m_dijet = bookHisto1D("dijet_mass", 100, 500, 2000, "signal_dijet_mass", "$m\textsubscript{j_1j_2}$", "$\frac{1}{\sigma}*\frac{d\sigma}{dm\textsubscript{j_1j_2}}$");
      _h_photon_pt  = bookHisto1D("photon_pT", 100, 0, 200, "signal_photon_pT", "photon_pT", "1/sigma*d_sigma/dpT");
      _h_m_bb = bookHisto1D("bjets_mass", 100, 80, 180, "signal_mass_b1b2", "m_b1b2", "1/sigma*dsigma/dmb1b2"); 
      _h_b1_pt = bookHisto1D("bjet1_pT", 100, 0, 500, "signal_b1_pT", "b1_pT", "1/sigma*d_sigma/dbb1pT");
      _h_b2_pt = bookHisto1D("bjet2_pT", 100, 30, 230, "signal_b2_pT", "b2_pT", "1/sigma*d_sigma/dbb2pT");
      _h_j1_pt = bookHisto1D("jet1_pT", 100, 40, 180 , "signal_j1_pT", "j1_pT", "1/sigma*d_sigma/dj1pT");
      _h_j2_pt = bookHisto1D("jet2_pT", 100, 40, 140 , "signal_j2_pT", "j2_pT", "1/sigma*d_sigma/dj2pT");
      _h_photon_eta = bookHisto1D("photon_eta", 30, 0, 3, "signal_photon_eta", "photon_eta", "1/sigma*d_sigma/dphoton_eta");
      _h_jj_deltaeta = bookHisto1D("jj_deltaeta", 30, 0, 6, "signal_deltaetajj", "delta_etaj1j2", "1/sigma*dsigma/d(deltaetaj1j2)");
      _h_bb_deltaeta = bookHisto1D("bb_deltaeta", 30, 0, 3, "signal_deltaetabb", "deltaeta_bb", "1/sigma*dsigma/d(deltaetab1b2)");
      _h_deltarap_jj = bookHisto1D("deltarap_jj", 40, 0, 4, "signal_deltarap_jj", "deltarap_jj", "1/sigma*dsigma/d(deltarap_jj)");
      _h_jet_mult = bookHisto1D("jet_multiplicity", 10, -0.5, 2.5, "signal_jet_mult", "N_jet", "1/sigma*dsigma/d(N_jet)");
      _h_deltaphi12 = bookHisto1D("deltaphi_jj", 20, 0, PI, "signal_dphi_jj", "dphi_jj", "1/sigma*dsigma/d(dphi)");
      _h_bjetmult = bookHisto1D("bjetmult", 10, -0.5, 2.5, "signal_bjetmult", "bjetmult", "1/sigma*dsigma/d(bjetmult)");
      _h_cut_flow = bookHisto1D("cut_flow",8, -0.5,7.5);
      _h_Zeppenfeld = bookHisto1D ("Zeppenfeld_variable", 10, -5, 5, "Zeppenfeld_variable", "Zeppenfeld(Z)", "$1/\sigma * d\sigma/dZ$");  
    }


    /// Perform the per-event analysis
    void analyze(const Event & evnt) {
      const double weight = evnt.weight();
      _h_cut_flow->fill(0, weight);
     // Get the photon
      const Particles& photons = apply<PromptFinalState>(evnt, "photons").particlesByPt(Cuts::abseta < 2.5);
      if (photons.empty())  vetoEvent;
      const  FourMomentum photon = photons[0].momentum();
      double photoneta = photons[0].abseta();


      Jets all_jets = apply<FastJets>(evnt, "JET").jetsByPt();
      Jets bjets;
      Jets extrajets;
	for ( Jet& j : all_jets){
                if ( j.bTagged() && j.pT()>30*GeV && j.abseta()<2.5) {bjets.push_back(j);}
                                             
                       else if (j.pT()>40*GeV && j.abseta()<4.5) {extrajets.push_back(j);} 
       };
 
       _h_cut_flow->fill(1, weight);
  
       if (bjets.size() < 2 || extrajets.size() < 2 || bjets.size() > 2 || extrajets.size() > 2 ) {vetoEvent;}

        _h_cut_flow->fill(2, weight);
       

for(Jet& b : bjets){
                    if(deltaR(b.eta(), b.phi(), photons[0].eta(), photons[0].phi()) < 0.4 ) {vetoEvent;}};

for(Jet& J : extrajets){
                    if(deltaR(J.eta(), J.phi(), photons[0].eta(), photons[0].phi()) < 0.4 ) {vetoEvent;}};
  
       _h_cut_flow->fill(3, weight);

     
      double bjet1pt = bjets[0].pT();
      double bjet2pt = bjets[1].pT();
      double extrajet1pt = extrajets[0].pT();
      double extrajet2pt = extrajets[1].pT();

      double bjet1eta = bjets[0].abseta();
      double bjets2eta = bjets[1].abseta();
      double bjetseta = abs(bjet1eta - bjets2eta);
      double extrajet1eta = extrajets[0].abseta();
      double extrajet2eta = extrajets[1].abseta();
      double deltaetaj1j2 = abs(extrajet1eta - extrajet2eta);
          
      const FourMomentum bj1(bjets[0].momentum());
      const FourMomentum bj2(bjets[1].momentum());
      const FourMomentum j1(extrajets[0].momentum());
      const FourMomentum j2(extrajets[1].momentum());
      const FourMomentum H(j1 + j2);
               double y_star= H.rapidity();
               if (y_star > 1/2*deltaetaj1j2) {vetoEvent;}
               _h_cut_flow->fill(4, weight);

      double mbb = FourMomentum(bj1+bj2).mass();
      double mjj = FourMomentum(j1+j2).mass();

      if (mjj < 600*GeV) {vetoEvent;}
       _h_cut_flow->fill(5, weight);

      if (mbb < 100*GeV || mbb > 140*GeV) {vetoEvent;}
        _h_cut_flow->fill(6, weight);

      double deltarap12 = j1.absrap() - j2.absrap();
      double deltaphi12 = deltaPhi(extrajets[0], extrajets[1]);

      // Fill histos
      double photon_pt = photons[0].pT();
      _h_photon_pt->fill(photon_pt, weight);
      _h_m_dijet->fill(mjj, weight);
      _h_m_bb->fill(mbb, weight);
      _h_b1_pt->fill(bjet1pt, weight);
      _h_b2_pt->fill(bjet2pt, weight);
      _h_j1_pt->fill(extrajet1pt, weight);
      _h_j2_pt->fill(extrajet2pt, weight);
      _h_photon_eta->fill(photoneta, weight);
      _h_jj_deltaeta->fill(deltaetaj1j2, weight);
      _h_bb_deltaeta->fill(bjetseta, weight);
      _h_jet_mult->fill(extrajets.size(), weight);
      _h_deltarap_jj->fill(deltarap12, weight);
      _h_bjetmult->fill(bjets.size(), weight);
      _h_deltaphi12->fill(deltaphi12, weight);
      _h_Zeppenfeld->fill(y_star, weight);
    }


    /// Finalize
    void finalize() {
      scale(_h_m_dijet, crossSection()/sumOfWeights());
      scale(_h_photon_pt, crossSection()/sumOfWeights());
      scale(_h_m_bb, crossSection()/sumOfWeights());
      scale(_h_b1_pt, crossSection()/sumOfWeights());
      scale(_h_b2_pt, crossSection()/sumOfWeights());
      scale(_h_j1_pt, crossSection()/sumOfWeights());
      scale(_h_j2_pt, crossSection()/sumOfWeights()); 
      scale(_h_photon_eta, crossSection()/sumOfWeights());
      scale(_h_jj_deltaeta, crossSection()/sumOfWeights());
      scale(_h_bb_deltaeta, crossSection()/sumOfWeights());
      scale(_h_jet_mult, crossSection()/sumOfWeights());
      scale( _h_deltarap_jj, crossSection()/sumOfWeights());
      scale(_h_deltaphi12,  crossSection()/sumOfWeights());
      scale(_h_bjetmult, crossSection()/sumOfWeights()); 
      scale(_h_cut_flow, crossSection()/sumOfWeights()); 
      scale(_h_Zeppenfeld, crossSection()/sumOfWeights());


}
    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_photon_pt;
    Histo1DPtr _h_m_dijet;
    Histo1DPtr _h_m_bb;
    Histo1DPtr _h_b1_pt;
    Histo1DPtr _h_b2_pt;
    Histo1DPtr _h_j1_pt;
    Histo1DPtr _h_j2_pt;
    Histo1DPtr _h_photon_eta;
    Histo1DPtr _h_jj_deltaeta;
    Histo1DPtr _h_bb_deltaeta;
    Histo1DPtr _h_jet_mult;
    Histo1DPtr _h_deltarap_jj;
    Histo1DPtr _h_deltaphi12;
    Histo1DPtr _h_bjetmult;
    Histo1DPtr _h_cut_flow;
    Histo1DPtr _h_Zeppenfeld;

    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(dijets);


}
