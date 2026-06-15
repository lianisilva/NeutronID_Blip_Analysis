  #pragma once
  
  // XSecAnalyzer includes
  #include "XSecAnalyzer/Selections/SelectionBase.hh"
  
  class Blip_n_Np : public SelectionBase {
  
  public:
  
  Blip_n_Np();

  virtual void define_constants();
  virtual int categorize_event( AnalysisEvent* event );
  virtual bool selection( AnalysisEvent* event );
  virtual bool define_signal( AnalysisEvent* event );
  virtual void compute_reco_observables( AnalysisEvent* event );
  virtual void compute_true_observables( AnalysisEvent* event );
  virtual void define_output_branches();
  virtual void define_category_map();
  virtual void reset();

private:

  // Preselection flags
  bool sel_reco_vertex_in_FV_      = false;
  bool sel_topo_cut_passed_        = false;
  bool sel_crt_veto_passed_        = false;
  bool sel_cosmic_ip_cut_passed_   = false;

  // Muon candidate flag
  bool sel_has_muon_candidate_     = false;
  bool sel_nu_mu_cc_               = false;

  int  muon_candidate_idx_         = -1;

};
