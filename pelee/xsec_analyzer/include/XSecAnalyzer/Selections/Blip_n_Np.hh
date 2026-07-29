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

  // numuCC preselection flags
  bool sel_reco_vertex_in_FV_      = false;
  bool sel_topo_cut_passed_        = false;
  bool sel_crt_veto_passed_        = false;
  bool sel_cosmic_ip_cut_passed_   = false;

  // muon candidate flag
  bool sel_has_muon_candidate_     = false;

  bool sel_nu_mu_cc_               = false;

  // blip preselection flags
  bool blip_not_touching_trk_	   = false;
  bool blip_far_from_trk_	   = false;
  bool blip_not_by_dead_wire_	   = false;
  bool blip_inside_FV_		   = false;
  bool blip_far_from_vtx_	   = false;
  bool blip_outside_shr_	   = false;

  int  muon_candidate_idx_         = -1;
  MyPointer<std::vector<int>> blip_idx_v_;

};
