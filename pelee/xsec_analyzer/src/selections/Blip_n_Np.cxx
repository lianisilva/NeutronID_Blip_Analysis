// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/TreeUtils.hh"

#include "XSecAnalyzer/Selections/EventCategoriesBlip_n_Np.hh"   // NEED to make this file
#include "XSecAnalyzer/Selections/Blip_n_Np.hh"

Blip_n_Np::Blip_n_Np() : SelectionBase( "Blip_n_Np" ) {
  this->define_category_map();
  this->define_constants();
}

void Blip_n_Np::define_constants() {
  // Define reco and true fiducial volumes alongside any other constants used within selection cuts
  // x_min, x_max, y_min, y_max, z_min, z_max
  //this->define_true_FV( 10., 246.35, -106.5, 106.5, 10., 1026.8 );
  //this->define_reco_FV( 10., 246.35, -106.5, 106.5, 10., 1026.8 );
  this->define_FV( "exact", 0., 256., -116.25, 116.25, 0., 1036.8 );
  this->define_FV( "buffer", 10., 246.35, -106.5, 106.5, 10., 1026.8 );
  this->define_FV_with_dead_zone( "bufferAndDeadZone",
				  5., 251., -110., 110., 20., 986.,   //box
				  675., 775. );			      // dead zone in Z
}

void Blip_n_Np::compute_reco_observables( AnalysisEvent* event ) {
  // Calculate reconstructed kinematic variables to be saved in the output

}

void Blip_n_Np::compute_true_observables( AnalysisEvent* event ) {
  // Calculate true kinematic variables to be saved in the output

}

int Blip_n_Np::categorize_event( AnalysisEvent* event ) {
  // Assign the event category of the selected event
  return 0;
}

bool Blip_n_Np::define_signal( AnalysisEvent* event ) {
  // Determine whether an input MC event fulfills the signal definition.
  // Only truth information should be used to determine the answer.

  // NumuCC Inclusive Selection
  bool in_fv   = point_inside_FV( this->FV("bufferAndDeadZone"), event->mc_nu_vx_,
                                  event->mc_nu_vy_, event->mc_nu_vz_ );
  bool is_numu = ( event->mc_nu_pdg_ == MUON_NEUTRINO );
  bool is_cc   = ( event->mc_nu_ccnc_ == CHARGED_CURRENT );

  // we may need to do blip-level truth selection here but I don't think it's possible since, like PFPs, blips are reco. objects

  return in_fv && is_numu && is_cc;

}


bool Blip_n_Np::selection( AnalysisEvent* event ) {
  // Determine whether an input event satisfies the selection criteria.
  // Only reco information should be used to determine the answer.
    TVector3 reco_nu_vtx( event->reco_nu_vtx_sce_x_,
                          event->reco_nu_vtx_sce_y_,
                          event->reco_nu_vtx_sce_z_ );

    // Preselection
    sel_reco_vertex_in_FV_    = ( event->slice_orig_pass_id_ == 1 )
                                && point_inside_FV( this->FV("bufferAndDeadZone"),
                                                    event->reco_nu_vtx_sce_x_,
                                                    event->reco_nu_vtx_sce_y_,
                                                    event->reco_nu_vtx_sce_z_ );
    sel_topo_cut_passed_      = ( event->slice_orig_topo_score_ > 0.06 );
    sel_crt_veto_passed_      = ( event->crtveto_ != 1 );
    sel_cosmic_ip_cut_passed_ = ( event->_closestNuCosmicDist_ > 5.0 );

//std::cout << "FV=" << sel_reco_vertex_in_FV_
//          << " topo=" << sel_topo_cut_passed_
//          << " crt=" << sel_crt_veto_passed_
//          << " ip=" << sel_cosmic_ip_cut_passed_
//          << " pass_id=" << event->slice_orig_pass_id_
//          << " topo_score=" << event->slice_orig_topo_score_
//          << " crtveto=" << event->crtveto_
//          << " cosmicdist=" << event->_closestNuCosmicDist_
//          << "\n";

    if ( !sel_reco_vertex_in_FV_ || !sel_topo_cut_passed_
      || !sel_crt_veto_passed_   || !sel_cosmic_ip_cut_passed_ ) {
        return false;
    }

    // PFP loop — find muon candidate
    muon_candidate_idx_     = -1;
    sel_has_muon_candidate_ = false;

    for ( int p = 0; p < (int)event->n_pfps_; ++p ) {

        if ( event->trk_score_v_->at(p) <= 0.8 ) continue;//{std::cout<<"Continued after trk score check\n";continue;}	
        TVector3 trkStart( event->trk_sce_start_x_v_->at(p),
                           event->trk_sce_start_y_v_->at(p),
                           event->trk_sce_start_z_v_->at(p) );
        TVector3 trkEnd(   event->trk_sce_end_x_v_->at(p),
                           event->trk_sce_end_y_v_->at(p),
                           event->trk_sce_end_z_v_->at(p) );

        float trk_length  = (trkEnd - trkStart).Mag();
        float dist_to_vtx = std::min( (trkStart - reco_nu_vtx).Mag(),
                                      (trkEnd   - reco_nu_vtx).Mag() );

        if ( !point_inside_FV( this->FV("buffer"),
                               trkStart.X(), trkStart.Y(), trkStart.Z() ) ) continue;//{std::cout<<"Continued after FV check\n";continue;}
        if ( !point_inside_FV( this->FV("buffer"),
                               trkEnd.X(),   trkEnd.Y(),   trkEnd.Z()   ) ) continue;//{std::cout<<"Continued after FV check\n";continue;}
        if ( dist_to_vtx >= 4 )                                           continue;//{std::cout<<"Continued after vtx check\n";continue;}
        if ( event->track_llr_pid_score_->at(p) <= 0.2 )                    continue;//{std::cout<<"Continued after PID check\n";continue;}
        if ( trk_length <= 10 )                                           continue;//{std::cout<<"Continued after trk length check\n"; continue;}

        sel_has_muon_candidate_ = true;
        muon_candidate_idx_     = p;
        break;
    }

    sel_nu_mu_cc_ = sel_has_muon_candidate_;
    return sel_nu_mu_cc_;
}

void Blip_n_Np::define_output_branches() {
  // Call set_branch() for every new variable to be saved to the output TTree
  set_branch( &sel_reco_vertex_in_FV_,    "reco_vertex_in_FV" );//,    kBool );
  set_branch( &sel_topo_cut_passed_,      "topo_cut_passed" );//,      kBool );
  set_branch( &sel_crt_veto_passed_,      "crt_veto_passed" );//,      kBool );
  set_branch( &sel_cosmic_ip_cut_passed_, "cosmic_ip_cut_passed" );//, kBool );
  set_branch( &sel_has_muon_candidate_,   "has_muon_candidate" );//,   kBool );
  set_branch( &sel_nu_mu_cc_,             "nu_mu_cc" );//,             kBool );
  set_branch( &muon_candidate_idx_,       "muon_candidate_idx" );//,   kInteger );
}

void Blip_n_Np::reset() {
  // Set variables managed by this class to their default values. This function
  // is called in preparation for analyzing each new input event.
  sel_reco_vertex_in_FV_    = false;
  sel_topo_cut_passed_      = false;
  sel_crt_veto_passed_      = false;
  sel_cosmic_ip_cut_passed_ = false;
  sel_has_muon_candidate_   = false;
  sel_nu_mu_cc_             = false;
  muon_candidate_idx_       = -1;
}

void Blip_n_Np::define_category_map() {
  // Define the mapping between each integer event category and
  // a string label / color integer code pair
  // The color codes are documented at
  // https://root.cern.ch/doc/master/classTColor.html
  categ_map_ = Blip_n_Np_MAP;
}


