
  //#ifndef OUTPUT_VARIABLES_H
  //#define OUTPUT_VARIABLES_H

  //#include <vector>
  //#include <string>

  const int kMaxBlips = 1000;
  const int kMaxG4    = 5000;
  const int kMaxEDeps = 5000;

  // Event-level variables (scalars)
  int event;
  int run;
  int subrun;
  int nparticles;
  int nblips;
  int nedeps;
  float nu_reco_vtx_x;
  float nu_reco_vtx_y;
  float nu_reco_vtx_z;
  float mctruth_nu_vtx_x;
  float mctruth_nu_vtx_y;
  float mctruth_nu_vtx_z;
  int mctruth_nu_pdg;
  int mctruth_nu_ccnc;

  int out_event;
  int out_run;
  int out_subrun;
  int out_nparticles;
  int out_nblips;
  float out_nu_reco_vtx_x;
  float out_nu_reco_vtx_y;
  float out_nu_reco_vtx_z;
  float out_mctruth_nu_vtx_x;
  float out_mctruth_nu_vtx_y;
  float out_mctruth_nu_vtx_z;
  int out_mctruth_nu_pdg;
  int out_mctruth_nu_ccnc;

  // Particle-level variables (vectors - one entry per particle)
  std::vector<std::string> *part_process;
  bool part_isPrimary[kMaxG4];
  int part_g4id[kMaxG4];
  int part_pdg[kMaxG4];
  int part_nDaughters[kMaxG4];
  int part_mother[kMaxG4];
  float part_KE[kMaxG4];
  float part_Px[kMaxG4];
  float part_Py[kMaxG4];
  float part_Pz[kMaxG4];
  float part_startPointx[kMaxG4];
  float part_startPointy[kMaxG4];
  float part_startPointz[kMaxG4];
  float part_endPointx[kMaxG4];
  float part_endPointy[kMaxG4];
  float part_endPointz[kMaxG4];
  float part_pathlen[kMaxG4];
  float part_depEnergy[kMaxG4];

  std::vector<std::string> *out_part_process;
  bool out_part_isPrimary[kMaxG4];
  int out_part_g4id[kMaxG4];
  int out_part_pdg[kMaxG4];
  int out_part_nDaughters[kMaxG4];
  int out_part_mother[kMaxG4];
  float out_part_KE[kMaxG4];
  float out_part_Px[kMaxG4];
  float out_part_Py[kMaxG4];
  float out_part_Pz[kMaxG4];
  float out_part_startPointx[kMaxG4];
  float out_part_startPointy[kMaxG4];
  float out_part_startPointz[kMaxG4];
  float out_part_endPointx[kMaxG4];
  float out_part_endPointy[kMaxG4];
  float out_part_endPointz[kMaxG4];
  float out_part_pathlen[kMaxG4];
  float out_part_depEnergy[kMaxG4];

  // Edep-level variables (vectors - one entry per edep)
  float edep_energy[kMaxEDeps];
  bool edep_isPrimary[kMaxEDeps];
  float edep_x[kMaxEDeps];
  float edep_y[kMaxEDeps];
  float edep_z[kMaxEDeps];
  int edep_blipid[kMaxEDeps];

  int out_nedeps[kMaxEDeps];
  float out_edep_energy[kMaxEDeps];
  bool out_edep_isPrimary[kMaxEDeps];
  float out_edep_x[kMaxEDeps];
  float out_edep_y[kMaxEDeps];
  float out_edep_z[kMaxEDeps];
  int out_edep_blipid[kMaxEDeps];

  // Blip-level variables (vectors - one entry per blip)
  int   blip_id[kMaxBlips];
  float blip_x[kMaxBlips];
  float blip_y[kMaxBlips];
  float blip_z[kMaxBlips];
  float blip_size[kMaxBlips];
  float blip_energy[kMaxBlips];
  int blip_charge[kMaxBlips];
  int blip_edepid[kMaxBlips];
  int blip_g4id[kMaxBlips];
  float blip_dx[kMaxBlips];
  float blip_dw[kMaxBlips];
  float blip_proxtrkdist[kMaxBlips];
  bool blip_touchtrk[kMaxBlips];
  bool blip_pl0_bydeadwire[kMaxBlips];
  bool blip_pl1_bydeadwire[kMaxBlips];
  bool blip_pl2_bydeadwire[kMaxBlips];

  float out_blip_x[kMaxBlips];
  float out_blip_y[kMaxBlips];
  float out_blip_z[kMaxBlips];
  float out_blip_size[kMaxBlips];
  float out_blip_energy[kMaxBlips];
  int out_blip_charge[kMaxBlips];
  int out_blip_edepid[kMaxBlips];
  int out_blip_g4id[kMaxBlips];
  float out_blip_dx[kMaxBlips];
  float out_blip_dw[kMaxBlips];
  float out_blip_proxtrkdist[kMaxBlips];
  bool out_blip_touchtrk[kMaxBlips];
  bool out_blip_pl0_bydeadwire[kMaxBlips];
  bool out_blip_pl1_bydeadwire[kMaxBlips];
  bool out_blip_pl2_bydeadwire[kMaxBlips];
 
  void reset_vars() {

    event=0;
    run=0;
    subrun=0;
    nparticles=0;
    nblips=0;
    nu_reco_vtx_x=0;
    nu_reco_vtx_y=0;
    nu_reco_vtx_z=0;
    mctruth_nu_vtx_x=0;
    mctruth_nu_vtx_y=0;
    mctruth_nu_vtx_z=0;
    mctruth_nu_pdg=0;
    mctruth_nu_ccnc=0;

    out_event=0;
    out_run=0;
    out_subrun=0;
    out_nparticles=0;
    out_nblips=0;
    out_nu_reco_vtx_x=0;
    out_nu_reco_vtx_y=0;
    out_nu_reco_vtx_z=0;
    out_mctruth_nu_vtx_x=0;
    out_mctruth_nu_vtx_y=0;
    out_mctruth_nu_vtx_z=0;
    out_mctruth_nu_pdg=0;
    out_mctruth_nu_ccnc=0;

    for (int i=0; i<kMaxG4; i++) {
      //part_process->at(i)="";
      part_isPrimary[i]=0;
      part_g4id[i]=0;
      part_pdg[i]=0;
      part_nDaughters[i]=0;
      part_mother[i]=0;
      part_KE[i]=0;
      part_Px[i]=0;
      part_Py[i]=0;
      part_Pz[i]=0;
      part_startPointx[i]=0;
      part_startPointy[i]=0;
      part_startPointz[i]=0;
      part_endPointx[i]=0;
      part_endPointy[i]=0;
      part_endPointz[i]=0;
      part_pathlen[i]=0;
      part_depEnergy[i]=0;

      //out_part_process->at(i)="";
      out_part_isPrimary[i]=0;
      out_part_g4id[i]=0;
      out_part_pdg[i]=0;
      out_part_nDaughters[i]=0;
      out_part_mother[i]=0;
      out_part_KE[i]=0;
      out_part_Px[i]=0;
      out_part_Py[i]=0;
      out_part_Pz[i]=0;
      out_part_startPointx[i]=0;
      out_part_startPointy[i]=0;
      out_part_startPointz[i]=0;
      out_part_endPointx[i]=0;
      out_part_endPointy[i]=0;
      out_part_endPointz[i]=0;
      out_part_pathlen[i]=0;
      out_part_depEnergy[i]=0;
    }

    for (int i=0; i<kMaxBlips; i++) {
      blip_x[i]=0;
      blip_y[i]=0;
      blip_z[i]=0;
      blip_size[i]=0;
      blip_energy[i]=0;
      blip_charge[i]=0;
      blip_edepid[i]=0;
      blip_g4id[i]=0;
      blip_dx[i]=0;
      blip_dw[i]=0;
      blip_proxtrkdist[i]=0;
      blip_touchtrk[i]=0;
      blip_pl0_bydeadwire[i]=0;
      blip_pl1_bydeadwire[i]=0;
      blip_pl2_bydeadwire[i]=0;

      out_blip_x[i]=0;
      out_blip_y[i]=0;
      out_blip_z[i]=0;
      out_blip_size[i]=0;
      out_blip_energy[i]=0;
      out_blip_charge[i]=0;
      out_blip_edepid[i]=0;
      out_blip_g4id[i]=0;
      out_blip_dx[i]=0;
      out_blip_dw[i]=0;
      out_blip_proxtrkdist[i]=0;
      out_blip_touchtrk[i]=0;
      out_blip_pl0_bydeadwire[i]=0;
      out_blip_pl1_bydeadwire[i]=0;
      out_blip_pl2_bydeadwire[i]=0;
    }
  }

  // Helper function to clear vectors for new event
  //void clear_branch_vars(); 

  //#endif

