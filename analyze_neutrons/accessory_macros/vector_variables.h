
  //#ifndef OUTPUT_VARIABLES_H
  //#define OUTPUT_VARIABLES_H

  //#include <vector>
  //#include <string>

  // Event-level variables (scalars)
  int event;
  int run;
  int subrun;
  int nparticles;
  int nblips;
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
  std::vector<std::string> part_process;
  std::vector<bool> part_isPrimary;
  std::vector<int> part_g4id;
  std::vector<int> part_pdg;
  std::vector<int> part_nDaughters;
  std::vector<int> part_mother;
  std::vector<float> part_KE;
  std::vector<float> part_Px;
  std::vector<float> part_Py;
  std::vector<float> part_Pz;
  std::vector<float> part_startPointx;
  std::vector<float> part_startPointy;
  std::vector<float> part_startPointz;
  std::vector<float> part_endPointx;
  std::vector<float> part_endPointy;
  std::vector<float> part_endPointz;
  std::vector<float> part_pathlen;
  std::vector<float> part_depEnergy;

  std::vector<std::string> out_part_process;
  std::vector<bool> out_part_isPrimary;
  std::vector<int> out_part_g4id;
  std::vector<int> out_part_pdg;
  std::vector<int> out_part_nDaughters;
  std::vector<int> out_part_mother;
  std::vector<float> out_part_KE;
  std::vector<float> out_part_Px;
  std::vector<float> out_part_Py;
  std::vector<float> out_part_Pz;
  std::vector<float> out_part_startPointx;
  std::vector<float> out_part_startPointy;
  std::vector<float> out_part_startPointz;
  std::vector<float> out_part_endPointx;
  std::vector<float> out_part_endPointy;
  std::vector<float> out_part_endPointz;
  std::vector<float> out_part_pathlen;
  std::vector<float> out_part_depEnergy;

  // Blip-level variables (vectors - one entry per blip)
  std::vector<float> blip_x;
  std::vector<float> blip_y;
  std::vector<float> blip_z;
  std::vector<float> blip_size;
  std::vector<float> blip_energy;
  std::vector<int> blip_charge;
  std::vector<int> blip_edepid;
  std::vector<int> blip_g4id;
  std::vector<float> blip_dx;
  std::vector<float> blip_dw;
  std::vector<float> blip_proxtrkdist;
  std::vector<bool> blip_touchtrk;
  std::vector<int> blip_pl0_bydeadwire;
  std::vector<int> blip_pl1_bydeadwire;
  std::vector<int> blip_pl2_bydeadwire;

  std::vector<float> out_blip_x;
  std::vector<float> out_blip_y;
  std::vector<float> out_blip_z;
  std::vector<float> out_blip_size;
  std::vector<float> out_blip_energy;
  std::vector<int> out_blip_charge;
  std::vector<int> out_blip_edepid;
  std::vector<int> out_blip_g4id;
  std::vector<float> out_blip_dx;
  std::vector<float> out_blip_dw;
  std::vector<float> out_blip_proxtrkdist;
  std::vector<bool> out_blip_touchtrk;
  std::vector<int> out_blip_pl0_bydeadwire;
  std::vector<int> out_blip_pl1_bydeadwire;
  std::vector<int> out_blip_pl2_bydeadwire;

  // Helper function to clear vectors for new event
  void clear_branch_vars();

  //int out_nedeps;
  //float out_edep_energy;
  //bool out_edep_isPrimary;
  //float out_edep_x;
  //float out_edep_y;
  //float out_edep_z;
  //int out_edep_blipid; 

  //#endif

