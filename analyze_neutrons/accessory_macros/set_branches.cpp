
  //#include "in_out_variables.h"
  #include "../variables.h"
  #define SET_BRANCH(tree, name) tree->SetBranchAddress(#name, &name)

  //void clear_branch_vars() {
  //  out_part_process.clear();
  //  out_part_isPrimary.clear();
  //  out_part_g4id.clear();
  //  out_part_pdg.clear();
  //  out_part_nDaughters.clear();
  //  out_part_mother.clear();
  //  out_part_KE.clear();
  //  out_part_Px.clear();
  //  out_part_Py.clear();
  //  out_part_Pz.clear();
  //  out_part_startPointx.clear();
  //  out_part_startPointy.clear();
  //  out_part_startPointz.clear();
  //  out_part_endPointx.clear();
  //  out_part_endPointy.clear();
  //  out_part_endPointz.clear();
  //  out_part_pathlen.clear();
  //  out_part_depEnergy.clear();
  //  
  //  out_blip_x.clear();
  //  out_blip_y.clear();
  //  out_blip_z.clear();
  //  out_blip_size.clear();
  //  out_blip_energy.clear();
  //  out_blip_charge.clear();
  //  out_blip_edepid.clear();
  //  out_blip_g4id.clear();
  //  out_blip_dx.clear();
  //  out_blip_dw.clear();
  //  out_blip_proxtrkdist.clear();
  //  out_blip_touchtrk.clear();
  //  out_blip_pl0_bydeadwire.clear();
  //  out_blip_pl1_bydeadwire.clear();
  //  out_blip_pl2_bydeadwire.clear();
  //}

  // Set branch address for all variables you plan to use
  // (if a branch address isn't set, the variables will 
  // not be filled when each event is read). For a description 
  // of each variable, see the comments in 'variables.h'
  void get_branches(TTree* tree, bool isOutTree) {
    //if(isOutTree) {
    //  SET_BRANCH(tree, out_event);
    //  SET_BRANCH(tree, out_run);
    //  SET_BRANCH(tree, out_subrun);
    //  SET_BRANCH(tree, out_nparticles);
    //  SET_BRANCH(tree, out_part_process);
    //  SET_BRANCH(tree, out_part_isPrimary);
    //  SET_BRANCH(tree, out_part_g4id);
    //  SET_BRANCH(tree, out_part_pdg);
    //  SET_BRANCH(tree, out_part_mother);
    //  SET_BRANCH(tree, out_part_KE);
    //  SET_BRANCH(tree, out_part_pathlen);
    //  SET_BRANCH(tree, out_part_depEnergy);
    //  SET_BRANCH(tree, out_part_nDaughters);
    //  SET_BRANCH(tree, out_part_Px);
    //  SET_BRANCH(tree, out_part_Py);
    //  SET_BRANCH(tree, out_part_Pz);
    //  SET_BRANCH(tree, out_part_startPointx);
    //  SET_BRANCH(tree, out_part_startPointy);
    //  SET_BRANCH(tree, out_part_startPointz);
    //  SET_BRANCH(tree, out_part_endPointx);
    //  SET_BRANCH(tree, out_part_endPointy);
    //  SET_BRANCH(tree, out_part_endPointz);
    //  SET_BRANCH(tree, out_nblips);
    //  SET_BRANCH(tree, out_blip_x);
    //  SET_BRANCH(tree, out_blip_y);
    //  SET_BRANCH(tree, out_blip_z);
    //  SET_BRANCH(tree, out_blip_size);
    //  SET_BRANCH(tree, out_blip_energy);
    //  SET_BRANCH(tree, out_blip_charge);
    //  SET_BRANCH(tree, out_blip_edepid);
    //  SET_BRANCH(tree, out_blip_g4id);
    //  SET_BRANCH(tree, out_blip_dx);
    //  SET_BRANCH(tree, out_blip_dw);
    //  SET_BRANCH(tree, out_blip_proxtrkdist);
    //  SET_BRANCH(tree, out_blip_touchtrk);
    //  SET_BRANCH(tree, out_blip_pl0_bydeadwire);
    //  SET_BRANCH(tree, out_blip_pl1_bydeadwire);
    //  SET_BRANCH(tree, out_blip_pl2_bydeadwire);
    //  SET_BRANCH(tree, out_nu_reco_vtx_x);
    //  SET_BRANCH(tree, out_nu_reco_vtx_y);
    //  SET_BRANCH(tree, out_nu_reco_vtx_z);
    //  SET_BRANCH(tree, out_mctruth_nu_vtx_x);
    //  SET_BRANCH(tree, out_mctruth_nu_vtx_y);
    //  SET_BRANCH(tree, out_mctruth_nu_vtx_z);
    //  SET_BRANCH(tree, out_mctruth_nu_pdg);
    //  SET_BRANCH(tree, out_mctruth_nu_ccnc);

    //} else {                
      tree->SetBranchAddress("event",		&event);
      tree->SetBranchAddress("run",		&run);
      tree->SetBranchAddress("subrun",		&subrun);
      tree->SetBranchAddress("nparticles",	&nparticles);
      tree->SetBranchAddress("part_process",	&part_process);
      tree->SetBranchAddress("part_isPrimary",	&part_isPrimary);
      tree->SetBranchAddress("part_isContained",&part_isContained);
      tree->SetBranchAddress("part_g4id",	&part_g4id);
      tree->SetBranchAddress("part_pdg",	&part_pdg);
      tree->SetBranchAddress("part_mother",	&part_mother);
      tree->SetBranchAddress("part_KE",		&part_KE);
      tree->SetBranchAddress("part_pathlen",	&part_pathlen);
      tree->SetBranchAddress("part_depEnergy",	&part_depEnergy);
      tree->SetBranchAddress("part_nDaughters",	&part_nDaughters);
      tree->SetBranchAddress("part_Px",		&part_Px);
      tree->SetBranchAddress("part_Py",		&part_Py);
      tree->SetBranchAddress("part_Pz",		&part_Pz);
      tree->SetBranchAddress("part_startPointx",&part_startPointx);
      tree->SetBranchAddress("part_startPointy",&part_startPointy);
      tree->SetBranchAddress("part_startPointz",&part_startPointz);
      tree->SetBranchAddress("part_endPointx",	&part_endPointx);
      tree->SetBranchAddress("part_endPointy",	&part_endPointy);
      tree->SetBranchAddress("part_endPointz",	&part_endPointz);
      //tree->SetBranchAddress("part_startT",     &part_startT);
      tree->SetBranchAddress("nedeps",		&nedeps);
      tree->SetBranchAddress("edep_energy",	&edep_energy);
      tree->SetBranchAddress("edep_isPrimary",	&edep_isPrimary);
      tree->SetBranchAddress("edep_x",		&edep_x);
      tree->SetBranchAddress("edep_y",		&edep_y);
      tree->SetBranchAddress("edep_z",		&edep_z);
      tree->SetBranchAddress("edep_blipid",	&edep_blipid);
      tree->SetBranchAddress("nblips",		&nblips);
      tree->SetBranchAddress("blip_id",		&blip_id);
      tree->SetBranchAddress("blip_x",		&blip_x);
      tree->SetBranchAddress("blip_y",		&blip_y);
      tree->SetBranchAddress("blip_z",		&blip_z);
      tree->SetBranchAddress("blip_size",	&blip_size);
      tree->SetBranchAddress("blip_energy",	&blip_energy);
      tree->SetBranchAddress("blip_charge",	&blip_charge);
      tree->SetBranchAddress("blip_edepid",	&blip_edepid);
      tree->SetBranchAddress("blip_g4id",	&blip_g4id);
      tree->SetBranchAddress("blip_dx",		&blip_dx);
      tree->SetBranchAddress("blip_dw",		&blip_dw);
      tree->SetBranchAddress("blip_proxtrkdist",&blip_proxtrkdist);
      tree->SetBranchAddress("blip_touchtrk",	&blip_touchtrk);
      tree->SetBranchAddress("blip_pl0_bydeadwire",&blip_pl0_bydeadwire);
      tree->SetBranchAddress("blip_pl1_bydeadwire",&blip_pl1_bydeadwire);
      tree->SetBranchAddress("blip_pl2_bydeadwire",&blip_pl2_bydeadwire);
      tree->SetBranchAddress("ntrks",           &ntrks);
      tree->SetBranchAddress("trk_length",      &trk_length);
      tree->SetBranchAddress("trk_g4id",        &trk_g4id);
      tree->SetBranchAddress("nu_crt_closestCosDist",&nu_crt_closestCosDist);
      tree->SetBranchAddress("nu_reco_vtx_x",	&nu_reco_vtx_x);
      tree->SetBranchAddress("nu_reco_vtx_y",	&nu_reco_vtx_y);
      tree->SetBranchAddress("nu_reco_vtx_z",	&nu_reco_vtx_z);
      tree->SetBranchAddress("mctruth_nu_vtx_x",&mctruth_nu_vtx_x);
      tree->SetBranchAddress("mctruth_nu_vtx_y",&mctruth_nu_vtx_y);
      tree->SetBranchAddress("mctruth_nu_vtx_z",&mctruth_nu_vtx_z);
      tree->SetBranchAddress("mctruth_nu_pdg",	&mctruth_nu_pdg);
      tree->SetBranchAddress("mctruth_nu_ccnc",	&mctruth_nu_ccnc);
      tree->SetBranchAddress("nhits",		&nhits);
      tree->SetBranchAddress("hit_blipid",	&hit_blipid);
      tree->SetBranchAddress("hit_plane",	&hit_plane);
      tree->SetBranchAddress("hit_wire",	&hit_wire);
      tree->SetBranchAddress("hit_peakT",	&hit_peakT);
      tree->SetBranchAddress("hit_driftT",	&hit_driftT);
      tree->SetBranchAddress("hit_rms",		&hit_rms);
      tree->SetBranchAddress("hit_amp",		&hit_amp);
      //tree->SetBranchAddress("hit_plane",     &hit_plane);
    //}
  }


  // function for when you make your own
  // tree to fill, not just use the branches
  // of a premade TTree like above
  //void define_branches(TTree* tree, bool isOutTree) {
  //  if( isOutTree ) {
  //    tree->Branch("event",&event,"event/I");
  //    tree->Branch("run",&run,"run/I");
  //    tree->Branch("subrun",&subrun,"subrun/I");
  //    tree->Branch("nparticles",&nparticles,"nparticles/I");
  //    tree->Branch("nblips",&nblips,"nblips/I");
  //    tree->Branch("nu_reco_vtx_x",&nu_reco_vtx_x,"nu_reco_vtx_x/F");
  //    tree->Branch("nu_reco_vtx_y",&nu_reco_vtx_y,"nu_reco_vtx_y/F");
  //    tree->Branch("nu_reco_vtx_z",&nu_reco_vtx_z,"nu_reco_vtx_z/F");
  //    tree->Branch("mctruth_nu_vtx_x",&mctruth_nu_vtx_x,"mctruth_nu_vtx_x/F");
  //    tree->Branch("mctruth_nu_vtx_y",&mctruth_nu_vtx_y,"mctruth_nu_vtx_y/F");
  //    tree->Branch("mctruth_nu_vtx_z",&mctruth_nu_vtx_z,"mctruth_nu_vtx_z/F");
  //    tree->Branch("mctruth_nu_pdg",&mctruth_nu_pdg,"mctruth_nu_pdg/I");
  //    tree->Branch("mctruth_nu_ccnc",&mctruth_nu_ccnc,"mctruth_nu_ccnc/I");

  //    tree->Branch("part_process",&part_process);   // do I need to do something special for this one?
  //    tree->Branch("part_isPrimary",part_isPrimary,"part_isPrimary[nparticles]/O"); 
  //    tree->Branch("part_g4id",part_g4id,"part_g4id[nparticles]/I");
  //    tree->Branch("part_pdg",part_pdg,"part_pdg[nparticles]/I");
  //    tree->Branch("part_mother",part_mother,"part_mother[nparticles]/I");
  //    tree->Branch("part_KE",part_KE,"part_KE[nparticles]/F");
  //    tree->Branch("part_pathlen",part_pathlen,"part_pathlen[nparticles]/F");
  //    tree->Branch("part_depEnergy",part_depEnergy,"part_depEnergy[nparticles]/F");
  //    tree->Branch("part_nDaughters",part_nDaughters,"part_nDaughters[nparticles]/I");
  //    tree->Branch("part_Px",part_Px,"part_Px[nparticles]/F");
  //    tree->Branch("part_Py",part_Py,"part_Py[nparticles]/F");
  //    tree->Branch("part_Pz",part_Pz,"part_Pz[nparticles]/F");
  //    tree->Branch("part_startPointx",part_startPointx,"part_startPointx[nparticles]/F");
  //    tree->Branch("part_startPointy",part_startPointy,"part_startPointy[nparticles]/F");
  //    tree->Branch("part_startPointz",part_startPointz,"part_startPointz[nparticles]/F");
  //    tree->Branch("part_endPointx",part_endPointx,"part_endPointx[nparticles]/F");
  //    tree->Branch("part_endPointy",part_endPointy,"part_endPointy[nparticles]/F");
  //    tree->Branch("part_endPointz",part_endPointz,"part_endPointz[nparticles]/F");

  //    tree->Branch("blip_x",blip_x,"blip_x[nblips]/F");
  //    tree->Branch("blip_y",blip_y,"blip_y[nblips]/F");
  //    tree->Branch("blip_z",blip_z,"blip_z[nblips]/F");
  //    tree->Branch("blip_size",blip_size,"blip_size[nblips]/F");
  //    tree->Branch("blip_energy",blip_energy,"blip_energy[nblips]/F");
  //    tree->Branch("blip_charge",blip_charge,"blip_charge[nblips]/I");
  //    tree->Branch("blip_edepid",blip_edepid,"blip_edepid[nblips]/I");
  //    tree->Branch("blip_g4id",blip_g4id,"blip_g4id[nblips]/I");
  //    tree->Branch("blip_dx",blip_dx,"blip_dx[nblips]/F");
  //    tree->Branch("blip_dw",blip_dw,"blip_dw[nblips]/F"); 
  //    tree->Branch("blip_proxtrkdist",blip_proxtrkdist,"blip_proxtrkdist[nblips]/F");
  //    tree->Branch("blip_touchtrk",blip_touchtrk,"blip_touchtrk[nblips]/O");
  //    tree->Branch("blip_pl0_bydeadwire",blip_pl0_bydeadwire,"blip_pl0_bydeadwire[nblips]/I");
  //    tree->Branch("blip_pl1_bydeadwire",blip_pl1_bydeadwire,"blip_pl1_bydeadwire[nblips]/I");
  //    tree->Branch("blip_pl2_bydeadwire",blip_pl2_bydeadwire,"blip_pl2_bydeadwire[nblips]/I");

  //  } else {

  //    tree->Branch("out_event",&out_event,"out_event/I");
  //    tree->Branch("out_run",&out_run,"out_run/I");
  //    tree->Branch("out_subrun",&out_subrun,"out_subrun/I");
  //    tree->Branch("out_nparticles",&out_nparticles,"out_nparticles/I");
  //    tree->Branch("out_nblips",&out_nblips,"out_nblips/I");
  //    tree->Branch("out_nu_reco_vtx_x",&out_nu_reco_vtx_x,"out_nu_reco_vtx_x/F");
  //    tree->Branch("out_nu_reco_vtx_y",&out_nu_reco_vtx_y,"out_nu_reco_vtx_y/F");
  //    tree->Branch("out_nu_reco_vtx_z",&out_nu_reco_vtx_z,"out_nu_reco_vtx_z/F");
  //    tree->Branch("out_mctruth_nu_vtx_x",&out_mctruth_nu_vtx_x,"out_mctruth_nu_vtx_x/F");
  //    tree->Branch("out_mctruth_nu_vtx_y",&out_mctruth_nu_vtx_y,"out_mctruth_nu_vtx_y/F");
  //    tree->Branch("out_mctruth_nu_vtx_z",&out_mctruth_nu_vtx_z,"out_mctruth_nu_vtx_z/F");
  //    tree->Branch("out_mctruth_nu_pdg",&out_mctruth_nu_pdg,"out_mctruth_nu_pdg/I");
  //    tree->Branch("out_mctruth_nu_ccnc",&out_mctruth_nu_ccnc,"out_mctruth_nu_ccnc/I");

  //    tree->Branch("out_part_process",&out_part_process);
  //    tree->Branch("out_part_isPrimary",out_part_isPrimary,"out_part_isPrimary[out_nparticles]/O");
  //    tree->Branch("out_part_g4id",out_part_g4id,"out_part_g4id[out_nparticles]/I");
  //    tree->Branch("out_part_pdg",out_part_pdg,"out_part_pdg[out_nparticles]/I");
  //    tree->Branch("out_part_mother",out_part_mother,"out_part_mother[out_nparticles]/I");
  //    tree->Branch("out_part_KE",out_part_KE,"out_part_KE[out_nparticles]/F");
  //    tree->Branch("out_part_pathlen",out_part_pathlen,"out_part_pathlen[out_nparticles]/F");
  //    tree->Branch("out_part_depEnergy",out_part_depEnergy,"out_part_depEnergy[out_nparticles]/F");
  //    tree->Branch("out_part_nDaughters",out_part_nDaughters,"out_part_nDaughters[out_nparticles]/I");
  //    tree->Branch("out_part_Px",out_part_Px,"out_part_Px[out_nparticles]/F");
  //    tree->Branch("out_part_Py",out_part_Py,"out_part_Py[out_nparticles]/F");
  //    tree->Branch("out_part_Pz",out_part_Pz,"out_part_Pz[out_nparticles]/F");
  //    tree->Branch("out_part_startPointx",out_part_startPointx,"out_part_startPointx[out_nparticles]/F");
  //    tree->Branch("out_part_startPointy",out_part_startPointy,"out_part_startPointy[out_nparticles]/F");
  //    tree->Branch("out_part_startPointz",out_part_startPointz,"out_part_startPointz[out_nparticles]/F");
  //    tree->Branch("out_part_endPointx",out_part_endPointx,"out_part_endPointx[out_nparticles]/F");
  //    tree->Branch("out_part_endPointy",out_part_endPointy,"out_part_endPointy[out_nparticles]/F");
  //    tree->Branch("out_part_endPointz",out_part_endPointz,"out_part_endPointz[out_nparticles]/F");
  //    
  //    tree->Branch("out_blip_x",out_blip_x,"out_blip_x[out_nblips]/F");
  //    tree->Branch("out_blip_y",out_blip_y,"out_blip_y[out_nblips]/F");
  //    tree->Branch("out_blip_z",out_blip_z,"out_blip_z[out_nblips]/F");
  //    tree->Branch("out_blip_size",out_blip_size,"out_blip_size[out_nblips]/F");
  //    tree->Branch("out_blip_energy",out_blip_energy,"out_blip_energy[out_nblips]/F");
  //    tree->Branch("out_blip_charge",out_blip_charge,"out_blip_charge[out_nblips]/I");
  //    tree->Branch("out_blip_edepid",out_blip_edepid,"out_blip_edepid[out_nblips]/I");
  //    tree->Branch("out_blip_g4id",out_blip_g4id,"out_blip_g4id[out_nblips]/I");
  //    tree->Branch("out_blip_dx",out_blip_dx,"out_blip_dx[out_nblips]/F");
  //    tree->Branch("out_blip_dw",out_blip_dw,"out_blip_dw[out_nblips]/F");
  //    tree->Branch("out_blip_proxtrkdist",out_blip_proxtrkdist,"out_blip_proxtrkdist[out_nblips]/F");
  //    tree->Branch("out_blip_touchtrk",out_blip_touchtrk,"out_blip_touchtrk[out_nblips]/O");
  //    tree->Branch("out_blip_pl0_bydeadwire", out_blip_pl0_bydeadwire, "out_blip_pl0_bydeadwire[out_nblips]/I");
  //    tree->Branch("out_blip_pl1_bydeadwire", out_blip_pl1_bydeadwire, "out_blip_pl1_bydeadwire[out_nblips]/I");
  //    tree->Branch("out_blip_pl2_bydeadwire", out_blip_pl2_bydeadwire, "out_blip_pl2_bydeadwire[out_nblips]/I");

  //    //tree->Branch("out_nedeps", &out_nedeps);
  //    //tree->Branch("out_edep_energy", &out_edep_energy);
  //    //tree->Branch("out_edep_isPrimary", &out_edep_isPrimary);
  //    //tree->Branch("out_edep_x", &out_edep_x);
  //    //tree->Branch("out_edep_y", &out_edep_y);
  //    //tree->Branch("out_edep_z", &out_edep_z);
  //    //tree->Branch("out_edep_blipid", &out_edep_blipid); 

  //  }
  //}

  //void set_branch_vars( bool isOutTree, int blip_idx, std::map<int, int> map_g4id_idx, float Nblips, bool eventLevel ) {   // input blip index
  //  if( isOutTree ) {
  //    if( eventLevel ) {
  //      event		 = out_event;
  //      run		 = out_run;
  //      subrun		 = out_subrun;
  //      nu_reco_vtx_x	 = out_nu_reco_vtx_x;
  //      nu_reco_vtx_y	 = out_nu_reco_vtx_y;
  //      nu_reco_vtx_z	 = out_nu_reco_vtx_z;
  //      mctruth_nu_vtx_x = out_mctruth_nu_vtx_x;
  //      mctruth_nu_vtx_y = out_mctruth_nu_vtx_y;
  //      mctruth_nu_vtx_z = out_mctruth_nu_vtx_z;
  //      mctruth_nu_pdg	 = out_mctruth_nu_pdg;
  //      mctruth_nu_ccnc	 = out_mctruth_nu_ccnc;
  //      nparticles	 = Nblips;
  //      nblips		 = Nblips;
  //    } else {
  //      int p = map_g4id_idx[blip_g4id[blip_idx]];
  //      part_process->at(p)=out_part_process->at(p);
  //      part_isPrimary[p]=out_part_isPrimary[p];
  //      part_g4id[p]=out_part_g4id[p];
  //      part_pdg[p]=out_part_pdg[p];
  //      part_mother[p]=out_part_mother[p];
  //      part_KE[p]=out_part_KE[p];
  //      part_pathlen[p]=out_part_pathlen[p];
  //      part_depEnergy[p]=out_part_depEnergy[p];
  //      part_nDaughters[p]=out_part_nDaughters[p];
  //      part_Px[p]=out_part_Px[p];
  //      part_Py[p]=out_part_Py[p];
  //      part_Pz[p]=out_part_Pz[p];
  //      part_startPointx[p]=out_part_startPointx[p];
  //      part_startPointy[p]=out_part_startPointy[p];
  //      part_startPointz[p]=out_part_startPointz[p];
  //      part_endPointx[p]=out_part_endPointx[p];
  //      part_endPointy[p]=out_part_endPointy[p];
  //      part_endPointz[p]=out_part_endPointz[p];
  //      blip_x[blip_idx]=out_blip_x[blip_idx];
  //      blip_y[blip_idx]=out_blip_y[blip_idx];
  //      blip_z[blip_idx]=out_blip_z[blip_idx];
  //      blip_size[blip_idx]=out_blip_size[blip_idx];
  //      blip_energy[blip_idx]=out_blip_energy[blip_idx];
  //      blip_charge[blip_idx]=out_blip_charge[blip_idx];
  //      blip_edepid[blip_idx]=out_blip_edepid[blip_idx];
  //      blip_dx[blip_idx]=out_blip_dx[blip_idx];
  //      blip_dw[blip_idx]=out_blip_dw[blip_idx];
  //      blip_proxtrkdist[blip_idx]=out_blip_proxtrkdist[blip_idx];
  //      blip_touchtrk[blip_idx]=out_blip_touchtrk[blip_idx];
  //      blip_pl0_bydeadwire[blip_idx]=out_blip_pl0_bydeadwire[blip_idx];
  //      blip_pl1_bydeadwire[blip_idx]=out_blip_pl1_bydeadwire[blip_idx];
  //      blip_pl2_bydeadwire[blip_idx]=out_blip_pl2_bydeadwire[blip_idx];
  //    }
  //  } else {
  //    if( eventLevel ) {
  //      out_event	     = event;
  //      out_run		     = run;
  //      out_subrun	     = subrun;
  //      out_nu_reco_vtx_x    = nu_reco_vtx_x;
  //      out_nu_reco_vtx_y    = nu_reco_vtx_y;
  //      out_nu_reco_vtx_z    = nu_reco_vtx_z;
  //      out_mctruth_nu_vtx_x = mctruth_nu_vtx_x;
  //      out_mctruth_nu_vtx_y = mctruth_nu_vtx_y;
  //      out_mctruth_nu_vtx_z = mctruth_nu_vtx_z;
  //      out_mctruth_nu_pdg   = mctruth_nu_pdg;
  //      out_mctruth_nu_ccnc  = mctruth_nu_ccnc;
  //      out_nparticles	     = Nblips;
  //      out_nblips	     = Nblips;
  //    } else {
  //      int p = map_g4id_idx[blip_g4id[blip_idx]];
  //      //out_part_process->at(p)=part_process->at(p);
  //      out_part_isPrimary[p]=part_isPrimary[p];
  //      out_part_g4id[p]=part_g4id[p];
  //      out_part_pdg[p]=part_pdg[p];
  //      out_part_nDaughters[p]=part_nDaughters[p];
  //      out_part_mother[p]=part_mother[p];
  //      out_part_KE[p]=part_KE[p];
  //      out_part_Px[p]=part_Px[p];
  //      out_part_Py[p]=part_Py[p];
  //      out_part_Pz[p]=part_Pz[p];
  //      out_part_startPointx[p]=part_startPointx[p];
  //      out_part_startPointy[p]=part_startPointy[p];
  //      out_part_startPointz[p]=part_startPointz[p];
  //      out_part_endPointx[p]=part_endPointx[p];
  //      out_part_endPointy[p]=part_endPointy[p];
  //      out_part_endPointz[p]=part_endPointz[p];
  //      out_part_pathlen[p]=part_pathlen[p];
  //      out_part_depEnergy[p]=part_depEnergy[p];
  //      
  //      //out_nedeps	[p]=nedeps;
  //      //out_edep_energy[p]=edep_energy;
  //      //out_edep_isPrimary[p]=edep_isPrimary;
  //      //out_edep_x	[p]=edep_x;
  //      //out_edep_y	[p]=edep_y;
  //      //out_edep_z	[p]=edep_z;
  //      //out_edep_blipid[p]=edep_blipid;
  //      
  //      out_blip_x[blip_idx]=blip_x[blip_idx];
  //      out_blip_y[blip_idx]=blip_y[blip_idx];
  //      out_blip_z[blip_idx]=blip_z[blip_idx];
  //      out_blip_size[blip_idx]=blip_size[blip_idx];
  //      out_blip_energy[blip_idx]=blip_energy[blip_idx];
  //      out_blip_charge[blip_idx]=blip_charge[blip_idx];
  //      out_blip_edepid[blip_idx]=blip_edepid[blip_idx];
  //      out_blip_g4id[blip_idx]=blip_g4id[blip_idx];
  //      out_blip_dx[blip_idx]=blip_dx[blip_idx];
  //      out_blip_dw[blip_idx]=blip_dw[blip_idx];
  //      out_blip_proxtrkdist[blip_idx]=blip_proxtrkdist[blip_idx];
  //      out_blip_touchtrk[blip_idx]=blip_touchtrk[blip_idx];
  //      out_blip_pl0_bydeadwire[blip_idx]=blip_pl0_bydeadwire[blip_idx];
  //      out_blip_pl1_bydeadwire[blip_idx]=blip_pl1_bydeadwire[blip_idx];
  //      out_blip_pl2_bydeadwire[blip_idx]=blip_pl2_bydeadwire[blip_idx];
  //    }
  //  }
  //}

