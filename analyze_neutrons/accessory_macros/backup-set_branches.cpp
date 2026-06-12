
  #include "in_out_variables.h"
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
    if(isOutTree) {
      SET_BRANCH(tree, out_event);
      SET_BRANCH(tree, out_run);
      SET_BRANCH(tree, out_subrun);
      SET_BRANCH(tree, out_nparticles);
      SET_BRANCH(tree, out_part_process);
      SET_BRANCH(tree, out_part_isPrimary);
      SET_BRANCH(tree, out_part_g4id);
      SET_BRANCH(tree, out_part_pdg);
      SET_BRANCH(tree, out_part_mother);
      SET_BRANCH(tree, out_part_KE);
      SET_BRANCH(tree, out_part_pathlen);
      SET_BRANCH(tree, out_part_depEnergy);
      SET_BRANCH(tree, out_part_nDaughters);
      SET_BRANCH(tree, out_part_Px);
      SET_BRANCH(tree, out_part_Py);
      SET_BRANCH(tree, out_part_Pz);
      SET_BRANCH(tree, out_part_startPointx);
      SET_BRANCH(tree, out_part_startPointy);
      SET_BRANCH(tree, out_part_startPointz);
      SET_BRANCH(tree, out_part_endPointx);
      SET_BRANCH(tree, out_part_endPointy);
      SET_BRANCH(tree, out_part_endPointz);
      SET_BRANCH(tree, out_nblips);
      SET_BRANCH(tree, out_blip_x);
      SET_BRANCH(tree, out_blip_y);
      SET_BRANCH(tree, out_blip_z);
      SET_BRANCH(tree, out_blip_size);
      SET_BRANCH(tree, out_blip_energy);
      SET_BRANCH(tree, out_blip_charge);
      SET_BRANCH(tree, out_blip_edepid);
      SET_BRANCH(tree, out_blip_g4id);
      SET_BRANCH(tree, out_blip_dx);
      SET_BRANCH(tree, out_blip_dw);
      SET_BRANCH(tree, out_blip_proxtrkdist);
      SET_BRANCH(tree, out_blip_touchtrk);
      SET_BRANCH(tree, out_blip_pl0_bydeadwire);
      SET_BRANCH(tree, out_blip_pl1_bydeadwire);
      SET_BRANCH(tree, out_blip_pl2_bydeadwire);
      SET_BRANCH(tree, out_nu_reco_vtx_x);
      SET_BRANCH(tree, out_nu_reco_vtx_y);
      SET_BRANCH(tree, out_nu_reco_vtx_z);
      SET_BRANCH(tree, out_mctruth_nu_vtx_x);
      SET_BRANCH(tree, out_mctruth_nu_vtx_y);
      SET_BRANCH(tree, out_mctruth_nu_vtx_z);
      SET_BRANCH(tree, out_mctruth_nu_pdg);
      SET_BRANCH(tree, out_mctruth_nu_ccnc);

    } else {                
      SET_BRANCH(tree, event);
      SET_BRANCH(tree, run);
      SET_BRANCH(tree, subrun);
      SET_BRANCH(tree, nparticles);
      SET_BRANCH(tree, part_process);         //tree->SetBranchAddress("part_isContained",&part_isContained);
      SET_BRANCH(tree, part_isPrimary);
      SET_BRANCH(tree, part_g4id);
      SET_BRANCH(tree, part_pdg);
      SET_BRANCH(tree, part_mother);
      SET_BRANCH(tree, part_KE);
      SET_BRANCH(tree, part_pathlen);
      SET_BRANCH(tree, part_depEnergy);
      SET_BRANCH(tree, part_nDaughters);
      SET_BRANCH(tree, part_Px);
      SET_BRANCH(tree, part_Py);
      SET_BRANCH(tree, part_Pz);
      SET_BRANCH(tree, part_startPointx);
      SET_BRANCH(tree, part_startPointy);
      SET_BRANCH(tree, part_startPointz);
      SET_BRANCH(tree, part_endPointx);
      SET_BRANCH(tree, part_endPointy);
      SET_BRANCH(tree, part_endPointz);
      //tree->SetBranchAddress("part_startT",     &part_startT);
      SET_BRANCH(tree, nedeps);
      SET_BRANCH(tree, edep_energy);
      SET_BRANCH(tree, edep_isPrimary);
      SET_BRANCH(tree, edep_x);
      SET_BRANCH(tree, edep_y);
      SET_BRANCH(tree, edep_z);
      SET_BRANCH(tree, edep_blipid);
      SET_BRANCH(tree, nblips);
      SET_BRANCH(tree, blip_x);
      SET_BRANCH(tree, blip_y);
      SET_BRANCH(tree, blip_z);
      SET_BRANCH(tree, blip_size);
      SET_BRANCH(tree, blip_energy);
      SET_BRANCH(tree, blip_charge);
      SET_BRANCH(tree, blip_edepid);
      SET_BRANCH(tree, blip_g4id);
      SET_BRANCH(tree, blip_dx);
      SET_BRANCH(tree, blip_dw);
      SET_BRANCH(tree, blip_proxtrkdist);
      SET_BRANCH(tree, blip_touchtrk);
      SET_BRANCH(tree, blip_pl0_bydeadwire);
      SET_BRANCH(tree, blip_pl1_bydeadwire);
      SET_BRANCH(tree, blip_pl2_bydeadwire);
      //tree->SetBranchAddress("ntrks",           &ntrks);
      //SET_BRANCH(tree, trk_length);
      SET_BRANCH(tree, nu_reco_vtx_x);
      SET_BRANCH(tree, nu_reco_vtx_y);
      SET_BRANCH(tree, nu_reco_vtx_z);
      SET_BRANCH(tree, mctruth_nu_vtx_x);
      SET_BRANCH(tree, mctruth_nu_vtx_y);
      SET_BRANCH(tree, mctruth_nu_vtx_z);
      SET_BRANCH(tree, mctruth_nu_pdg);
      SET_BRANCH(tree, mctruth_nu_ccnc);
      //tree->SetBranchAddress("nhits",         &nhits);
      //tree->SetBranchAddress("hit_blipid",      &hit_blipid);
      //tree->SetBranchAddress("hit_plane",       &hit_plane);
    }
  }


  // function for when you make your own
  // tree to fill, not just use the branches
  // of a premade TTree like above
  void define_branches(TTree* tree, isOutTree) {
    if( isOutTree ) {
      tree->Branch("event", &event);
      tree->Branch("run", &run);
      tree->Branch("subrun", &subrun);
      tree->Branch("nparticles", &nparticles);
      tree->Branch("part_process", &part_process);
      tree->Branch("part_isPrimary", &part_isPrimary);
      tree->Branch("part_g4id", &part_g4id);
      tree->Branch("part_pdg", &part_pdg);
      tree->Branch("part_mother", &part_mother);
      tree->Branch("part_KE", &part_KE);
      tree->Branch("part_pathlen", &part_pathlen);
      tree->Branch("part_depEnergy", &part_depEnergy); 
      tree->Branch("part_nDaughters", &part_nDaughters);
      tree->Branch("part_Px", &part_Px);
      tree->Branch("part_Py", &part_Py);
      tree->Branch("part_Pz", &part_Pz);
      tree->Branch("part_startPointx", &part_startPointx);
      tree->Branch("part_startPointy", &part_startPointy);
      tree->Branch("part_startPointz", &part_startPointz);
      tree->Branch("part_endPointx", &part_endPointx);
      tree->Branch("part_endPointy", &part_endPointy);
      tree->Branch("part_endPointz", &part_endPointz);
      tree->Branch("nblips", &nblips);         
      tree->Branch("blip_x", &blip_x);
      tree->Branch("blip_y", &blip_y);
      tree->Branch("blip_z", &blip_z);
      tree->Branch("blip_size", &blip_size);
      tree->Branch("blip_energy", &blip_energy);
      tree->Branch("blip_charge", &blip_charge);
      tree->Branch("blip_edepid", &blip_edepid);
      tree->Branch("blip_g4id", &blip_g4id);
      tree->Branch("blip_dx", &blip_dx);
      tree->Branch("blip_dw", &blip_dw);
      tree->Branch("blip_proxtrkdist", &blip_proxtrkdist);
      tree->Branch("blip_touchtrk", &blip_touchtrk);
      tree->Branch("blip_pl0_bydeadwire", &blip_pl0_bydeadwire);
      tree->Branch("blip_pl1_bydeadwire", &blip_pl1_bydeadwire);
      tree->Branch("blip_pl2_bydeadwire", &blip_pl2_bydeadwire);
      tree->Branch("nu_reco_vtx_x", &nu_reco_vtx_x);
      tree->Branch("nu_reco_vtx_y", &nu_reco_vtx_y);
      tree->Branch("nu_reco_vtx_z", &nu_reco_vtx_z);
      tree->Branch("mctruth_nu_vtx_x", &mctruth_nu_vtx_x);
      tree->Branch("mctruth_nu_vtx_y", &mctruth_nu_vtx_y);
      tree->Branch("mctruth_nu_vtx_z", &mctruth_nu_vtx_z);
      tree->Branch("mctruth_nu_pdg", &mctruth_nu_pdg);
      tree->Branch("mctruth_nu_ccnc", &mctruth_nu_ccnc);
    } else {
      tree->Branch("out_event", &out_event);
      tree->Branch("out_run", &out_run);
      tree->Branch("out_subrun", &out_subrun);
      tree->Branch("out_nparticles", &out_nparticles);
      tree->Branch("out_part_process", &out_part_process);
      tree->Branch("out_part_isPrimary", &out_part_isPrimary);
      tree->Branch("out_part_g4id", &out_part_g4id);
      tree->Branch("out_part_pdg", &out_part_pdg);
      tree->Branch("out_part_mother", &out_part_mother);
      tree->Branch("out_part_KE", &out_part_KE);
      tree->Branch("out_part_pathlen", &out_part_pathlen);
      tree->Branch("out_part_depEnergy", &out_part_depEnergy); 
      tree->Branch("out_part_nDaughters", &out_part_nDaughters);
      tree->Branch("out_part_Px", &out_part_Px);
      tree->Branch("out_part_Py", &out_part_Py);
      tree->Branch("out_part_Pz", &out_part_Pz);
      tree->Branch("out_part_startPointx", &out_part_startPointx);
      tree->Branch("out_part_startPointy", &out_part_startPointy);
      tree->Branch("out_part_startPointz", &out_part_startPointz);
      tree->Branch("out_part_endPointx", &out_part_endPointx);
      tree->Branch("out_part_endPointy", &out_part_endPointy);
      tree->Branch("out_part_endPointz", &out_part_endPointz);
      //tree->Branch("out_nedeps", &out_nedeps);
      //tree->Branch("out_edep_energy", &out_edep_energy);
      //tree->Branch("out_edep_isPrimary", &out_edep_isPrimary);
      //tree->Branch("out_edep_x", &out_edep_x);
      //tree->Branch("out_edep_y", &out_edep_y);
      //tree->Branch("out_edep_z", &out_edep_z);
      //tree->Branch("out_edep_blipid", &out_edep_blipid); 
      tree->Branch("out_nblips", &out_nblips);         
      tree->Branch("out_blip_x", &out_blip_x);
      tree->Branch("out_blip_y", &out_blip_y);
      tree->Branch("out_blip_z", &out_blip_z);
      tree->Branch("out_blip_size", &out_blip_size);
      tree->Branch("out_blip_energy", &out_blip_energy);
      tree->Branch("out_blip_charge", &out_blip_charge);
      tree->Branch("out_blip_edepid", &out_blip_edepid);
      tree->Branch("out_blip_g4id", &out_blip_g4id);
      tree->Branch("out_blip_dx", &out_blip_dx);
      tree->Branch("out_blip_dw", &out_blip_dw);
      tree->Branch("out_blip_proxtrkdist", &out_blip_proxtrkdist);
      tree->Branch("out_blip_touchtrk", &out_blip_touchtrk);
      tree->Branch("out_blip_pl0_bydeadwire", &out_blip_pl0_bydeadwire);
      tree->Branch("out_blip_pl1_bydeadwire", &out_blip_pl1_bydeadwire);
      tree->Branch("out_blip_pl2_bydeadwire", &out_blip_pl2_bydeadwire);
      tree->Branch("out_nu_reco_vtx_x", &out_nu_reco_vtx_x);
      tree->Branch("out_nu_reco_vtx_y", &out_nu_reco_vtx_y);
      tree->Branch("out_nu_reco_vtx_z", &out_nu_reco_vtx_z);
      tree->Branch("out_mctruth_nu_vtx_x", &out_mctruth_nu_vtx_x);
      tree->Branch("out_mctruth_nu_vtx_y", &out_mctruth_nu_vtx_y);
      tree->Branch("out_mctruth_nu_vtx_z", &out_mctruth_nu_vtx_z);
      tree->Branch("out_mctruth_nu_pdg", &out_mctruth_nu_pdg);
      tree->Branch("out_mctruth_nu_ccnc", &out_mctruth_nu_ccnc);
    }
  }

  void set_branch_vars( bool isOutTree, int blip_idx, std::map<int, int> map_g4id_idx, bool eventLevel ) {   // input blip index
    if( isOutTree ) {
      if( eventLevel ) {
	event		 = out_event;
	run		 = out_run;
	subrun		 = out_subrun;
	nu_reco_vtx_x	 = out_nu_reco_vtx_x;
	nu_reco_vtx_y	 = out_nu_reco_vtx_y;
	nu_reco_vtx_z	 = out_nu_reco_vtx_z;
	mctruth_nu_vtx_x = out_mctruth_nu_vtx_x;
	mctruth_nu_vtx_y = out_mctruth_nu_vtx_y;
	mctruth_nu_vtx_z = out_mctruth_nu_vtx_z;
	mctruth_nu_pdg	 = out_mctruth_nu_pdg;
	mctruth_nu_ccnc	 = out_mctruth_nu_ccnc;
	nparticles	 = part_g4id.size();
	nblips		 = blip_g4id.size();
      } else {
	part_process.push_back(out_part_process);
	part_isPrimary.push_back(out_part_isPrimary);
	part_g4id.push_back(out_part_g4id);
	part_pdg.push_back(out_part_pdg);
	part_mother.push_back(out_part_mother);
	part_KE.push_back(out_part_KE);
	part_pathlen.push_back(out_part_pathlen);
	part_depEnergy.push_back(out_part_depEnergy);
	part_nDaughters.push_back(out_part_nDaughter);
	part_Px.push_back(out_part_Px);
	part_Py.push_back(out_part_Py);
	part_Pz.push_back(out_part_Pz);
	part_startPointx.push_back(out_part_startPointx);
	part_startPointy.push_back(out_part_startPointy);
	part_startPointz.push_back(out_part_startPointz);
	part_endPointx.push_back(out_part_endPointx);
	part_endPointy.push_back(out_part_endPointy);
	part_endPointz.push_back(out_part_endPointz);
	blip_x.push_back(out_blip_x);
	blip_y.push_back(out_blip_y);
	blip_z.push_back(out_blip_z);
	blip_size.push_back(out_blip_size);
	blip_energy.push_back(out_blip_energy);
	blip_charge.push_back(out_blip_charge);
	blip_edepid.push_back(out_blip_edepid);
	blip_dx.push_back(out_blip_dx);
	blip_dw.push_back(out_blip_dw);
	blip_proxtrkdist.push_back(out_blip_proxtrkdist);
	blip_touchtrk.push_back(out_blip_touchtrk);
	blip_pl0_bydeadwire.push_back(out_blip_pl0_bydeadwire);
	blip_pl1_bydeadwire.push_back(out_blip_pl1_bydeadwire);
	blip_pl2_bydeadwire.push_back(out_blip_pl2_bydeadwire);
      }
    } else {
      if( eventLevel ) {
	out_event	     = event;
	out_run		     = run;
	out_subrun	     = subrun;
	out_nu_reco_vtx_x    = nu_reco_vtx_x;
	out_nu_reco_vtx_y    = nu_reco_vtx_y;
	out_nu_reco_vtx_z    = nu_reco_vtx_z;
	out_mctruth_nu_vtx_x = mctruth_nu_vtx_x;
	out_mctruth_nu_vtx_y = mctruth_nu_vtx_y;
	out_mctruth_nu_vtx_z = mctruth_nu_vtx_z;
	out_mctruth_nu_pdg   = mctruth_nu_pdg;
	out_mctruth_nu_ccnc  = mctruth_nu_ccnc;
	out_nparticles	     = out_part_g4id.size(); //nparticles;
	out_nblips	     = out_blip_g4id.size();
      } else {
	int p = map_g4id_idx[blip_g4id[blip_idx]];
	out_part_process.push_back(part_process->at(p));
	out_part_isPrimary.push_back(part_isPrimary[p]);
	out_part_g4id.push_back(part_g4id[p]);
	out_part_pdg.push_back(part_pdg[p]);
	out_part_nDaughters.push_back(part_nDaughters[p]);
	out_part_mother.push_back(part_mother[p]);
	out_part_KE.push_back(part_KE[p]);
	out_part_Px.push_back(part_Px[p]);
	out_part_Py.push_back(part_Py[p]);
	out_part_Pz.push_back(part_Pz[p]);
	out_part_startPointx.push_back(part_startPointx[p]);
	out_part_startPointy.push_back(part_startPointy[p]);
	out_part_startPointz.push_back(part_startPointz[p]);
	out_part_endPointx.push_back(part_endPointx[p]);
	out_part_endPointy.push_back(part_endPointy[p]);
	out_part_endPointz.push_back(part_endPointz[p]);
	out_part_pathlen.push_back(part_pathlen[p]);
	out_part_depEnergy.push_back(part_depEnergy[p]);
	
	//out_nedeps	.push_back(nedeps;
	//out_edep_energy.push_back(edep_energy;
	//out_edep_isPrimary.push_back(edep_isPrimary;
	//out_edep_x	.push_back(edep_x;
	//out_edep_y	.push_back(edep_y;
	//out_edep_z	.push_back(edep_z;
	//out_edep_blipid.push_back(edep_blipid;
	
	out_blip_x.push_back(blip_x[blip_idx]);
	out_blip_y.push_back(blip_y[blip_idx]);
	out_blip_z.push_back(blip_z[blip_idx]);
	out_blip_size.push_back(blip_size[blip_idx]);
	out_blip_energy.push_back(blip_energy[blip_idx]);
	out_blip_charge.push_back(blip_charge[blip_idx]);
	out_blip_edepid.push_back(blip_edepid[blip_idx]);
	out_blip_g4id.push_back(blip_g4id[blip_idx]);
	out_blip_dx.push_back(blip_dx[blip_idx]);
	out_blip_dw.push_back(blip_dw[blip_idx]);
	out_blip_proxtrkdist.push_back(blip_proxtrkdist[blip_idx]);
	out_blip_touchtrk.push_back(blip_touchtrk[blip_idx]);
	out_blip_pl0_bydeadwire.push_back(blip_pl0_bydeadwire[blip_idx]);
	out_blip_pl1_bydeadwire.push_back(blip_pl1_bydeadwire[blip_idx]);
	out_blip_pl2_bydeadwire.push_back(blip_pl2_bydeadwire[blip_idx]);
      }
    }
  }

