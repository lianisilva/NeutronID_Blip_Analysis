  #include "TCanvas.h"
  #include "TView3D.h"
  #include "TPolyLine3D.h"
  #include "TPolyMarker3D.h"
  #include "TText.h"
  #include "TLegend.h"
  #include "TPaveText.h"
  
  class Visualizer3D {
  private:
    TCanvas* canvas;
    std::vector<TPolyLine3D*> part_trks;
    std::vector<TPolyLine3D*> momentum_arrows;
    TPolyMarker3D* nu_vertex;
    TPaveText* info_box;
    
    // Fixed track length for uniform appearance
    double fixed_track_length;
    
    // Store particle end positions for chaining
    std::map<int, double> particle_end_x;
    std::map<int, double> particle_end_y;
    std::map<int, double> particle_end_z;
      
  public:
    Visualizer3D(const char* name = "3D", const char* title = "3D Particle Visualization", 
                 double track_length = 100.0) {
      canvas = new TCanvas(name, title, 1200, 800);
      canvas->SetFillColor(0);
      fixed_track_length = track_length;
      
      // Create neutrino vertex marker
      nu_vertex = new TPolyMarker3D(1, 29);   // one point, star
      nu_vertex->SetMarkerColor(kRed);
      nu_vertex->SetMarkerSize(3.0);
      
      // Info box for displaying particle information
      info_box = new TPaveText(0.03, 0.03, 0.3, 0.3, "NDC");
      info_box->SetFillColor(kWhite);
      info_box->SetBorderSize(1);
      info_box->SetTextAlign(12);
      info_box->SetTextSize(0.02);
    }
    
    ~Visualizer3D() {
        // Clean up memory
        for(auto obj : part_trks) delete obj;
        for(auto obj : momentum_arrows) delete obj;
        delete nu_vertex;
        delete info_box;
        delete canvas;
    }
    
    void SetNeutrinoVertex(double x, double y, double z) {
      nu_vertex->SetPoint(0, x, y, z);
    }
    
    // Add particle with automatic chaining based on mother g4id
    void AddParticle(int g4id, int pdg, double Px, double Py, double Pz, double energy, 
                     bool is_primary = false, int mother_g4id = -1) {
      if(nu_vertex->GetN() == 0) {
        std::cerr << "Warning: Set neutrino vertex before adding particles!" << std::endl;
        return;
      }
      
      // Determine start position
      double start_x, start_y, start_z;
      
      if(is_primary || mother_g4id == -1) {
        // Start from neutrino vertex
        nu_vertex->GetPoint(0, start_x, start_y, start_z);
      } else {
        // Start from mother particle's end position
        if(particle_end_x.find(mother_g4id) != particle_end_x.end()) {
          start_x = particle_end_x[mother_g4id];
          start_y = particle_end_y[mother_g4id];
          start_z = particle_end_z[mother_g4id];
        } else {
          // Fallback to vertex if mother not found
          nu_vertex->GetPoint(0, start_x, start_y, start_z);
        }
      }
      
      AddParticleFromPosition(g4id, pdg, start_x, start_y, start_z, Px, Py, Pz, energy, is_primary);
    }
    
    // Original method - now stores end position for chaining
    void AddParticleFromPosition(int g4id, int pdg, float start_x, float start_y, float start_z,
                                double Px, double Py, double Pz, double energy, 
                                bool is_primary = false) {
      Color_t color = GetParticleColor(pdg);
      
      // Calculate direction from momentum
      double P_mag = sqrt(Px*Px + Py*Py + Pz*Pz);
      if(P_mag == 0) return;
      
      // Normalize momentum direction
      double dir_x = Px / P_mag;
      double dir_y = Py / P_mag;
      double dir_z = Pz / P_mag;
      
      // Create uniform track from start position in momentum direction
      double end_x = start_x + dir_x * fixed_track_length;
      double end_y = start_y + dir_y * fixed_track_length;
      double end_z = start_z + dir_z * fixed_track_length;
      
      // Store end position for potential daughter particles
      particle_end_x[g4id] = end_x;
      particle_end_y[g4id] = end_y;
      particle_end_z[g4id] = end_z;
      
      // Draw particle track
      TPolyLine3D* track = new TPolyLine3D(2);
      track->SetPoint(0, start_x, start_y, start_z);
      track->SetPoint(1, end_x, end_y, end_z);
      track->SetLineColor(color);
      track->SetLineWidth(is_primary ? 4 : 2);
      part_trks.push_back(track);
      
      // Add arrowhead at end
      AddArrowHead(end_x, end_y, end_z, dir_x, dir_y, dir_z, 
                  fixed_track_length * 0.15, color);
      
      // Add particle info to info box
      AddParticleInfo(pdg, energy, P_mag, is_primary);
    }
    
    void Draw(double range = 200.0) {
      canvas->cd();
      
      if(nu_vertex->GetN() == 0) {
        std::cerr << "Error: No neutrino vertex set!" << std::endl;
        return;
      }
      
      // Get vertex position for centering
      double nu_x, nu_y, nu_z;
      nu_vertex->GetPoint(0, nu_x, nu_y, nu_z);
      
      // Set up 3D view centered on vertex
      TView* tempView = TView3D::CreateView(1);
      TView3D* view = dynamic_cast<TView3D*>(tempView);
      if (!view) { std::cerr << "Error: Failed to create TView3D" << std::endl; return; }
      
      view->SetRange(nu_x - range, nu_y - range, nu_z - range, 
                     nu_x + range, nu_y + range, nu_z + range);
      
      bool first = true;
      
      // Draw vertex and tracks
      nu_vertex->Draw(first ? "" : "same"); 
      first = false;
      
      for(auto track : part_trks) { 
        track->Draw("same"); 
      }
      
      for(auto arrow : momentum_arrows) { 
        arrow->Draw("same"); 
      }
      
      // Add legend and info
      AddLegend();
      info_box->Draw("same");
      
      canvas->Update();
    }
  
    std::unordered_map<int, Color_t> colorMap = {
      {12, kRed}, {2212, kBlue}, {2112, kViolet}, {11, kOrange}, {22, kYellow}
    };
  
    void SaveAs(const char* filename) { canvas->SaveAs(filename); }
    TCanvas* GetCanvas() { return canvas; }
  
  private:
    Color_t GetParticleColor(int pdg) {
      auto it = colorMap.find(abs(pdg));
      return (it != colorMap.end()) ? it->second : kBlack;
    }
    
    void AddArrowHead(double tip_x, double tip_y, double tip_z,
                     double dir_x, double dir_y, double dir_z,
                     double head_length, Color_t color) {
        
      // Create simple arrowhead with 4 lines forming a cone
      double perp1_x, perp1_y, perp1_z;
      double perp2_x, perp2_y, perp2_z;
      
      // Find perpendicular vectors
      if(fabs(dir_x) < 0.9) {
          perp1_x = 0; perp1_y = dir_z; perp1_z = -dir_y;
      } else {
          perp1_x = -dir_z; perp1_y = 0; perp1_z = dir_x;
      }
      
      // Normalize
      double norm = sqrt(perp1_x*perp1_x + perp1_y*perp1_y + perp1_z*perp1_z);
      if(norm > 0) {
          perp1_x /= norm; perp1_y /= norm; perp1_z /= norm;
      }
      
      // Second perpendicular (cross product)
      perp2_x = dir_y*perp1_z - dir_z*perp1_y;
      perp2_y = dir_z*perp1_x - dir_x*perp1_z;
      perp2_z = dir_x*perp1_y - dir_y*perp1_x;
      
      // Create 4 arrowhead lines
      for(int i = 0; i < 4; i++) {
          double angle = i * TMath::Pi() / 2.0;
          double px = perp1_x * cos(angle) + perp2_x * sin(angle);
          double py = perp1_y * cos(angle) + perp2_y * sin(angle);
          double pz = perp1_z * cos(angle) + perp2_z * sin(angle);
          
          TPolyLine3D* head_line = new TPolyLine3D(2);
          head_line->SetPoint(0, tip_x, tip_y, tip_z);
          head_line->SetPoint(1, tip_x - dir_x * head_length + px * head_length * 0.4,
                                 tip_y - dir_y * head_length + py * head_length * 0.4,
                                 tip_z - dir_z * head_length + pz * head_length * 0.4);
          head_line->SetLineColor(color);
          head_line->SetLineWidth(2);
          momentum_arrows.push_back(head_line);
      }
    }
    
    void AddParticleInfo(int pdg, double energy, double momentum, bool is_primary) {
      TString particle_name = GetParticleName(pdg);
      TString info = Form("%s: E=%.1f MeV, p=%.1f MeV/c", 
                         particle_name.Data(), energy, momentum);
      if(is_primary) info += " (Primary)";
      info_box->AddText(info);
  
      TText* lastText = (TText*)info_box->GetListOfLines()->Last();
      lastText->SetTextColor(GetParticleColor(pdg));
    }
    
    TString GetParticleName(int pdg) {
      switch(abs(pdg)) {
        case 2212: return "Proton";
        case 2112: return "Neutron";
        case 11:   return (pdg > 0) ? "e+" : "e-";
        case 12:   return "Neutrino";
        case 22:   return "Gamma";
        default:   return Form("PDG_%d", pdg);
      }
    }
    
    void AddLegend() {
      TLegend* legend = new TLegend(0.75, 0.75, 0.98, 0.98);
      legend->SetFillColor(kWhite);
      legend->SetBorderSize(1);
      legend->AddEntry((TObject*)0, "Particle Directions", "");
      legend->AddEntry((TObject*)0, "(Equal length tracks)", "");
      legend->AddEntry((TObject*)0, "", "");
      legend->AddEntry(nu_vertex, "Neutrino Vertex", "p");
      
      // Create dummy markers for legend
      TPolyLine3D *p_line = new TPolyLine3D(2);
      p_line->SetLineColor(kBlue);
      p_line->SetLineWidth(3);
      legend->AddEntry(p_line, "Proton", "l");
      
      TPolyLine3D *n_line = new TPolyLine3D(2);
      n_line->SetLineColor(kViolet);
      n_line->SetLineWidth(3);
      legend->AddEntry(n_line, "Neutron", "l");
      
      TPolyLine3D *e_line = new TPolyLine3D(2);
      e_line->SetLineColor(kOrange);
      e_line->SetLineWidth(3);
      legend->AddEntry(e_line, "Electron", "l");
      
      TPolyLine3D *g_line = new TPolyLine3D(2);
      g_line->SetLineColor(kYellow);
      g_line->SetLineWidth(3);
      legend->AddEntry(g_line, "Gamma", "l");
      
      legend->Draw("same");
    }
  };
