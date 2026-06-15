// Standard library includes
#include <algorithm>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"

// STV analysis includes
#include "XSecAnalyzer/FilePropertiesManager.hh"
#include "XSecAnalyzer/MCC9SystematicsCalculator.hh"
#include "XSecAnalyzer/PlotUtils.hh"
#include "XSecAnalyzer/SliceBinning.hh"
#include "XSecAnalyzer/SliceHistogram.hh"

using NFT = NtupleFileType;

void scale_by_bin_width(SliceHistogram* pSlice)
{
    int num_slice_bins = pSlice->hist_->GetNbinsX();
    TMatrixD trans_mat( num_slice_bins, num_slice_bins );
    for ( int b = 0; b < num_slice_bins; ++b ) {
      const auto width = pSlice->hist_->GetBinWidth( b + 1 );
      // width *= other_var_width;
      trans_mat( b, b ) = 1 / width;
    }
    pSlice->transform(trans_mat);
}

void set_mc_histogram_style( int event_category, TH1* mc_hist, int color ) {
    mc_hist->SetFillColor( color );
    mc_hist->SetLineColor( color );
    mc_hist->SetStats( false );
  }

void set_ext_histogram_style( TH1* ext_hist ) {
    ext_hist->SetFillColor( 28 );
    ext_hist->SetLineColor( 28 );
    ext_hist->SetLineWidth( 2 );
    ext_hist->SetFillStyle( 3005 );
    ext_hist->SetStats( false );
  }

void variable_stack_plots() {

  bool normaliseByBinWidth = false;

  std::cout << "HERE -1\n";

  auto* syst_ptr = new MCC9SystematicsCalculator(

  // LLR PID
  // "/Users/patrick/Documents/MicroBooNE/CrossSections/NCpi_Analysis/univmake_outputs/univmake-output-NCpi-HighPurity_Run4b_sidebands.root",
  // "/Users/patrick/Documents/MicroBooNE/CrossSections/NCpi_Analysis/univmake_outputs/univmake-output-NCpi-HighPurity_nuwro_fake_data.root",


  // Sideband LLR PID
  "/Users/patrick/Documents/MicroBooNE/CrossSections/NCpi_Analysis/univmake_outputs/univmake-output-CCpiSideband-Run4cd.root",


    "/Users/patrick/Documents/MicroBooNE/CrossSections/NCpi_Analysis/xsec_analyzer/configs/systcalc.conf" );

  syst_ptr->set_syst_mode(syst_ptr->SystMode::VaryBackgroundAndSignalDirectly);
  //syst_ptr->set_syst_mode(syst_ptr->SystMode::ForXSec);

  auto& syst = *syst_ptr;

  std::cout << "HERE 0\n";

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();  
  TH2D* category_hist = syst.cv_universe().hist_categ_.get();

  std::cout << "HERE 1\n";

  // Total MC+EXT prediction in reco bin space. Start by getting EXT.
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

   std::cout << "HERE 1a\n";

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

   std::cout << "HERE 1b\n";
  
  // Keys are covariance matrix types, values are CovMatrix objects that represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();

  std::cout << "HERE 1c\n";

  auto& matrix_map = *matrix_map_ptr;

   std::cout << "HERE 1d\n";

  //auto* sb_ptr = new SliceBinning( "../configs/NC1pi_slice_config.txt" );
  auto* sb_ptr = new SliceBinning( "../configs/CCpiSideband_slice_config.txt" );

  std::cout << "HERE 2\n";

  auto& sb = *sb_ptr;

  unsigned int sl_idx = 0;
  const auto& slice = sb.slices_.at( sl_idx ); // only considering single slice at a time

  // We now have all of the reco bin space histograms that we need as input.
  // Use them to make new histograms in slice space.
  SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
    *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );
  
  SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
    *reco_ext_hist, slice, &matrix_map.at("EXTstats") );
  
  SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
    *reco_mc_plus_ext_hist, slice, &matrix_map.at("PredTotal") );

  auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
  std::cout << "\u03C7\u00b2 = "
    << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bins,"
    << " p-value = " << chi2_result.p_value_ << '\n';

  // Build a stack of categorized central-value MC predictions plus the
  // extBNB contribution in slice space
  set_ext_histogram_style( slice_ext->hist_.get() );

  THStack* slice_pred_stack = new THStack( "mc+ext", "" );

  if (normaliseByBinWidth) scale_by_bin_width(slice_ext);
  slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB
  
  const auto& sel_for_cat = syst.get_selection_for_categories();
  const auto& cat_map = sel_for_cat.category_map();

  // Legend
  TLegend *leg = new TLegend(0.1,0.84,0.9,0.99);
  leg->SetLineColor(kWhite);
  leg->SetTextAlign(12);
  leg->SetNColumns(4);

  // Go in reverse so that signal ends up on top. Note that this index is
  // one-based to match the ROOT histograms
  int cat_bin_index = cat_map.size();
  
  for ( auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter )
  {
    int cat = iter->first;
    int color = iter->second.second;
    std::string label = iter->second.first;

    TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
      cat_bin_index, cat_bin_index );
    temp_mc_hist->SetDirectory( nullptr );

    SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
      *temp_mc_hist, slice  );

    set_mc_histogram_style( cat, temp_slice_mc->hist_.get(), color );

    if (normaliseByBinWidth) scale_by_bin_width(temp_slice_mc);

    slice_pred_stack->Add( temp_slice_mc->hist_.get() );

    std::string cat_col_prefix = "MC" + std::to_string( cat );

    // print number of entries for each category
    double entries = temp_slice_mc->hist_->Integral();
    std::cout << "Category " << cat << " (" << label << "): " << entries << " entries\n";

    --cat_bin_index;
  }

  // Second loop to construct legend in desired order
  cat_bin_index = 1;
  for ( auto iter = cat_map.begin(); iter != cat_map.end(); ++iter )
  {
    int cat = iter->first;
    int color = iter->second.second;
    std::string label = iter->second.first;

    if (label == "Unknown" || label == "Other") {
      ++cat_bin_index;
      continue;
    }

    TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
      cat_bin_index, cat_bin_index );
    temp_mc_hist->SetDirectory( nullptr );

    SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
      *temp_mc_hist, slice  );

    set_mc_histogram_style( cat, temp_slice_mc->hist_.get(), color );

    leg->AddEntry(temp_slice_mc->hist_.get(), label.c_str(), "f");

    ++cat_bin_index;
  }

  leg->AddEntry(slice_ext->hist_.get(), "EXT", "f");
  
  TCanvas* c1 = new TCanvas("", "", 1200, 1000);
  //TCanvas* c1 = new TCanvas("", "", 995, 1000); // ??
  TPad *upperPad = new TPad("upperPad", "Upper Pad", 0.01, 0.25, 0.99, 0.99);
  TPad *lowerPad = new TPad("lowerPad", "Lower Pad", 0.01, 0.01, 0.99, 0.24);
  upperPad->Draw();
  lowerPad->Draw();

  upperPad->cd();  // Switch to the upper pad
  gPad->SetBottomMargin(0.0125);
  gPad->SetTopMargin(0.18);

  slice_bnb->hist_->SetLineColor( kBlack );
  slice_bnb->hist_->SetLineWidth( 3 );
  slice_bnb->hist_->SetMarkerStyle( kFullCircle );
  slice_bnb->hist_->SetMarkerSize( 0.8 );
  slice_bnb->hist_->SetStats( false );

  if (normaliseByBinWidth) scale_by_bin_width(slice_bnb);
  if (normaliseByBinWidth) scale_by_bin_width(slice_mc_plus_ext);

  double ymax = std::max( slice_bnb->hist_->GetMaximum(),
    slice_mc_plus_ext->hist_->GetMaximum() ) * 1.8;
  slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., ymax );
  slice_bnb->hist_->SetTitle("");
  
  slice_bnb->hist_->GetXaxis()->SetLabelOffset(999); // Hide X-axis labels
  slice_bnb->hist_->GetXaxis()->SetTitleOffset(999); // Hide X-axis labels
  slice_bnb->hist_->GetXaxis()->SetTickLength(0.01);

  slice_bnb->hist_->GetYaxis()->SetTitleSize( 0.062 );
  slice_bnb->hist_->GetYaxis()->SetTitleOffset( 0.855 );
  slice_bnb->hist_->GetYaxis()->SetLabelSize( 0.048 );  //0.052 

  //slice_bnb->hist_->GetYaxis()->SetMaxDigits(3);
  //slice_bnb->hist_->GetYaxis()->SetExponentOffset(-0.01, 0.01);

  slice_bnb->hist_->Draw( "e" );

  slice_pred_stack->Draw( "hist same" );

  slice_mc_plus_ext->hist_->SetLineWidth( 3 ); // 3

  slice_mc_plus_ext->hist_->SetMarkerColor(kBlack);
  slice_mc_plus_ext->hist_->SetLineColor(kBlack);
  slice_mc_plus_ext->hist_->DrawCopy( "hist same" );

  slice_mc_plus_ext->hist_->SetFillColor(kBlack);
  slice_mc_plus_ext->hist_->SetFillStyle(3004);
  slice_mc_plus_ext->hist_->Draw("e2 same");

  slice_bnb->hist_->Draw( "same e" );

  // Create the label text with the value
  TLatex label;
  label.SetTextAlign(12); // Set text alignment (left-aligned)
  label.SetNDC(); // Set position in normalized coordinates
  char labelText1[100];
  char labelText2[100];
  char labelText3[100];
  char labelText4[100];
  char labelText5[100];
  //sprintf(labelText1, "NuWro Fake Data");
  sprintf(labelText1, "MicroBooNE BNB Data");
  //sprintf(labelText2, "Run 4 3.5#times10^{19} POT");
  //sprintf(labelText2, "Run 4b 1.375#times10^{20} POT");
  sprintf(labelText2, "Run 4cd 1.397#times10^{20} POT");
  //sprintf(labelText2, "Run 4d 0.493#times10^{20} POT");
  ////sprintf(labelText2, "Run 5 1.469#times10^{20} POT");
  sprintf(labelText3, "#chi^{2} = %.2f / %d Bins", chi2_result.chi2_, chi2_result.num_bins_);
  sprintf(labelText4, "p-value = %.2f", chi2_result.p_value_);
  sprintf(labelText5, "%.1f #sigma", TMath::Sqrt( TMath::ChisquareQuantile( 1 - chi2_result.p_value_, 1 ) ));
  label.SetTextSize(0.05);
  label.DrawLatex(0.55, 0.76, labelText1);
  label.DrawLatex(0.55, 0.70, labelText2);
  label.DrawLatex(0.55, 0.64, labelText3);
  label.DrawLatex(0.55, 0.58, labelText4);
  //label.DrawLatex(0.55, 0.52, labelText5);

  // draw legend
  leg->Draw("Same");

  gPad->RedrawAxis();

  std::cout << "Data events: " << slice_bnb->hist_.get()->Integral() << std::endl;
  std::cout << "EXT events: " << slice_ext->hist_.get()->Integral() << std::endl;
  std::cout << "MC+EXT events: " << slice_mc_plus_ext->hist_->Integral() << std::endl;

  // Create and draw the ratio plot
  TH1D *h_ratio = (TH1D*)slice_bnb->hist_.get()->Clone("h_ratio");
  TH1D *h_ratio_error = (TH1D*)slice_mc_plus_ext->hist_.get()->Clone("h_ratio_error");
  TH1D *h_ratio_values = (TH1D*)slice_mc_plus_ext->hist_.get()->Clone("h_ratio_values");
  
  for(int i=0 ; i <= h_ratio_values->GetNbinsX() ; i++){
    h_ratio_values->SetBinError(i, 0);
  }
  
  //h_ratio->Divide(slice_mc_plus_ext->hist_.get());
  h_ratio->Divide(h_ratio_values);

  h_ratio_error->Divide(h_ratio_values);

  lowerPad->cd();  // Switch to the lower pad
  gPad->SetBottomMargin(0.35);

  h_ratio->GetYaxis()->SetRangeUser( 0, 2 );
  h_ratio->GetYaxis()->SetTitle("Ratio");
  h_ratio->GetYaxis()->SetTitleSize(0.20); // 20
  h_ratio->GetYaxis()->SetTitleOffset(0.235);
  h_ratio->GetYaxis()->SetLabelSize(0.16);
  h_ratio->GetYaxis()->SetNdivisions(505);

  h_ratio->SetBit(TH1::kNoTitle);

  h_ratio->GetXaxis()->SetLabelSize(0.16);  // Adjust X-axis label size
  h_ratio->GetXaxis()->SetLabelOffset(0.01);
  h_ratio->GetXaxis()->SetTitleSize(0.18);
  h_ratio->GetXaxis()->SetTitleOffset(0.9);
  h_ratio->GetXaxis()->SetTickLength(0.03);

  //h_ratio->GetXaxis()->SetLabelOffset(999); // Hide X-axis labels

  h_ratio->Draw();
  h_ratio_error->Draw("e2 same");


  //std::cout << "Bin 1 error: " << h_ratio->GetBinError(1) << std::endl;


  // Get the binning and axis labels for the current slice by cloning the
  // (empty) histogram owned by the Slice object
  TH1* slice_hist = dynamic_cast< TH1* >(
    slice.hist_->Clone("slice_hist") );

  slice_hist->SetDirectory( nullptr );

  // Keys are labels, values are fractional uncertainty histograms
  auto* fr_unc_hists = new std::map< std::string, TH1* >();
  auto& frac_uncertainty_hists = *fr_unc_hists;

  // Show fractional uncertainties computed using these covariance matrices
  // in the ROOT plot. All configured fractional uncertainties will be
  // included in the output pgfplots file regardless of whether they appear
  // in this vector.

  const std::vector< std::string > cov_mat_keys = { "total",
    "detVar_total", "flux", "flux_beamline", "reint", "xsec_total", "POT", "numTargets", "dirtNorm",
    "MCstats", "EXTstats", "BNBstats", "NuWroGenie"
  };
  
  // show detvars
  /*
  const std::vector< std::string > cov_mat_keys = {"total","detVar_total",
    "detVarLYdown", "detVarLYrayl", "detVarLYatten", "detVarRecomb2", "detVarSCE",
    "detVarWMX", "detVarWMYZ"
  };
  */
  
  /*
  // show xsec uncertainties
  const std::vector< std::string > cov_mat_keys = {"total", "xsec_multi", "xsec_unisim", "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_NormCCCOH",
    "xsec_NormNCCOH", "xsec_RPA_CCQE", "xsec_ThetaDelta2NRad", "xsec_Theta_Delta2Npi", "xsec_VecFFCCQEshape", "xsec_XSecShape_CCMEC", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC" };
  */
  // Loop over the various systematic uncertainties
  int color = 0;
  for ( const auto& pair : matrix_map ) {

    const auto& key = pair.first;
    const auto& cov_matrix = pair.second;

    //std::cout << key << std::endl;

   SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &cov_matrix );

    // The SliceHistogram object already set the bin errors appropriately
    // based on the slice covariance matrix. Just change the bin contents
    // for the current histogram to be fractional uncertainties. Also set
    // the "uncertainties on the uncertainties" to zero.
    // TODO: revisit this last bit, possibly assign bin errors here
    for ( const auto& bin_pair : slice.bin_map_ ) {
      int global_bin_idx = bin_pair.first;
      double y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
      double err = slice_for_syst->hist_->GetBinError( global_bin_idx );
      double frac = 0.;
      if ( y > 0. ) frac = err / y;

      slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
      slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
    }

    // Check whether the current covariance matrix name is present in
    // the vector defined above this loop. If it isn't, don't bother to
    // plot it, and just move on to the next one.
    auto cbegin = cov_mat_keys.cbegin();
    auto cend = cov_mat_keys.cend();
    auto iter = std::find( cbegin, cend, key );
    if ( iter == cend ) continue;

    frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

    //std::cout << key << std::endl;

    if ( color <= 9 ) ++color;
    if ( color == 5 ) ++color;
    if ( color >= 10 ) color += 10;

    slice_for_syst->hist_->SetLineColor( color );
    slice_for_syst->hist_->SetLineWidth( 4 );

    if (key == "BNBstats") {
      slice_for_syst->hist_->SetLineStyle( 9 );
      slice_for_syst->hist_->SetLineWidth( 5 );
      slice_for_syst->hist_->SetLineColor( kBlack );
    }
    if (key == "NuWroGenie") slice_for_syst->hist_->SetLineColor( kOrange );
  }

  TCanvas* c2 = new TCanvas;
  TLegend* lg2 = new TLegend( 0.2, 0.7, 0.9, 0.9 );
  lg2->SetTextAlign(12);
  lg2->SetNColumns(3);

  auto* total_frac_err_hist = frac_uncertainty_hists.at( "total" );
  total_frac_err_hist->SetStats( false );
  total_frac_err_hist->GetYaxis()->SetRangeUser( 0., 0.5);
  //total_frac_err_hist->GetMaximum() * 1.05 );
  total_frac_err_hist->SetLineColor( kBlack );
  total_frac_err_hist->SetLineWidth( 5 );
  //total_frac_err_hist->SetTitle("Opening Angle");
  //total_frac_err_hist->SetTitle("Full NC 1pi- Selection");
  total_frac_err_hist->SetTitle("CC 1pi+ Sideband");
  total_frac_err_hist->GetYaxis()->SetTitle("Fractional Uncertainty");
  total_frac_err_hist->Draw( "hist" );

  stringstream ss; ss.precision(1);
  ss << std::fixed << "Total";
  //ss << ": " << total_frac_err_hist->GetBinContent( 1 )*100. << "%";
  string label_str = ss.str();

  lg2->AddEntry( total_frac_err_hist, label_str.c_str(), "l" );

  for ( auto& pair : frac_uncertainty_hists ) {
    const auto& name = pair.first;
    TH1* hist = pair.second;
    // We already plotted the "total" one above
    if ( name == "total" ) continue;

    // refine labels
    std::string name_alt = "";
    if ( name == "BNBstats") name_alt = "Data Stat"; 
    if ( name == "EXTstats") name_alt = "EXT Stat";
    if ( name == "MCstats") name_alt = "MC Stat";
    if ( name == "detVar_total") name_alt = "DetVar Total";
    if ( name == "flux") name_alt = "Flux";
    if ( name == "numTargets") name_alt = "Num Targets";
    if ( name == "xsec_total") name_alt = "X-Sec Total";
    if ( name == "reint") name_alt = "Reinteractions";

    stringstream ss; ss.precision(1);
    ss << std::fixed; 
    if (name_alt != "") ss << name_alt; 
    else ss << name;
    //ss << ": " << hist->GetBinContent( 1 )*100. << "%";
    string label_str = ss.str();
    lg2->AddEntry( hist, label_str.c_str(), "l" );
    hist->Draw( "same hist" );

    std::cout << name << " frac err in bin #1 = "
      << hist->GetBinContent( 1 )*100. << "%\n";
  }

  lg2->Draw( "same" );

  std::cout << "Total frac error in bin #1 = "
    << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";
}

int main() {
  variable_stack_plots();
  return 0;
}
