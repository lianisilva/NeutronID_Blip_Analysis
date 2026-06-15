// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/DummySelection.hh"

DummySelection::DummySelection() : SelectionBase("DummySel") {
}

void DummySelection::define_constants() {
  //Define Reco&True FV, alongside any other constants used within selection cuts
}

void DummySelection::compute_reco_observables(AnalysisEvent* Event) {
  //Define the reconstructed kinematic variables the xsec measurement will be provided in/of interest
}

void DummySelection::compute_true_observables(AnalysisEvent* Event) {
  //Define the true kinematic variables of interest
}

int DummySelection::categorize_event(AnalysisEvent* Event) {
  //Define the event category of the selected event
  return 0;
}

bool DummySelection::define_signal(AnalysisEvent* Event) {
  //Define the MC True signal definition
  return false;
}

bool DummySelection::selection(AnalysisEvent* Event) {
  //Apply the selection cuts on Reco variables
  return false;
}

void DummySelection::define_output_branches() {
  //Save any additional variables to output TTree
}

void DummySelection::reset() {}

void DummySelection::define_category_map() {
  // Define the mapping between each integer event category and
  // a string label / color integer code pair
  std::map< int, std::pair< std::string, int > >
    temp_map = { { 0, { "Unknown", 0 } } };
}
