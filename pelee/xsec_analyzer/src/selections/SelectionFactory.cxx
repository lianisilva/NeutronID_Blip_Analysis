
// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/CC1muNp0piNn.hh"
#include "XSecAnalyzer/Selections/Blip_n_Np.hh"
#include "XSecAnalyzer/Selections/DummySelection.hh"
#include "XSecAnalyzer/Selections/SelectionFactory.hh"

SelectionFactory::SelectionFactory() {
}

SelectionBase* SelectionFactory::CreateSelection(
  const std::string& selection_name )
{
  SelectionBase* sel;
  if ( selection_name == "CC1muNp0piNn" ) {
    CC1muNp0piNn* CC1muNp0piNnSel = new CC1muNp0piNn();
    sel = (SelectionBase*)CC1muNp0piNnSel;
  }
  else if ( selection_name == "Blip_n_Np" ) {
    sel = new Blip_n_Np;
  }
  else if ( selection_name == "Dummy" ) {
    sel = new DummySelection;
  }
  else {
    std::cerr << "Selection name requested: " << selection_name
      << " is not implemented in " << __FILE__ << '\n';
    throw;
  }

  // Ensure that the owned map of category definitions is set up
  sel->define_category_map();

  return sel;
}
