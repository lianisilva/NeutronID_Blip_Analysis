#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
using namespace std;

// Input file name
std::string fFileName = "blipana_neutron_0to200MeV_overlay_20250210.root";

void tree2text(){
  // Read in the TFile and grab the TTree
  TFile* file = new TFile(fFileName.c_str(),"READ");
  TTree* tree = (TTree*)file->Get("blipana/anatree");
  tr->Scan(); // prints the content on the screen

  float a,b,c; // create variables of the same type as the branches you want to access

  tr->SetBranchAddress("TS",&a); // for all the TTree branches you need this
  tr->SetBranchAddress("ns",&b);
  tr->SetBranchAddress("nserr",&c);

  ofstream myfile;
  myfile.open ("example.txt");
  myfile << "TS ns nserr\n";

  for (int i=0;i<tr->GetEntries();i++){
    // loop over the tree
    tr->GetEntry(i);
    cout << a << " " << b << " "<< c << endl; //print to the screen
    myfile << a << " " << b << " "<< c<<"\n"; //write to file
  }
  myfile.close();
}
