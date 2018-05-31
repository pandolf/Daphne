#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include "../interface/DaphneCommon.h"



int main( int argc, char* argv[]) {


  if( argc < 2 ) {
    std::cout << "-> USAGE: ./drawGenTree [prodName] [dataset=\'QCD_Pt_15to30\']" << std::endl;
    exit(1);
  }


  std::string prodName(argv[1]);
  std::string dataset("QCD_Pt_15to30");
  if( argc>2 ) {
    dataset = std::string(argv[2]);
  }

  TFile* file = TFile::Open( Form("%s/%s/mergedTree.root", prodName.c_str(), dataset.c_str()) );
  std::cout << " -> Opened file: " << file->GetName() << std::endl;
  TTree* tree = (TTree*)file->Get("gentree");

  int nMC;
  tree->SetBranchAddress( "nMC", &nMC );
  float vertR[300];
  tree->SetBranchAddress( "vertR", vertR );
  float pMC[300];
  tree->SetBranchAddress( "pMC", pMC );
  float etaMC[300];
  tree->SetBranchAddress( "etaMC", etaMC );
  int decayMode[300];
  tree->SetBranchAddress( "decayMode", decayMode );
  //int decayMode[300];
  //tree->SetBranchAddress( "decayMode[nMC]", decayMode );

  TH1D* h1_m3_nuclint = new TH1D( "m3_nuclint", "", 100, 0.1, 0.6 );
  TH1D* h1_m3_ppp     = new TH1D( "m3_ppp"    , "", 100, 0.1, 0.6 );
  TH1D* h1_m3_pee     = new TH1D( "m3_pee"    , "", 100, 0.1, 0.6 );
  

  int nentries = tree->GetEntries();

  std::cout << " -> Starting loop on entries" << std::endl;
 
  for( unsigned iEntry = 0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( iEntry % 100000 == 0 ) std::cout << " Entry: " << iEntry << " / " << nentries << std::endl;

    for( unsigned iMC=0; iMC<nMC; ++iMC ) {

      if( decayMode[iMC]!=0 && decayMode[iMC]!=7 && decayMode[iMC]!=8 ) continue;

    } // for nMC

  } // for entries

  return 0;

}
