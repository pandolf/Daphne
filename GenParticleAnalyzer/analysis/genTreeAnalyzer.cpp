#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"





int main() {

  TFile* file = TFile::Open( "../genTree_QCD_Pt15to30.root" );
  TTree* tree = (TTree*)file->Get("demo/gentree");


  int nMC;
  tree->SetBranchAddress( "nMC", &nMC );
  float vertR[300];
  tree->SetBranchAddress( "vertR", vertR );
  float pMC[300];
  tree->SetBranchAddress( "pMC", pMC );
  float etaMC[300];
  tree->SetBranchAddress( "etaMC", etaMC );

  int nGoodEta=0;
  int nGoodEtaVert=0;
  int nGoodEtaVertP=0;

  int nentries = tree->GetEntries();

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 200 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

    tree->GetEntry(iEntry);


    for( unsigned i=0; i<nMC; ++i ) {

      if( fabs(etaMC[i])<2.4 ) {

        nGoodEta += 1;

        if( vertR[i] < 60. ) {

          nGoodEtaVert += 1;

          if( pMC[i] < 1.2 ) {

            nGoodEtaVertP += 1;

          } // if pt

        } // if vert

      } // if eta

    } // for MC

  } // for entries

  std::cout << "Found these kappas in " << nentries << " events: " << std::endl;
  std::cout << " |eta|<2.4   : " << nGoodEta      << " (" << 100.*nGoodEta     /nentries << "%)" << std::endl;
  std::cout << "+ vertR<60cm : " << nGoodEtaVert  << " (" << 100.*nGoodEtaVert /nentries << "%)" << std::endl;
  std::cout << "+ p<1.2 GeV  : " << nGoodEtaVertP << " (" << 100.*nGoodEtaVertP/nentries << "%)" << std::endl;

  return 0;

}
