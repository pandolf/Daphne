#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"

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
  int nCharged[300];
  tree->SetBranchAddress( "nCharged", nCharged );
  //int decayMode[300];
  //tree->SetBranchAddress( "decayMode[nMC]", decayMode );

  float ptDau0[300];
  tree->SetBranchAddress( "ptDau0", ptDau0 );
  float etaDau0[300];
  tree->SetBranchAddress( "etaDau0", etaDau0 );
  float phiDau0[300];
  tree->SetBranchAddress( "phiDau0", phiDau0 );
  float mDau0[300];
  tree->SetBranchAddress( "mDau0", mDau0 );
  int pdgIdDau0[300];
  tree->SetBranchAddress( "pdgIdDau0", pdgIdDau0 );

  float ptDau1[300];
  tree->SetBranchAddress( "ptDau1", ptDau1 );
  float etaDau1[300];
  tree->SetBranchAddress( "etaDau1", etaDau1 );
  float phiDau1[300];
  tree->SetBranchAddress( "phiDau1", phiDau1 );
  float mDau1[300];
  tree->SetBranchAddress( "mDau1", mDau1 );
  int pdgIdDau1[300];
  tree->SetBranchAddress( "pdgIdDau1", pdgIdDau1 );

  float ptDau2[300];
  tree->SetBranchAddress( "ptDau2", ptDau2 );
  float etaDau2[300];
  tree->SetBranchAddress( "etaDau2", etaDau2 );
  float phiDau2[300];
  tree->SetBranchAddress( "phiDau2", phiDau2 );
  float mDau2[300];
  tree->SetBranchAddress( "mDau2", mDau2 );
  int pdgIdDau2[300];
  tree->SetBranchAddress( "pdgIdDau2", pdgIdDau2 );



  TH1D* h1_nTracks_nuclint = new TH1D( "nTracks_nuclint", "", 9, -0.5, 8.5 );
  TH1D* h1_nTracks_ppp     = new TH1D( "nTracks_ppp"    , "", 9, -0.5, 8.5 );
  TH1D* h1_nTracks_pee     = new TH1D( "nTracks_pee"    , "", 9, -0.5, 8.5 );

  TH1D* h1_m3_long_nuclint = new TH1D( "m3_long_nuclint", "", 50, 0., 4. );
  TH1D* h1_m3_long_ppp     = new TH1D( "m3_long_ppp"    , "", 50, 0., 4. );
  TH1D* h1_m3_long_pee     = new TH1D( "m3_long_pee"    , "", 50, 0., 4. );

  TH1D* h1_m3_ppp     = new TH1D( "m3_ppp"    , "", 100, 0.1, 0.6 );
  TH1D* h1_m3_pee     = new TH1D( "m3_pee"    , "", 100, 0.1, 0.6 );
  

  int nentries = tree->GetEntries();
  nentries = 1000000;

  std::cout << " -> Starting loop on entries" << std::endl;
 
  for( unsigned iEntry = 0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( iEntry % 100000 == 0 ) std::cout << " Entry: " << iEntry << " / " << nentries << std::endl;

    for( unsigned iMC=0; iMC<nMC; ++iMC ) {

      if( decayMode[iMC]!=0 && decayMode[iMC]!=7 && decayMode[iMC]!=8 ) continue;

      std::vector<TLorentzVector*> dau;

      float m_pi = 0.1396;
      
      TLorentzVector* dau0 = new TLorentzVector();
      dau0->SetPtEtaPhiM( ptDau0[iMC], etaDau0[iMC], phiDau0[iMC], m_pi );
      dau.push_back(dau0);

      TLorentzVector* dau1 = new TLorentzVector();
      dau1->SetPtEtaPhiM( ptDau1[iMC], etaDau1[iMC], phiDau1[iMC], m_pi );
      dau.push_back(dau1);

      TLorentzVector* dau2 = new TLorentzVector();
      dau2->SetPtEtaPhiM( ptDau2[iMC], etaDau2[iMC], phiDau2[iMC], m_pi );
      dau.push_back(dau2);

      TLorentzVector threeTracks = (*dau0) + (*dau1) + (*dau2);

      if( decayMode[iMC]==0 ) {

        h1_nTracks_nuclint->Fill( nCharged[iMC] );
        h1_m3_long_nuclint->Fill( threeTracks.M() );

      } else if( decayMode[iMC]==7 ) {

        h1_nTracks_ppp->Fill( nCharged[iMC] );
        h1_m3_long_ppp->Fill( threeTracks.M() );

      } // if decay mode

    } // for nMC

  } // for entries


  TFile* outfile = TFile::Open( Form("%s/%s/histos_temp.root", prodName.c_str(), dataset.c_str()), "RECREATE" );
  outfile->cd();

  h1_nTracks_nuclint->Write();
  h1_nTracks_ppp->Write();

  h1_m3_long_nuclint->Write();
  h1_m3_long_ppp->Write();

  outfile->Close();

  std::cout << "-> Find your histos in: " << outfile->GetName() << std::endl;

  return 0;

}
