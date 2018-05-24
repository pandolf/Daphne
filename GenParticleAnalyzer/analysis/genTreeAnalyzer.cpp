#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"





int main() {

  //TFile* file = TFile::Open( "../genTree_Flat15_3000.root" );
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
  int decayMode[300];
  tree->SetBranchAddress( "decayMode", decayMode );
  int nLowP[300];
  tree->SetBranchAddress( "nLowP", nLowP );
  int nCharged[300];
  tree->SetBranchAddress( "nCharged", nCharged );
  int nDau[300];
  tree->SetBranchAddress( "nDau", nDau );
  float ptDau[300][30];
  tree->SetBranchAddress( "ptDau", ptDau );
  float etaDau[300][30];
  tree->SetBranchAddress( "etaDau", etaDau );
  float phiDau[300][30];
  tree->SetBranchAddress( "phiDau", phiDau );
  float mDau[300][30];
  tree->SetBranchAddress( "mDau", mDau );
  int pdgIdDau[300][30];
  tree->SetBranchAddress( "pdgIdDau", pdgIdDau );


  TFile* outfile = TFile::Open( "histos.root", "RECREATE" );
  outfile->cd();

  TH1D* h1_cutflow = new TH1D( "cutflow", "", 6, -0.5, 5.5 );
  TH1D* h1_nCharged_nuclint = new TH1D( "nCharged_nuclint", "", 21, -0.5, 20.5 );
  TH1D* h1_mass_ppp = new TH1D( "mass_ppp", "", 50., 0., 1.);
  TH1D* h1_mass_pee = new TH1D( "mass_pee", "", 50., 0., 1.);
  TH1D* h1_mass_nuclint = new TH1D( "mass_nuclint", "", 50., 0., 1.);

  int nGoodEta=0;
  int nGoodEtaVert=0;
  int nGoodEtaVertP=0;
  int nGoodEtaVertP_pee=0;
  int nGoodEtaVertP_peeLowP=0;

  int nentries = tree->GetEntries();

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 200 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

    tree->GetEntry(iEntry);


    for( unsigned i=0; i<nMC; ++i ) {

      h1_cutflow->Fill( 0 );

      if( fabs(etaMC[i])<2.4 ) {

        h1_cutflow->Fill( 1 );
        nGoodEta += 1;

        if( vertR[i] < 60. && vertR[i]>0. ) {

          h1_cutflow->Fill( 2 );
          nGoodEtaVert += 1;

          if( pMC[i] < 1.05 ) {

            h1_cutflow->Fill( 3 );
            nGoodEtaVertP += 1;

            if( decayMode[i] == 0 ) { // nuclear interactions
              h1_nCharged_nuclint->Fill( nCharged[i] );
            }


            if( nCharged[i]==3 ) {
            //if( decayMode[i]==8 || decayMode[i]==7 ) {

              h1_cutflow->Fill( 4 );
              nGoodEtaVertP_pee += 1;

              std::vector<TLorentzVector> daughters;
              TLorentzVector tot(0.,0.,0.,0.);
std::cout << std::endl;

              for( unsigned iD=0; iD<nDau[i]; ++iD ) {

                TLorentzVector thisDau;
                thisDau.SetPtEtaPhiM( ptDau[i][iD], etaDau[i][iD], phiDau[i][iD], mDau[i][iD] );
std::cout << "pdgid: " << pdgId[i][iD] << " pt: " << ptDau[i][iD] << " eta: " << etaDau[i][iD] << " phi: " << phiDau[i][iD] << " m: " << mDau[i][iD] << std::endl;
                tot += thisDau;
                daughters.push_back(thisDau);

              } // for daughters

              float mass = tot.M();

              if( decayMode[i]==7 ) { // pi+ pi- pi+

                h1_mass_ppp->Fill( mass );

              } else if( decayMode[i]==8 ) { // pi+ e+ e-

                h1_mass_pee->Fill( mass );

              } else if( decayMode[i]==0 ) { // nuclear interactions

                h1_mass_nuclint->Fill( mass );

              }



              if( nLowP[i]>=2 ) {

                h1_cutflow->Fill( 5 );
                nGoodEtaVertP_peeLowP += 1;

              } // if low P

            } // if pee

          } // if p

        } // if vert

      } // if eta

    } // for MC

  } // for entries


  h1_cutflow->Scale(1./nentries);


  std::cout << "Found these kappas in " << nentries << " events: " << std::endl;
  std::cout << " |eta|<2.4   : " << nGoodEta        << " (" << (float)nGoodEta       /nentries << " per event)" << std::endl;
  std::cout << "+ vertR<60cm : " << nGoodEtaVert    << " (" << (float)nGoodEtaVert   /nentries << " per event)" << std::endl;
  std::cout << "+ p<1.05 GeV : " << nGoodEtaVertP   << " (" << (float)nGoodEtaVertP  /nentries << " per event)" << std::endl;
  std::cout << "+ decays to eep  : " << nGoodEtaVertP_pee << " (" << (float)nGoodEtaVertP_pee/nentries << " per event)" << std::endl;
  std::cout << "+ 2 daughters below 0.15  : " << nGoodEtaVertP_peeLowP << " (" << (float)nGoodEtaVertP_peeLowP/nentries << " per event)" << std::endl;

  outfile->cd();

  h1_cutflow->Write();
  h1_mass_ppp->Write();
  h1_mass_pee->Write();
  h1_mass_nuclint->Write();
  h1_nCharged_nuclint->Write();

  outfile->Close();

  std::cout << "-> Find your histos in: " << outfile->GetName() << std::endl;

  return 0;

}
