#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"

bool TEST = false;



int main( int argc, char* argv[]) {

  TFile* file;
  TTree* tree;

  if( argc < 2 || TEST ) {

    file = TFile::Open( "../genTree_QCD_Pt15to30.root" );
    tree = (TTree*)file->Get("demo/gentree");

    if( tree ) {
      std::cout << "-> You didn't pass enough arguments, so running on test tree: ../genTree_QCD_Pt15to30.root" << std::endl;
      TEST = true;
    } else {
      std::cout << "-> USAGE: ./drawGenTree [prodName] [dataset=\'QCD_Pt_15to30\']" << std::endl;
      exit(1);
    }

  }


  std::string prodName(argv[1]);
  std::string dataset("QCD_Pt_15to30");
  if( argc>2 ) {
    dataset = std::string(argv[2]);
  }

  file = TFile::Open( Form("%s/%s/mergedTree.root", prodName.c_str(), dataset.c_str()) );
  std::cout << " -> Opened file: " << file->GetName() << std::endl;
  tree = (TTree*)file->Get("gentree");


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

  float m_ppp[300];
  tree->SetBranchAddress( "m_ppp", m_ppp );
  float m_pee0[300];
  tree->SetBranchAddress( "m_pee0", m_pee0 );
  float m_pee1[300];
  tree->SetBranchAddress( "m_pee1", m_pee1 );
  float m_pee2[300];
  tree->SetBranchAddress( "m_pee2", m_pee2 );


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


  TFile* outfile = TFile::Open( "histos.root", "RECREATE" );
  outfile->cd();

  TH1D* h1_cutflow = new TH1D( "cutflow", "", 6, -0.5, 5.5 );
  TH1D* h1_nCharged_nuclint = new TH1D( "nCharged_nuclint", "", 21, -0.5, 20.5 );

  TH1D* h1_mPPP_d0       = new TH1D( "mPPP_d0"       , "", 50., 0., 1.);
  TH1D* h1_mPEE_d0       = new TH1D( "mPEE_d0"       , "", 50., 0., 1.);

  TH1D* h1_mPPP_d7       = new TH1D( "mPPP_d7"       , "", 50., 0., 1.);
  TH1D* h1_mPEE_d7       = new TH1D( "mPEE_d7"       , "", 50., 0., 1.);

  TH1D* h1_mPPP_d8       = new TH1D( "mPPP_d8"       , "", 50., 0., 1.);
  TH1D* h1_mPEE_d8       = new TH1D( "mPEE_d8"       , "", 50., 0., 1.);
  TH1D* h1_mPEE_d8_right = new TH1D( "mPEE_d8_right" , "", 50., 0., 1.);
  TH1D* h1_mPEE_d8_wrong = new TH1D( "mPEE_d8_wrong" , "", 50., 0., 1.);


  int nGoodEta=0;
  int nGoodEtaVert=0;
  int nGoodEtaVertP=0;
  int nGoodEtaVertP_pee=0;
  int nGoodEtaVertP_peeLowP=0;

  int nentries = tree->GetEntries();
nentries = 100000;

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

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


            //if( nCharged[i]==3 ) {
            if( (decayMode[i]==0 && (nCharged[i]==3 || nCharged[i]==4)) || decayMode[i]==7 || decayMode[i]==8 ) {

              h1_cutflow->Fill( 4 );
              nGoodEtaVertP_pee += 1;


              if( decayMode[i]==0 ) { // nuclear interactions

                h1_mPPP_d8->Fill( m_ppp [i] );

                h1_mPEE_d8->Fill( m_pee0[i] );
                h1_mPEE_d8->Fill( m_pee1[i] );
                h1_mPEE_d8->Fill( m_pee2[i] );

              } else if( decayMode[i]==7 ) { // pi+ pi- pi+

                h1_mPPP_d7->Fill( m_ppp[i] );

                h1_mPEE_d7->Fill( m_pee0[i] );
                h1_mPEE_d7->Fill( m_pee1[i] );
                h1_mPEE_d7->Fill( m_pee2[i] );

              } else if( decayMode[i]==8 ) { // pi+ e- e+

                h1_mPPP_d8->Fill( m_ppp [i] );

                h1_mPEE_d8->Fill( m_pee0[i] );
                h1_mPEE_d8->Fill( m_pee1[i] );
                h1_mPEE_d8->Fill( m_pee2[i] );

                if(        abs(pdgIdDau0[i])==11  && abs(pdgIdDau1[i])==11  && abs(pdgIdDau2[i])==211 ) {

                  h1_mPEE_d8_wrong->Fill( m_pee0[i] );
                  h1_mPEE_d8_wrong->Fill( m_pee1[i] );
                  h1_mPEE_d8_right->Fill( m_pee2[i] );

                } else if( abs(pdgIdDau0[i])==11  && abs(pdgIdDau1[i])==211 && abs(pdgIdDau2[i])==11  ) {

                  h1_mPEE_d8_wrong->Fill( m_pee0[i] );
                  h1_mPEE_d8_right->Fill( m_pee1[i] );
                  h1_mPEE_d8_wrong->Fill( m_pee2[i] );

                } else if( abs(pdgIdDau0[i])==211 && abs(pdgIdDau1[i])==11  && abs(pdgIdDau2[i])==11  ) {

                  h1_mPEE_d8_right->Fill( m_pee0[i] );
                  h1_mPEE_d8_wrong->Fill( m_pee1[i] );
                  h1_mPEE_d8_wrong->Fill( m_pee2[i] );

                } else { // shouldnt be possible

                  std::cout << "THERE IS AN ERROR THIS SHOULDN'T BE POSSIBLE" << std::endl;
                  h1_mPEE_d8_wrong->Fill( m_pee0[i] );
                  h1_mPEE_d8_wrong->Fill( m_pee1[i] );
                  h1_mPEE_d8_wrong->Fill( m_pee2[i] );

                }

              } // if decaymode


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

  h1_mPPP_d0->Write();
  h1_mPEE_d0->Write();
  h1_mPPP_d7->Write();
  h1_mPEE_d7->Write();
  h1_mPPP_d8->Write();
  h1_mPEE_d8->Write();
  h1_mPEE_d8_right->Write();
  h1_mPEE_d8_wrong->Write();

  h1_nCharged_nuclint->Write();

  outfile->Close();

  std::cout << "-> Find your histos in: " << outfile->GetName() << std::endl;

  return 0;

}
