#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "../interface/DaphneCommon.h"




int main( int argc, char* argv[]) {


  if( argc < 2 ) {
    std::cout << "-> USAGE: ./drawGenHistos [prodName] [dataset=\'QCD_Pt_15to30\']" << std::endl;
    exit(1);
  }

  DaphneCommon::setStyle();

  TPaveText* labelTop = DaphneCommon::getLabelTopSimulation();

  std::string prodName(argv[1]);
  std::string dataset("QCD_Pt_15to30");
  if( argc>2 ) {
    dataset = std::string(argv[2]);
  }


  std::string outdir( Form("%s/%s/plots", prodName.c_str(), dataset.c_str()) );
  system( Form("mkdir -p %s", outdir.c_str()) );

  TFile* file = TFile::Open( Form("%s/%s/histos_FAST.root", prodName.c_str(), dataset.c_str()) );
  //TFile* file = TFile::Open( Form("%s/%s/histos.root", prodName.c_str(), dataset.c_str()) );
  std::cout << " -> Opened file: " << file->GetName() << std::endl;


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH1D* h1_nCharged = (TH1D*)file->Get("nCharged_nuclint");

  h1_nCharged->SetXTitle( "Number of Vertex Tracks" );
  h1_nCharged->SetYTitle( "Normalized to Unity" );

  h1_nCharged->SetFillStyle(3004);
  h1_nCharged->SetFillColor(46);
  h1_nCharged->SetLineColor(46);
  h1_nCharged->SetLineWidth(2);
  h1_nCharged->DrawNormalized();

  labelTop->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/nCharged_nuclint.eps", outdir.c_str()) );
  c1->SaveAs( Form("%s/nCharged_nuclint.pdf", outdir.c_str()) );

  c1->Clear();


  TH1D* h1_mPPP_long_d0 = (TH1D*)file->Get("mPPP_long_d0");
  TH1D* h1_mPPP_long_d7 = (TH1D*)file->Get("mPPP_long_d7");
  
  h1_mPPP_long_d7->SetXTitle("M(3 tracks) [GeV]");
  h1_mPPP_long_d7->SetYTitle("Normalized to Unity");

  h1_mPPP_long_d7->SetFillColor(38);
  h1_mPPP_long_d7->SetLineColor(38);
  h1_mPPP_long_d7->SetLineWidth(2);
  h1_mPPP_long_d7->SetFillStyle(3004);

  h1_mPPP_long_d0->SetFillColor(46);
  h1_mPPP_long_d0->SetLineColor(46);
  h1_mPPP_long_d0->SetLineWidth(2);
  h1_mPPP_long_d0->SetFillStyle(3005);

  h1_mPPP_long_d7->DrawNormalized();
  h1_mPPP_long_d0->DrawNormalized("same");

  TLegend* legend = new TLegend( 0.6, 0.7, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_mPPP_long_d7, "K #rightarrow #pi #pi #pi", "F" );
  legend->AddEntry( h1_mPPP_long_d0, "Nucl. Int.", "F" );
  legend->Draw("same");

  labelTop->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/mPPP_d07.eps", outdir.c_str()) );
  c1->SaveAs( Form("%s/mPPP_d07.pdf", outdir.c_str()) );


  return 0;

}
