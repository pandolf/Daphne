#include "../interface/DaphneCommon.h"

#include "TCanvas.h"
#include "TH2D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMath.h"







int main() {

  DaphneCommon::setStyle();

  float xMax = 370.;

  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax, 10, 1.5E-7, 1. );
  h2_axes->SetXTitle( "m_{A'} [MeV]" );
  h2_axes->SetYTitle( "BR( X ) / #epsilon ^{2}" );

  TF1* br_pizero = new TF1("br_pizero", "TMath::Max( 0., 0.21*TMath::Power( 1.-x*x/(135.*135.) , 3. ))", 0., xMax);
  br_pizero->SetLineColor(46); 
  br_pizero->SetLineWidth(2); 

  // the following is taken from eq. 12 of http://arxiv.org/abs/0811.1030v1
  TF1* br_kappa3 = new TF1("br_kappa", "8E-5*x*x/10000.", 0., xMax);
  br_kappa3->SetLineColor(38);
  br_kappa3->SetLineWidth(2);

  // the following is taken from eq. 14 of https://arxiv.org/pdf/1309.5084.pdf
  TF1* br_kappa3_2 = new TF1("br_kappa_2", "x*x/(1024.*137.*97.4*494.*5.3E-2)*(3.+6.*x*x/(494.*494.))*TMath::Max( TMath::Power(1.+(140.*140.-x*x)*(140.*140.-x*x)/(494.*494.*494.*494.) - 2.*(140.*140. + x*x)/(494.*494.), 1.5), 0.)", 1., 400.);
  br_kappa3_2->SetLineColor(38);
  //br_kappa3_2->SetLineColor(kRed+3);
  br_kappa3_2->SetLineWidth(2);
  //br_kappa3_2->SetLineStyle(2);

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  c1->SetLogy();

  h2_axes->Draw();

  br_pizero->Draw("L same" );
  //br_kappa3->Draw("L same" );
  br_kappa3_2->Draw("L same" );

  TLegend* legend = new TLegend( 0.5, 0.75, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( br_pizero, "K^{#pm} #rightarrow #pi^{#pm}#pi^{0} #rightarrow #pi^{#pm}A'#gamma", "L" );
  legend->AddEntry( br_kappa3, "K^{#pm} #rightarrow #pi^{#pm} A'", "L" );
  legend->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs("branchingRatios.pdf");
  c1->SaveAs("branchingRatios.eps");

  c1->Clear();

  TF1* br_delta = new TF1("br_delta", "0.005*TMath::Power(1.- x*x/((1232.-938.)*(1232.-938.)), 1.5 )", 0., xMax);
  br_delta->SetLineColor(30);
  br_delta->SetLineWidth(2);

  h2_axes->Draw();
  br_pizero->Draw("L same" );
  br_kappa3->Draw("L same" );
  br_delta->Draw("L same");

  TLegend* legend2 = new TLegend( 0.5, 0.72, 0.9, 0.9 );
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.035);
  legend2->AddEntry( br_pizero, "K^{#pm} #rightarrow #pi^{#pm}#pi^{0} #rightarrow #pi^{#pm}A'#gamma", "L" );
  legend2->AddEntry( br_kappa3, "K^{#pm} #rightarrow #pi^{#pm} A'", "L" );
  legend2->AddEntry( br_delta , "#Delta #rightarrow p A'", "L" );
  legend2->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs("branchingRatios_all3.pdf");
  c1->SaveAs("branchingRatios_all3.eps");


  std::cout << "BR(m=200MeV) = " << br_kappa3_2->Eval(200.) << " * epsilon^2" << std::endl;

  return 0;

}
