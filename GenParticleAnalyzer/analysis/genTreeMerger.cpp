#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"




int main( int argc, char* argv[] ) {

  if( argc<3 ) {
    std::cout << "Usage: ./mergeTrees [productionName] [datasetName]" << std::endl;
    exit(1);
  }

  std::string prodName(argv[1]);
  std::string datasetName(argv[2]);

  std::string path("/eos/cms/store/group/phys_higgs/pandolf/Daphne/GenAnalysis/");
  std::string treePath(Form("%s/%s/%s", path.c_str(), prodName.c_str(), datasetName.c_str()));
  std::cout << "-> Starting: " << datasetName << std::endl;
  std::cout << "-> Looking for gentrees in: " << treePath << std::endl;

  TChain* chain = new TChain("gentree");



  chain->Add( Form("%s/genTree_*.root/demo/gentree", treePath.c_str()) );


  std::cout << "   Added trees. Total entries: " << chain->GetEntries() << std::endl;


  std::string outdir( Form("%s/%s", prodName.c_str(), datasetName.c_str()) );
  system( Form("mkdir -p %s", outdir.c_str()) );


  std::string outfileName( Form("%s/mergedTree.root", outdir.c_str()) );
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );
  outfile->cd();

  std::cout << "   Merging..." << std::endl;
  TTree* tree = chain->CloneTree();

  outfile->cd();
  tree->Write();
  outfile->Close();

  std::cout << "   Written merged tree in: " << outfileName.c_str() << std::endl;


  return 0;

}
