#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"




int main( int argc, char* argv[]) {


   if( argc < 3 ) {
     std::cout << "-> USAGE: ./drawGenTree [prodName] [dataset=\'QCD_Pt_15to30\']" << std::endl;
     exit(1);
   }


   std::string prodName(argv[1]);
   std::string dataset("QCD_Pt_15to30");
   if( argc>2 ) {
     dataset = std::string(argv[2]);
   }

   TFile* file = TFile::Open( Form("%s/%s/mergedTree.root", prodName.c_str(), dataset.c_str()) );
   TTree* tree = (TTree*)file->Get("gentree");






}
