#ifndef DAPHNEANALYZER_H
#define DAPHNEANALYZER_H



// system include files
#include <memory>
#include <iostream>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"




//
// class decleration
//

class DaphneAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

   public:
      explicit DaphneAnalyzer(const edm::ParameterSet&);
      ~DaphneAnalyzer();

      float computeMass( const std::vector<TLorentzVector>& pions, int index_ele=-1 );

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      //edm::InputTag trackTags_; 
      //std::string filename_; 
      edm::Service<TFileService> fs_;

      TTree * m_tree ;

      Int_t event;


      edm::EDGetTokenT< std::vector<reco::Vertex> > vertexContainerToken_;
  
};


#endif
