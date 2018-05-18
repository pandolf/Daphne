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


//
// class decleration
//

class GenParticleAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

   public:
      explicit GenParticleAnalyzer(const edm::ParameterSet&);
      ~GenParticleAnalyzer();

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

      Float_t ptHat;

      Int_t nMC;
      Int_t pdgIdMC[1000];
      Int_t statusMC[1000];
      Float_t massMC[1000];
      Int_t motherIDMC[1000];
      Float_t pxMC[1000];
      Float_t pyMC[1000];
      Float_t pzMC[1000];
      Float_t eMC[1000];
      Float_t etaMC[1000];
      Float_t phiMC[1000];

  
};

