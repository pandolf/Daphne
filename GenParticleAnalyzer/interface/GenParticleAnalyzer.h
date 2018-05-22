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


#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"


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
      Int_t pdgIdMC[300];
      Int_t statusMC[300];
      Float_t massMC[300];
      Int_t motherIDMC[300];
      Float_t pMC[300];
      Float_t ptMC[300];
      Float_t pzMC[300];
      Float_t eMC[300];
      Float_t etaMC[300];
      Float_t phiMC[300];
      Float_t vertR[300];

      edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionToken_;
      edm::EDGetTokenT<edm::SimTrackContainer > simTrackContainerToken_;
      edm::EDGetTokenT<edm::SimVertexContainer> simVertexContainerToken_;
  
};

