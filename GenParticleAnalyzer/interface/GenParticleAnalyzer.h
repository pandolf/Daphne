#ifndef GENPARTICLEANALYZER_H
#define GENPARTICLEANALYZER_H


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
      TH1D * m_h1_nGoodKappas ;

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
      Int_t trkIso[300];
      Int_t nLowP[300];
      Int_t nCharged[300];
      Int_t decayMode[300];
      Int_t   nDau[300];

      Float_t m_ppp[300];
      Float_t m_pee0[300];
      Float_t m_pee1[300];
      Float_t m_pee2[300];
      
      Float_t ptDau0[300];
      Float_t mDau0[300];
      Float_t etaDau0[300];
      Float_t phiDau0[300];
      Int_t   pdgIdDau0[300];

      Float_t ptDau1[300];
      Float_t mDau1[300];
      Float_t etaDau1[300];
      Float_t phiDau1[300];
      Int_t   pdgIdDau1[300];

      Float_t ptDau2[300];
      Float_t mDau2[300];
      Float_t etaDau2[300];
      Float_t phiDau2[300];
      Int_t   pdgIdDau2[300];

      Int_t nTrackable;
      Int_t nTrackableKappas;
      Int_t nTrackablePions;
      Int_t nTrackableProtons;
      Int_t nTrackableProtons_dEdx;
      Int_t nTrackableProtonsFromDelta;
      Int_t nTrackableProtonsFromDelta_dEdx;

      //Float_t ptDau[300][100];
      //Float_t mDau[300][100];
      //Float_t etaDau[300][100];
      //Float_t phiDau[300][100];
      //Float_t pdgIdDau[300][100];

      edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionToken_;
      edm::EDGetTokenT<edm::SimTrackContainer > simTrackContainerToken_;
      edm::EDGetTokenT<edm::SimVertexContainer> simVertexContainerToken_;
  
};


#endif
