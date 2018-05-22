// -*- C++ -*-
//
// Package:    Daphne/GenParticleAnalyzer
// Class:      GenParticleAnalyzer
// 
/**\class GenParticleAnalyzer GenParticleAnalyzer.cc Daphne/GenParticleAnalyzer/plugins/GenParticleAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesco Pandolfi
//         Created:  Thu, 17 May 2018 12:05:16 GMT
//
//


//
// constructors and destructor
//


#include <memory>

// user include files
#include "Daphne/GenParticleAnalyzer/interface/GenParticleAnalyzer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"


//#include "MagneticField/Engine/interface/MagneticField.h"
//#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
//#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
//#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
//#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
//#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
//#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
//#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
//#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
//#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
//#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"

//#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"


using namespace edm;
using namespace std;
using namespace reco;


GenParticleAnalyzer::GenParticleAnalyzer(const edm::ParameterSet& conf)
{
  //MCTruthCollection_ = iConfig.getUntrackedParameter<edm::InputTag>("MCTruthCollection");
  //trackTags_ = iConfig.getUntrackedParameter<edm::InputTag>("tracks");

  genParticleCollectionToken_ = consumes< GenParticleCollection >(edm::InputTag("genParticles"));
  PFJetCollectionToken_ = consumes<PFJetCollection>(edm::InputTag("ak4PFJets"));

  m_tree = fs_->make<TTree>("gentree","");

  m_tree->Branch("event",&event,"event/I");
  m_tree->Branch("ptHat",&ptHat,"ptHat/F");
  m_tree->Branch("nMC",&nMC,"nMC/I");
  m_tree->Branch("pdgIdMC",pdgIdMC,"pdgIdMC[nMC]/I");
  m_tree->Branch("statusMC",statusMC,"statusMC[nMC]/I");
  m_tree->Branch("pxMC ",pxMC ,"pxMC[nMC]/F");
  m_tree->Branch("pyMC ",pyMC ,"pyMC[nMC]/F");
  m_tree->Branch("pzMC ",pzMC ,"pzMC[nMC]/F");
  m_tree->Branch("eMC  ",eMC  ,"eMC[nMC]/F");
  m_tree->Branch("etaMC",etaMC,"etaMC[nMC]/F");
  m_tree->Branch("phiMC",phiMC,"phiMC[nMC]/F");

  event = 0;  
  
}


GenParticleAnalyzer::~GenParticleAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GenParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   nMC = 0;

   // get generated pt hat
   //Handle<double> genEventScale; 
   //iEvent.getByLabel( "genEventScale", genEventScale ); 
   //ptHat = genEventInfo->qScale();   

   //Handle<GenEventInfoProduct> genEventInfo; 
   //iEvent.getByLabel( "generator", genEventInfo ); 

   Handle< GenParticleCollection > genParticles;
   iEvent.getByToken( genParticleCollectionToken_, genParticles );

   //Handle<PFJetCollection> pfJets;
   //iEvent.getByToken( PFJetCollectionToken_, pfJets );


   /// get MC info from GenParticleCandidates 

   //std::cout << "pfjets: " << pfJets->size() << std::endl;
   std::cout << "genParticles: " << genParticles->size() << std::endl;

   std::vector< const GenParticle* > kappas;


   for ( GenParticleCollection::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p ) {
    
     if( nMC==1000 ) break;

     pdgIdMC[nMC] = p->pdgId();
     statusMC[nMC] = p->status();
     massMC[nMC] = p->mass();
     pxMC[nMC] = p->px();	 
     pyMC[nMC] = p->py();	 
     pzMC[nMC] = p->pz();	 
     eMC[nMC] = p->energy();	 
     etaMC[nMC] = p->eta();	 
     phiMC[nMC] = p->phi();	 

     //std::cout << "pdgIdMC[nMC] : " <<  pdgIdMC[nMC] << std::endl;
     //std::cout << "statusMC[nMC] : " << statusMC[nMC] << std::endl;
     //std::cout << "massMC[nMC] : " <<   massMC[nMC] << std::endl;
     //std::cout << "pxMC[nMC] : " <<     pxMC[nMC] << std::endl;
     //std::cout << "pyMC[nMC] : " <<     pyMC[nMC] << std::endl;
     //std::cout << "pzMC[nMC] : " <<     pzMC[nMC] << std::endl;
     //std::cout << "eMC[nMC] : " <<      eMC[nMC] << std::endl;
     //std::cout << "etaMC[nMC] : " <<    etaMC[nMC] << std::endl;
     //std::cout << "phiMC[nMC] : " <<    phiMC[nMC] << std::endl;

     const GenParticle* thisGenP =  (const GenParticle*)(&(*p));

     if( abs(pdgIdMC[nMC])==321 ) kappas.push_back( thisGenP );

     nMC++;

   }
   
   //Handle<SimTrackContainer> simTracks_h;
   //const SimTrackContainer* simTracks;
   //iEvent.getByLabel("g4SimHits", simTracks_h);
   //simTracks = simTracks_h.product();

   //Handle<SimVertexContainer> simVert_h;
   //const SimVertexContainer* simVertices;
   //iEvent.getByLabel("g4SimHits", simVert_h);
   //simVertices = simVert_h.product();

   //std::cout << simVertices << simTracks << std::endl;


   //// ------------------ BEGIN K0Short 
   //
   //std::vector<int> K0S_trackId;

   //for( int i=0; i<K0S.size(); ++i ) {

   //  K0S_trackId.push_back(-1);

   //  for(SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim) {

   //    if(iSim->genpartIndex() != -1) {

   //     const reco::Candidate* p = &(*genParticles)[iSim->genpartIndex()-1];

   //     if( (K0S[i]->pdgId() != p->pdgId()) ||
   //         (K0S[i]->eta() != p->eta()) ||
   //         (K0S[i]->phi() != p->phi()) ||
   //         (K0S[i]->energy() != p->energy()) ) continue;
   //         
   //       K0S_trackId[i] = iSim->trackId();

   //     } //if genpartindex

   //  } //for sim tracks

   //} //for K0s


   //for( int i=0; i<K0S_trackId.size(); ++i ) {

   //  bool first=true;

   //  int vertIndex = 0;

   //  Float_t remainingE = K0S[i]->energy();
   //  Float_t remainingPx = K0S[i]->px();
   //  Float_t remainingPy = K0S[i]->py();
   //  Float_t remainingPz = K0S[i]->pz();

   //  for(SimVertexContainer::const_iterator iVert = simVertices->begin(); iVert != simVertices->end(); ++iVert) {

   //    if( K0S_trackId[i] == iVert->parentIndex() ) {

   //      for(SimTrackContainer::const_iterator iTrack=simTracks->begin(); iTrack!=simTracks->end(); ++iTrack) {
   // 
   //        if( iTrack->vertIndex()==vertIndex ) {
   //       
   //          Float_t vertX = iVert->position().x();
   //          Float_t vertY = iVert->position().y();
   //          Float_t vertZ = iVert->position().z();
   //       
   //          Float_t eta = iTrack->momentum().eta();
   //          Float_t theta = 2.*atan( exp(-eta));

   //          Float_t R = sqrt(vertX*vertX+vertY*vertY);
   //          if( (first)&&(fabs(eta)<1.4) ) {
   //            h1_mfpK0S_itCone5->Fill(R*sin(theta));
   //            first=false;
   //          }
   //       
   //          bool decayedBeforeCalo=false;

   //          if( fabs(eta)<1.4 ) {
   //            if( R<129. ) decayedBeforeCalo=true;
   //          } else { 
   //            if( fabs(vertZ) < 304. ) decayedBeforeCalo=true;
   //          }

   //          double trackE =  iTrack->momentum().e();
   //          double trackPx = iTrack->momentum().px();
   //          double trackPy = iTrack->momentum().py();
   //          double trackPz = iTrack->momentum().pz();
   //          double trackP = sqrt( trackPx*trackPx + trackPy*trackPy + trackPz*trackPz );
   //          double trackPt = trackP*sin(theta);

   //          if( (decayedBeforeCalo) && (iTrack->type()==111) ) { //pizeros
   //              nPhotonsGen_itCone5[nJetGen_itCone5] += 2;  //boh
   //              ePhotonsGen_itCone5[nJetGen_itCone5]  += trackE;
   //              pxPhotonsGen_itCone5[nJetGen_itCone5] += trackPx;
   //              pyPhotonsGen_itCone5[nJetGen_itCone5] += trackPy;
   //              pzPhotonsGen_itCone5[nJetGen_itCone5] += trackPz;
   //          } else if( (R<30.) && (fabs(eta)<2.5) && (iTrack->charge()!=0) && (trackPt>0.1) ) { //min it tracking (4 step) pt = 100 mev? it's 75MeV, but it's fine like this
   //               nTracksGen_itCone5[nJetGen_itCone5] += 1;
   //               eTracksGen_itCone5[nJetGen_itCone5] += trackE;
   //              pxTracksGen_itCone5[nJetGen_itCone5] += trackPx;
   //              pyTracksGen_itCone5[nJetGen_itCone5] += trackPy;
   //              pzTracksGen_itCone5[nJetGen_itCone5] += trackPz;
   //          } else {
   //            nNeutralHadronsGen_itCone5[nJetGen_itCone5] += 1;
   //            eNeutralHadronsGen_itCone5[nJetGen_itCone5]  += trackE;
   //            pxNeutralHadronsGen_itCone5[nJetGen_itCone5] += trackPx;
   //            pyNeutralHadronsGen_itCone5[nJetGen_itCone5] += trackPy;
   //            pzNeutralHadronsGen_itCone5[nJetGen_itCone5] += trackPz;
   //          }

   //          remainingE  -= trackE;
   //          remainingPx -= trackPx;
   //          remainingPy -= trackPy;
   //          remainingPz -= trackPz;

   //        } //if vertindex

   //      } //for simtracks
   //    
   //    } // if found track ID

   //    ++vertIndex;

   //  } //for sim vertex

   //  if( remainingE > 0.001 ) {
   //     nNeutralHadronsGen_itCone5[nJetGen_itCone5] += 1;
   //     eNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingE;
   //    pxNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingPx;
   //    pyNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingPy;
   //    pzNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingPz;
   //  }
   //    
   //} //for K0S track ID

  
   event++;  
   m_tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
GenParticleAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenParticleAnalyzer::endJob() {

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenParticleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(GenParticleAnalyzer);
