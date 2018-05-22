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
//#include "DataFormats/JetReco/interface/PFJet.h"
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
  simTrackContainerToken_  = consumes<SimTrackContainer>(edm::InputTag("g4SimHits"));
  simVertexContainerToken_ = consumes<SimVertexContainer>(edm::InputTag("g4SimHits"));

  m_tree = fs_->make<TTree>("gentree","");

  m_tree->Branch("event",&event,"event/I");
  m_tree->Branch("ptHat",&ptHat,"ptHat/F");
  m_tree->Branch("nMC",&nMC,"nMC/I");
  m_tree->Branch("pdgIdMC",pdgIdMC,"pdgIdMC[nMC]/I");
  m_tree->Branch("statusMC",statusMC,"statusMC[nMC]/I");
  m_tree->Branch("pMC ",pMC ,"pMC[nMC]/F");
  m_tree->Branch("ptMC ",ptMC ,"ptMC[nMC]/F");
  m_tree->Branch("pzMC ",pzMC ,"pzMC[nMC]/F");
  m_tree->Branch("eMC  ",eMC  ,"eMC[nMC]/F");
  m_tree->Branch("etaMC",etaMC,"etaMC[nMC]/F");
  m_tree->Branch("phiMC",phiMC,"phiMC[nMC]/F");
  m_tree->Branch("vertR",vertR,"vertR[nMC]/F");

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

   Handle<SimTrackContainer> simTracks;
   iEvent.getByToken( simTrackContainerToken_, simTracks);

   Handle<SimVertexContainer> simVertices;
   iEvent.getByToken( simVertexContainerToken_, simVertices);

   //iEvent.getByLabel("g4SimHits", simVert_h);




   /// get MC info from GenParticleCandidates 

   //std::cout << "pfjets: " << pfJets->size() << std::endl;

   std::vector< const GenParticle* > kappas;


   for ( GenParticleCollection::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p ) {
    
     if( nMC==1000 ) break;

     if( abs(p->pdgId())!=321 ) continue;

     pdgIdMC[nMC] = p->pdgId();
     statusMC[nMC] = p->status();
     massMC[nMC] = p->mass();
     ptMC[nMC] = sqrt( p->px()*p->px() + p->py()*p->py() );
     pzMC[nMC] = p->pz();
     pMC[nMC] = sqrt( pzMC[nMC]*pzMC[nMC] + ptMC[nMC]*ptMC[nMC] );
     eMC[nMC] = p->energy();
     etaMC[nMC] = p->eta();
     phiMC[nMC] = p->phi();

     vertR[nMC] = -1.;

     const GenParticle* thisGenP =  (const GenParticle*)(&(*p));

     if( abs(pdgIdMC[nMC])==321 ) kappas.push_back( thisGenP );

     nMC++;

   }
   
   //std::cout << simVertices << simTracks << std::endl;


   std::vector<int> kappa_trackId;

   for( unsigned int i=0; i<kappas.size(); ++i ) {

     kappa_trackId.push_back(-1);

     for(SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim) {

       if(iSim->genpartIndex() != -1) {

        const reco::Candidate* p = &(*genParticles)[iSim->genpartIndex()-1];

        if( (kappas[i]->pdgId() != p->pdgId()) ||
            (kappas[i]->eta() != p->eta()) ||
            (kappas[i]->phi() != p->phi()) ||
            (kappas[i]->energy() != p->energy()) ) continue;
            
          kappa_trackId[i] = iSim->trackId();

        } //if genpartindex

     } //for sim tracks

   } //for kappa


   for( unsigned int i=0; i<kappa_trackId.size(); ++i ) {

     int vertIndex = 0;

     for(SimVertexContainer::const_iterator iVert = simVertices->begin(); iVert != simVertices->end(); ++iVert) {

       if( kappa_trackId[i] == iVert->parentIndex() ) {

         for(SimTrackContainer::const_iterator iTrack=simTracks->begin(); iTrack!=simTracks->end(); ++iTrack) {
    
           if( iTrack->vertIndex()==vertIndex ) {
          
             Float_t vertX = iVert->position().x();
             Float_t vertY = iVert->position().y();
             //Float_t vertZ = iVert->position().z();
          
             //Float_t eta = iTrack->momentum().eta();
             //Float_t theta = 2.*atan( exp(-eta));

             vertR[i] = sqrt(vertX*vertX+vertY*vertY);

             //if( (first)&&(fabs(eta)<1.4) ) {
             //  h1_mfpkappa_itCone5->Fill(R*sin(theta));
             //  first=false;
             //}
          
             //bool decayedBeforeCalo=false;

             //if( fabs(eta)<1.4 ) {
             //  if( R<129. ) decayedBeforeCalo=true;
             //} else { 
             //  if( fabs(vertZ) < 304. ) decayedBeforeCalo=true;
             //}

             //double trackE =  iTrack->momentum().e();
             //double trackPx = iTrack->momentum().px();
             //double trackPy = iTrack->momentum().py();
             //double trackPz = iTrack->momentum().pz();
             //double trackP = sqrt( trackPx*trackPx + trackPy*trackPy + trackPz*trackPz );
             //double trackPt = trackP*sin(theta);


           } //if vertindex

         } //for simtracks
       
       } // if found track ID

       ++vertIndex;

     } //for sim vertex

       
   } //for kappa track ID

  
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
