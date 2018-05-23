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

#include "SimDataFormats/Track/interface/CoreSimTrack.h"
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


//class PartLite {
//
// public:
//
//  int pdgId;
//  float p;
//
// private:
//
//
//};


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
  m_tree->Branch("pMC",pMC ,"pMC[nMC]/F");
  m_tree->Branch("ptMC",ptMC ,"ptMC[nMC]/F");
  m_tree->Branch("pzMC",pzMC ,"pzMC[nMC]/F");
  m_tree->Branch("eMC",eMC  ,"eMC[nMC]/F");
  m_tree->Branch("etaMC",etaMC,"etaMC[nMC]/F");
  m_tree->Branch("phiMC",phiMC,"phiMC[nMC]/F");
  m_tree->Branch("vertR",vertR,"vertR[nMC]/F");
  m_tree->Branch("trkIso",trkIso,"trkIso[nMC]/I");
  m_tree->Branch("nLowP",nLowP,"nLowP[nMC]/I");
  m_tree->Branch("decayMode",decayMode,"decayMode[nMC]/I");

  m_h1_nGoodKappas = fs_->make<TH1D>("nGoodKappas", "", 6, -0.5, 5.5 );

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


   //Handle<GenEventInfoProduct> genEventInfo; 
   //iEvent.getByLabel( "generator", genEventInfo ); 
   //ptHat = genEventInfo->qScale();   


   Handle<GenParticleCollection> genParticles;
   iEvent.getByToken( genParticleCollectionToken_, genParticles );

   Handle<SimTrackContainer> simTracks;
   iEvent.getByToken( simTrackContainerToken_, simTracks);

   Handle<SimVertexContainer> simVertices;
   iEvent.getByToken( simVertexContainerToken_, simVertices);




   /// get MC info from GenParticleCandidates 

   std::vector< const GenParticle* > kappas;
   int nGoodKappas=0;

   for ( GenParticleCollection::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p ) {

     if( nMC==300 ) break;

     if( abs(p->pdgId())!=321 ) continue; // K+/-

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
     trkIso[nMC] = 0;
     nLowP[nMC] = 0;
     decayMode[nMC] = -1;

     const GenParticle* thisGenP =  (const GenParticle*)(&(*p));

     kappas.push_back( thisGenP );

     nMC++;

   }
   

   for( unsigned i=0; i<kappas.size(); ++i ) {

     for ( GenParticleCollection::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p ) {

       if( p->pdgId()  == kappas[i]->pdgId() && 
           p->eta()    == kappas[i]->eta() && 
           p->phi()    == kappas[i]->phi() && 
           p->energy() == kappas[i]->energy() ) continue;

       if( p->charge()==0 ) continue;

       float deltaEta = p->eta() - kappas[i]->eta();
       float deltaPhi = p->phi() - kappas[i]->phi();
       float pi = 3.14159;
       while( deltaPhi >  pi ) deltaPhi -= 2.*pi;
       while( deltaPhi < -pi ) deltaPhi += 2.*pi;
       float deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

       if( deltaR<0.05 ) trkIso[i] += 1;

     } // for genparticles

   } // for kappas



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



   // look for sim vertex:

   for( unsigned int i=0; i<kappa_trackId.size(); ++i ) {

     //std::vector<PartLite> daughters;
     int nDaughters = 0;
     int nPiCharged = 0;
     int nPiNeutral = 0;
     int nBaryons = 0;
     int nElectrons = 0;
     int nMuons = 0;
     int nPhotons = 0;

     for(SimVertexContainer::const_iterator iVert = simVertices->begin(); iVert != simVertices->end(); ++iVert) {

       if( kappa_trackId[i] == iVert->parentIndex() ) {

         Float_t vertX = iVert->position().x();
         Float_t vertY = iVert->position().y();
         //Float_t vertZ = iVert->position().z();
         
         //Float_t eta = iTrack->momentum().eta();
         //Float_t theta = 2.*atan( exp(-eta));

         vertR[i] = sqrt(vertX*vertX+vertY*vertY);

         if( vertR[i] < 60. && pMC[i]<1.05 && fabs(etaMC[i])<2.4 ) nGoodKappas++;


         if( vertR[i]<60. && fabs(etaMC[i])<2.4 ) {

           //std::cout << std::endl;
           //std::cout << "PDGID: " << kappas[i]->pdgId() << " e: " << kappas[i]->energy() << " eta: " << kappas[i]->eta() << " phi: " << kappas[i]->phi() << " decayed to: " << std::endl;

           // loop again on simtracks to find decay products
           for(SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim) {

             if( iSim->noVertex() ) continue; // ignore prompt particle (vertIndex = -1) 

             if( abs(iSim->vertIndex()) == iVert->vertexId() ) {

               //float e = iSim->momentum().T();
               float x = iSim->momentum().X();
               float y = iSim->momentum().Y();
               float z = iSim->momentum().Z();
               float p = sqrt( x*x + y*y + z*z );
               //float mass = sqrt( e*e - p*p );
               int id = iSim->type();
               //std::cout << "  " << iSim->type() << " m: " << mass << " e: " << iSim->momentum().T() << " phi: " << iSim->momentum().Phi() << std::endl;
               //PartLite thisPart;
               //thisPart.pdgId = id;
               //thisPart.p = p;
               //daughters.push_back(thisPart);

               nDaughters++;
               if( abs(id) == 211 ) nPiCharged++;
               else if( id == 111 ) nPiNeutral++;
               else if( id > 1000 ) nBaryons++;
               else if( abs(id) == 11 ) nElectrons++;
               else if( abs(id) == 13 ) nMuons++;
               else if( id == 22 ) nPhotons++;

               if( p<0.15 && fabs(iSim->charge())>0.01 ) nLowP[i]++;

             } // if found vertex

           } // for simtracks

         } // if vertR<60.

       } // if found track ID

     } //for sim vertex

     
     if( nBaryons>0 )                                           decayMode[i] =  0; // nuclear interaction
     else if( nDaughters==1 && nElectrons==1 )                  decayMode[i] =  1; // e nu (neutrino not tracked)
     else if( nDaughters==1 && nMuons    ==1 )                  decayMode[i] =  2; // m nu (neutrino not tracked)
     else if( nDaughters==2 && nElectrons==1 && nPiNeutral==1 ) decayMode[i] =  3; // pi0 e nu
     else if( nDaughters==2 && nMuons    ==1 && nPiNeutral==1 ) decayMode[i] =  4; // pi0 m nu
     else if( nDaughters==2 && nPiCharged==1 && nPiNeutral==1 ) decayMode[i] =  5; // pi+ pi0
     else if( nDaughters==3 && nPiCharged==1 && nPiNeutral==2 ) decayMode[i] =  6; // pi+ pi0 pi0
     else if( nDaughters==3 && nPiCharged==3                  ) decayMode[i] =  7; // pi+ pi- pi+
     else if( nDaughters==3 && nPiCharged==1 && nElectrons==2 ) decayMode[i] =  8; // pi+ e+ e-
     else if( nDaughters==2 && nMuons    ==1 && nPhotons  ==1 ) decayMode[i] =  9; // mu nu gamma
     else if( nDaughters==2 && nElectrons==1 && nPhotons  ==1 ) decayMode[i] = 10; // e nu gamma
     else if( nDaughters==3 && nElectrons==2 && nMuons    ==1 ) decayMode[i] = 11; // m nu e+ e-
     else if( nDaughters==3 && nElectrons==3 )                  decayMode[i] = 12; // e nu e+ e-
     else if( nDaughters>0  )                                   decayMode[i] = 13; // other

     //if( decayMode[i] == 13 ) {
     //  std::cout << "nDaughters: " <<  nDaughters << std::endl;
     //  std::cout << "nPiCharged: " <<  nPiCharged << std::endl;
     //  std::cout << "nPiNeutral: " <<  nPiNeutral << std::endl;
     //  std::cout << "nBaryons  : " <<  nBaryons   << std::endl;
     //  std::cout << "nElectrons: " <<  nElectrons << std::endl;
     //  std::cout << "nMuons    : " <<  nMuons     << std::endl;
     //  std::cout << "nPhotons  : " <<  nPhotons   << std::endl;
     //}

       
   } //for kappa track ID


   m_h1_nGoodKappas->Fill( nGoodKappas );
  
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
