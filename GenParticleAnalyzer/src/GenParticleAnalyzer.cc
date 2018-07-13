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

#include "TLorentzVector.h"



using namespace edm;
using namespace std;
using namespace reco;


class PartLite {

 public:

  int pdgId;
  float pt;
  float eta;
  float phi;
  float m;

 private:


};


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
  m_tree->Branch("vertZ",vertZ,"vertZ[nMC]/F");
  m_tree->Branch("trkIso",trkIso,"trkIso[nMC]/I");
  m_tree->Branch("nLowP",nLowP,"nLowP[nMC]/I");
  m_tree->Branch("nCharged",nCharged,"nCharged[nMC]/I");
  m_tree->Branch("decayMode",decayMode,"decayMode[nMC]/I");
  m_tree->Branch("nDau",nDau,"nDau[nMC]/I");

  m_tree->Branch("m_ppp",m_ppp,"m_ppp[nMC]/F");

  m_tree->Branch("m_pee0",m_pee0,"m_pee0[nMC]/F");
  m_tree->Branch("m_pee1",m_pee1,"m_pee1[nMC]/F");
  m_tree->Branch("m_pee2",m_pee2,"m_pee2[nMC]/F");

  m_tree->Branch("ptDau0",ptDau0,"ptDau0[nMC]/F");
  m_tree->Branch("mDau0",mDau0,"mDau0[nMC]/F");
  m_tree->Branch("etaDau0",etaDau0,"etaDau0[nMC]/F");
  m_tree->Branch("phiDau0",phiDau0,"phiDau0[nMC]/F");
  m_tree->Branch("pdgIdDau0",pdgIdDau0,"pdgIdDau0[nMC]/I");

  m_tree->Branch("ptDau1",ptDau1,"ptDau1[nMC]/F");
  m_tree->Branch("mDau1",mDau1,"mDau1[nMC]/F");
  m_tree->Branch("etaDau1",etaDau1,"etaDau1[nMC]/F");
  m_tree->Branch("phiDau1",phiDau1,"phiDau1[nMC]/F");
  m_tree->Branch("pdgIdDau1",pdgIdDau1,"pdgIdDau1[nMC]/I");

  m_tree->Branch("ptDau2",ptDau2,"ptDau2[nMC]/F");
  m_tree->Branch("mDau2",mDau2,"mDau2[nMC]/F");
  m_tree->Branch("etaDau2",etaDau2,"etaDau2[nMC]/F");
  m_tree->Branch("phiDau2",phiDau2,"phiDau2[nMC]/F");
  m_tree->Branch("pdgIdDau2",pdgIdDau2,"pdgIdDau2[nMC]/I");

  m_tree->Branch("nTrackable"                     , &nTrackable                     , "nTrackable/I"                      );
  m_tree->Branch("nTrackable_pt1"                 , &nTrackable_pt1                 , "nTrackable_pt1/I"                  );
  m_tree->Branch("nTrackablePions"                , &nTrackablePions                , "nTrackablePions/I"                 );
  m_tree->Branch("nTrackableKappas"               , &nTrackableKappas               , "nTrackableKappas/I"                );
  m_tree->Branch("nTrackableProtons"              , &nTrackableProtons              , "nTrackableProtons/I"               );
  m_tree->Branch("nTrackableProtons_dEdx"         , &nTrackableProtons_dEdx         , "nTrackableProtons_dEdx/I"          );
  m_tree->Branch("nTrackableProtonsFromDelta"     , &nTrackableProtonsFromDelta     , "nTrackableProtonsFromDelta/I"      );
  m_tree->Branch("nTrackableProtonsFromDelta_dEdx", &nTrackableProtonsFromDelta_dEdx, "nTrackableProtonsFromDelta_dEdx/I" );

  //m_tree->Branch("ptDau",ptDau,"ptDau[nMC][nDau]/F");
  //m_tree->Branch("mDau",mDau,"mDau[nMC][nDau]/F");
  //m_tree->Branch("etaDau",etaDau,"etaDau[nMC][nDau]/F");
  //m_tree->Branch("phiDau",phiDau,"phiDau[nMC][nDau]/F");
  //m_tree->Branch("pdgIdDau",pdgIdDau,"pdgIdDau[nMC][nDau]/I");

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

       
   nTrackable = 0;
   nTrackable_pt1 = 0;
   nTrackableKappas = 0;
   nTrackablePions = 0;
   nTrackableProtons = 0;
   nTrackableProtons_dEdx = 0;
   nTrackableProtonsFromDelta = 0;
   nTrackableProtonsFromDelta_dEdx = 0;

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

   std::vector< const GenParticle* > trackedParticles;
   //std::vector<int> delta_indexes;

   int nGoodKappas=0;
   int index=0;

   //std::cout << "ISPROMPT:    PDGID:    MOM_PDGID:" << std::endl;
   for ( GenParticleCollection::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p ) {

     if( nMC==300 ) break;

     float pt = sqrt( p->px()*p->px() + p->py()*p->py() );
     bool isTrackable = (pt>0.5 && abs(p->eta()<2.4));

     if( abs(p->pdgId())==321 || abs(p->pdgId())==211 || (abs(p->pdgId())==2212 && p->mother()!=0) ) {
       if( isTrackable ) {
         nTrackable++;
         if( pt > 1.0 )
           nTrackable_pt1++;
         if( abs(p->pdgId())==321 ) {
           nTrackableKappas++;
         } else if( abs(p->pdgId())==211 ) {
           nTrackablePions++;
         }
       }
     }

     //if( p->mother()!=0 )
     //  std::cout << p->isPromptFinalState() << " " <<  p->pdgId() << " " << p->mother()->pdgId() << " " << std::endl;
     //else
     //  std::cout << p->isPromptFinalState() << " " <<  p->pdgId() << "  no MOM" << " " << std::endl;
     if( abs(p->pdgId())==321  || // K+/-
         abs(p->pdgId())==2224 || // Delta++
         abs(p->pdgId())==2214 || // Delta+
         abs(p->pdgId())==2114 || // Delta0
         abs(p->pdgId())==1114 || // Delta-
         abs(p->pdgId())==3222 || // Sigma+
         abs(p->pdgId())==3112    // Sigma-
       ) { // K+/-, Delta++, Delta+, Delta0, Delta-

       pdgIdMC[nMC] = p->pdgId();
       statusMC[nMC] = p->status();
       massMC[nMC] = p->mass();
       ptMC[nMC] = pt;
       pzMC[nMC] = p->pz();
       pMC[nMC] = sqrt( pzMC[nMC]*pzMC[nMC] + ptMC[nMC]*ptMC[nMC] );
       eMC[nMC] = p->energy();
       etaMC[nMC] = p->eta();
       phiMC[nMC] = p->phi();

       vertR[nMC] = -1.;
       vertZ[nMC] = -1.;
       trkIso[nMC] = 0;
       nLowP[nMC] = 0;
       nCharged[nMC] = 0;
       decayMode[nMC] = -1;
    
       nDau[nMC] = 0;

       const GenParticle* thisGenP =  (const GenParticle*)(&(*p));
       trackedParticles.push_back( thisGenP );

       //if( p->pdgId()== 2214 || p->pdgId()== 1114 ) delta_indexes.push_back(index);

       nMC++;

     } else if( abs(p->pdgId())==2212 ) {  // protons

       float pp = sqrt( pt*pt + p->pz()*p->pz() );
       bool is_dEdx = pp>0.1 && pp<1.7;

       if( isTrackable ) {
          
         nTrackableProtons++;
       
         int momID = (p->mother()!=0) ? p->mother()->pdgId() : 0;

         if( abs(momID)==2224 || abs(momID)==2214 || abs(momID)==2114 || abs(momID)==1114 ) {

           nTrackableProtonsFromDelta++;

           if( is_dEdx ) nTrackableProtonsFromDelta_dEdx++;

         } 

         if( is_dEdx ) nTrackableProtons_dEdx++;
  
       }  // if trackable

     }  // if protons

     index++;

   } // for genparticles
   

   for( unsigned i=0; i<trackedParticles.size(); ++i ) {

     for ( GenParticleCollection::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p ) {

       if( p->pdgId()  == trackedParticles[i]->pdgId() && 
           p->eta()    == trackedParticles[i]->eta() && 
           p->phi()    == trackedParticles[i]->phi() && 
           p->energy() == trackedParticles[i]->energy() ) continue;

       if( p->charge()==0 ) continue;

       float deltaEta = p->eta() - trackedParticles[i]->eta();
       float deltaPhi = p->phi() - trackedParticles[i]->phi();
       float pi = 3.14159;
       while( deltaPhi >  pi ) deltaPhi -= 2.*pi;
       while( deltaPhi < -pi ) deltaPhi += 2.*pi;
       float deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

       if( deltaR<0.05 ) trkIso[i] += 1;

     } // for genparticles

   } // for trackedParticles



   std::vector<int> simTrackId;

   for( unsigned int i=0; i<trackedParticles.size(); ++i ) {

     simTrackId.push_back(-1);

     if( abs(trackedParticles[i]->pdgId()) != 321  && // K+- 
         abs(trackedParticles[i]->pdgId()) != 3222 && // Sigma+
         abs(trackedParticles[i]->pdgId()) != 3112    // Sigma-
         ) continue;

     for(SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim) {

       if(iSim->genpartIndex() != -1) {

         const reco::Candidate* p = &(*genParticles)[iSim->genpartIndex()-1];

         if( (trackedParticles[i]->pdgId()  != p->pdgId()) ||
             (trackedParticles[i]->eta()    != p->eta())   ||
             (trackedParticles[i]->phi()    != p->phi())   ||
             (trackedParticles[i]->energy() != p->energy()) ) continue;
             
         simTrackId[i] = iSim->trackId();

       } //if genpartindex

     } //for sim tracks

   } //for trackedParticles



   // save daughters

   for( unsigned int i=0; i<trackedParticles.size(); ++i ) {


     // kappas first (here we need to go through the simvertex)

     if( abs(trackedParticles[i]->pdgId()) == 321 ) {

       std::vector<PartLite*> daughtersK;
       int nDaughters = 0;
       int nPiCharged = 0;
       int nPiNeutral = 0;
       int nBaryons = 0;
       int nElectrons = 0;
       int nMuons = 0;
       int nPhotons = 0;
       int nProtons = 0;

       for(SimVertexContainer::const_iterator iVert = simVertices->begin(); iVert != simVertices->end(); ++iVert) {

         if( simTrackId[i] == iVert->parentIndex() ) {

           Float_t vertx = iVert->position().x();
           Float_t verty = iVert->position().y();
           Float_t vertz = iVert->position().z();
           
           //Float_t eta = iTrack->momentum().eta();
           //Float_t theta = 2.*atan( exp(-eta));

           vertR[i] = sqrt(vertx*vertx+verty*verty);
           vertZ[i] = vertz;

           if( vertR[i] < 60. && pMC[i]<1.05 && fabs(etaMC[i])<2.4 ) nGoodKappas++;


           if( vertR[i]<60. && fabs(etaMC[i])<2.4 && fabs(vertZ[i])<2. ) {

             //std::cout << std::endl;
             //std::cout << "PDGID: " << kappas[i]->pdgId() << " e: " << kappas[i]->energy() << " eta: " << kappas[i]->eta() << " phi: " << kappas[i]->phi() << " decayed to: " << std::endl;

             // loop again on simtracks to find decay products
             for(SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim) {

               if( iSim->noVertex() ) continue; // ignore prompt particle (vertIndex = -1) 

               if( abs(iSim->vertIndex()) == iVert->vertexId() ) {

                 float e = iSim->momentum().T();
                 float x = iSim->momentum().X();
                 float y = iSim->momentum().Y();
                 float z = iSim->momentum().Z();
                 float p = sqrt( x*x + y*y + z*z );
                 //float mass = sqrt( e*e - p*p );
                 int id = iSim->type();
                 //std::cout << "  " << iSim->type() << " m: " << sqrt( e*e - p*p ) << " e: " << iSim->momentum().T() << " phi: " << iSim->momentum().Phi() << std::endl;

                 nDaughters++;
                 if( fabs(iSim->charge())>0.01 ) nCharged[i]++;

                 if( abs(id) == 211 ) nPiCharged++;
                 else if( id == 111 ) nPiNeutral++;
                 else if( abs(id)>1000 && abs(id)<9999 ) {
                   nBaryons++;
                   if( abs(id) == 2212 ) nProtons++;
                 }
                 else if( abs(id) == 11 ) nElectrons++;
                 else if( abs(id) == 13 ) nMuons++;
                 else if( id == 22 ) nPhotons++;


                 if( p<0.15 && fabs(iSim->charge())>0.01 ) nLowP[i]++;

                 //if( fabs(iSim->charge())>0.01 ) { // save only charged daughters
                   PartLite* thisPart = new PartLite();
                   thisPart->pdgId = id;
                   TLorentzVector thisp(x, y, z, e);
                   thisPart->pt  = thisp.Pt();
                   thisPart->eta = thisp.Eta();
                   thisPart->phi = thisp.Phi();
                   thisPart->m   = thisp.M();
                   //std::cout << "part:  " << thisPart->pdgId << " m: " << thisPart->m << " pt: " << thisPart->pt << " phi: " << thisPart->phi << std::endl;
                   daughtersK.push_back(thisPart);
                 //}


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


       if( (decayMode[i]==0 && (daughtersK.size()==3 || daughtersK.size()==4)) || decayMode[i] == 8 || decayMode[i] == 7 ) {

         nDau[i] = daughtersK.size();

         float m_pi = 0.1396;

         std::vector<TLorentzVector> pions;

         for( unsigned iD=0; iD<3; ++iD ) { // take three also in the case of 4 daughters (decayMode=0)
           
           TLorentzVector thisPion;
           thisPion.SetPtEtaPhiM( daughtersK[iD]->pt, daughtersK[iD]->eta, daughtersK[iD]->phi, m_pi );

           pions.push_back(thisPion);

         } // for daughters

         m_ppp [i] = computeMass( pions );
         m_pee0[i] = computeMass( pions, 0 );
         m_pee1[i] = computeMass( pions, 1 );
         m_pee2[i] = computeMass( pions, 2 );
         

         ptDau0   [i] = daughtersK[0]->pt;
         etaDau0  [i] = daughtersK[0]->eta;
         phiDau0  [i] = daughtersK[0]->phi;
         mDau0    [i] = daughtersK[0]->m;
         pdgIdDau0[i] = daughtersK[0]->pdgId;

         ptDau1   [i] = daughtersK[1]->pt;
         etaDau1  [i] = daughtersK[1]->eta;
         phiDau1  [i] = daughtersK[1]->phi;
         mDau1    [i] = daughtersK[1]->m;
         pdgIdDau1[i] = daughtersK[1]->pdgId;

         ptDau2   [i] = daughtersK[2]->pt;
         etaDau2  [i] = daughtersK[2]->eta;
         phiDau2  [i] = daughtersK[2]->phi;
         mDau2    [i] = daughtersK[2]->m;
         pdgIdDau2[i] = daughtersK[2]->pdgId;


       } // if interesting decayMode

       //if( decayMode[i] == 13 ) {
       //  std::cout << "nDaughters: " <<  nDaughters << std::endl;
       //  std::cout << "nPiCharged: " <<  nPiCharged << std::endl;
       //  std::cout << "nPiNeutral: " <<  nPiNeutral << std::endl;
       //  std::cout << "nBaryons  : " <<  nBaryons   << std::endl;
       //  std::cout << "nElectrons: " <<  nElectrons << std::endl;
       //  std::cout << "nMuons    : " <<  nMuons     << std::endl;
       //  std::cout << "nPhotons  : " <<  nPhotons   << std::endl;
       //}

       

     // and now deltas (here it's easier, we can use the genparticles)
     
     } else if( abs(trackedParticles[i]->pdgId())==2224 || abs(trackedParticles[i]->pdgId())==2214 || abs(trackedParticles[i]->pdgId())==2114 || abs(trackedParticles[i]->pdgId())==1114  ) { 

       std::vector<PartLite*> daughtersD;
     
       // second loop to find daughters
       for( GenParticleCollection::const_iterator p2 = genParticles->begin(); p2 != genParticles->end(); ++p2 ) { 

         //if( &(*p2) == trackedParticles[i] ) continue;
         if( p2->mother()==0 ) continue;
         if( p2->mother()->pdgId()  !=  trackedParticles[i]->pdgId() ) continue;

         if( p2->mother()->pdgId()  == trackedParticles[i]->pdgId()  &&
             p2->mother()->energy() == trackedParticles[i]->energy() &&
             p2->mother()->px()     == trackedParticles[i]->px()     &&
             p2->mother()->py()     == trackedParticles[i]->py()     &&
             p2->mother()->pz()     == trackedParticles[i]->pz()     ) {

           PartLite* thisPart = new PartLite();
           thisPart->pdgId = p2->pdgId();
           float x = p2->px();
           float y = p2->py();
           float z = p2->pz();
           float e = p2->energy();
           TLorentzVector thisp(x, y, z, e);
           thisPart->pt  = thisp.Pt();
           thisPart->eta = thisp.Eta();
           thisPart->phi = thisp.Phi();
           thisPart->m   = thisp.M();

           daughtersD.push_back(thisPart);

         }  // if found mother

       } // second loop to find daughters
       

       if( daughtersD.size()==2 ) {
         if(      abs(daughtersD[0]->pdgId * daughtersD[1]->pdgId) == 466732    ) decayMode[i]=1000; // p pi+-   (proton=2212, pi=211)
         else if( abs(daughtersD[0]->pdgId * daughtersD[1]->pdgId) == 48664     ) decayMode[i]=1001; // p gamma  (proton=2212, photon=22)
         else if( abs(daughtersD[0]->pdgId * daughtersD[1]->pdgId) == 245532    ) decayMode[i]=1002; // p pi0    (proton=2212, pizero=111)
         else if( abs(daughtersD[0]->pdgId)==2112 ||  abs(daughtersD[1]->pdgId) ) decayMode[i]=1003; // n X      (neutron=2112) 
         else {
           std::cout << "Unexpected decay mode!!" << std::endl;
           std::cout << "  " << daughtersD[0]->pdgId << " " << daughtersD[1]->pdgId << std::endl;
         }
       } else {
         std::cout << "Delta decaying to N!=2 particles??" << std::endl;
         for( unsigned iD=0; iD<daughtersD.size(); ++iD ) {
           std::cout << "   " << daughtersD[iD]->pdgId << std::endl;
         }
       }

       nDau[i] = daughtersD.size();

       ptDau0   [i] = daughtersD[0]->pt;
       etaDau0  [i] = daughtersD[0]->eta;
       phiDau0  [i] = daughtersD[0]->phi;
       mDau0    [i] = daughtersD[0]->m;
       pdgIdDau0[i] = daughtersD[0]->pdgId;

       ptDau1   [i] = daughtersD[1]->pt;
       etaDau1  [i] = daughtersD[1]->eta;
       phiDau1  [i] = daughtersD[1]->phi;
       mDau1    [i] = daughtersD[1]->m;
       pdgIdDau1[i] = daughtersD[1]->pdgId;

       ptDau2   [i] = 0.;
       etaDau2  [i] = 1000.;
       phiDau2  [i] = 0.;
       mDau2    [i] = 0.;
       pdgIdDau2[i] = 0;

     } // if deltas

   } //for trackedParticles


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


float GenParticleAnalyzer::computeMass( const std::vector<TLorentzVector>& pions, int index_ele ) {

  float m_ele = 0.000511;
  TLorentzVector sum(0.,0.,0.,0.);

  for( int i=0; i<(int)pions.size(); ++i ) {

    if( index_ele>=0 && i!=index_ele ) {

      TLorentzVector ele;
      ele.SetPtEtaPhiM( pions[i].Pt(), pions[i].Eta(), pions[i].Phi(), m_ele );

      sum += ele;

    } else {

      sum += pions[i];

    }

  } // for pions

  return sum.M();

}



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
