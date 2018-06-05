// -*- C++ -*-
//
// Package:    Daphne/DaphneAnalyzer
// Class:      DaphneAnalyzer
// 
/**\class DaphneAnalyzer DaphneAnalyzer.cc Daphne/DaphneAnalyzer/plugins/DaphneAnalyzer.cc

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
#include "Daphne/DaphneAnalyzer/interface/DaphneAnalyzer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TLorentzVector.h"



using namespace edm;
using namespace std;
using namespace reco;




DaphneAnalyzer::DaphneAnalyzer(const edm::ParameterSet& conf)
{

  vertexContainerToken_ = consumes< vector<reco::Vertex> >(edm::InputTag("offlinePrimaryVertices"));


  m_tree = fs_->make<TTree>("dtree","");

  m_tree->Branch("nVert"         ,&nVert         ,"nVert/I"         );
  m_tree->Branch("nTracksPV"     ,&nTracksPV     ,"nTracksPV/I"     );
  m_tree->Branch("nTracksPV_min0",&nTracksPV_min0,"nTracksPV_min0/I");
  m_tree->Branch("nTracksPU"     ,&nTracksPU     ,"nTracksPU/F"     );
  m_tree->Branch("nTracksPU_min0",&nTracksPU_min0,"nTracksPU_min0/F");

  event = 0;  
  
}


DaphneAnalyzer::~DaphneAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
DaphneAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


   Handle< vector<reco::Vertex> > vertices;
   iEvent.getByToken( vertexContainerToken_, vertices);


   nVert = vertices->size();

   nTracksPU      = 0.;
   nTracksPU_min0 = 0.;

   bool first = true;

   for( vector<reco::Vertex>::const_iterator v=vertices->begin(); v != vertices->end(); ++v ) {

     if( first ) {

       nTracksPV      = v->nTracks();
       nTracksPV_min0 = v->nTracks(0.);

       first = false;

     } else {
  
       nTracksPU      += v->nTracks();
       nTracksPU_min0 += v->nTracks(0.);

     }

   } // for vertices

   if( nVert>1 ) {

     nTracksPU      /= ( nVert-1. );
     nTracksPU_min0 /= ( nVert-1. );

   }

   event++;  
   m_tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
DaphneAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DaphneAnalyzer::endJob() {

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------



void
DaphneAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(DaphneAnalyzer);
