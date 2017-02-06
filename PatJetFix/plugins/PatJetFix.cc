// -*- C++ -*-
//
// Package:    ECALSlewFix/PatJetFix
// Class:      PatJetFix
// 
/**\class PatJetFix PatJetFix.cc ECALSlewFix/PatJetFix/plugins/PatJetFix.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rishi Patel
//         Created:  Sun, 05 Feb 2017 21:04:50 GMT
//
//


// system include files
#include <memory>
#include <algorithm>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TVector3.h"



//
// class declaration
//

class PatJetFix : public edm::stream::EDProducer<> {
   public:
      explicit PatJetFix(const edm::ParameterSet&);
      ~PatJetFix();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
	edm::InputTag met_;
	edm::InputTag metFix_;
	edm::InputTag pho_;
	edm::InputTag ele_;
	edm::InputTag jets_;
	edm::EDGetTokenT<edm::View<pat::Jet> > JetTok_;
	edm::EDGetTokenT<edm::View<pat::Electron> >  eleToken_;
	edm::EDGetTokenT<edm::View<pat::Photon> >  phoToken_;
	edm::EDGetTokenT<edm::View<pat::Electron> >  eleFixToken_;
	edm::EDGetTokenT<edm::View<pat::Photon>  >  phoFixToken_;
	edm::EDGetTokenT<edm::View<pat::PackedCandidate>  >  pCandToken_;
	edm::InputTag phoFix_;
	edm::InputTag eleFix_;
	edm::InputTag pCand_;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
PatJetFix::PatJetFix(const edm::ParameterSet& iConfig)
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
   //now do what ever initialization is needed
   eleFix_=iConfig.getParameter<edm::InputTag>("electronsFixed");
   ele_=iConfig.getParameter<edm::InputTag>("electrons");
   pho_=iConfig.getParameter<edm::InputTag>("photons");
   phoFix_=iConfig.getParameter<edm::InputTag>("photonsFixed");
   jets_=iConfig.getParameter<edm::InputTag>("jets");
   pCandToken_ = consumes< edm::View<pat::PackedCandidate> >(pCand_);
   phoToken_ = consumes< edm::View<pat::Photon> >(pho_);
   phoFixToken_ = consumes< edm::View<pat::Photon> >(phoFix_);
   eleToken_ = consumes< edm::View<pat::Electron> >(ele_);
   eleFixToken_ = consumes< edm::View<pat::Electron> >(eleFix_);
   JetTok_=consumes< edm::View<pat::Jet> >(jets_);
   produces<std::vector<pat::Jet> >("CorrPatJets");  
}


PatJetFix::~PatJetFix()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PatJetFix::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

Handle< View< pat::PackedCandidate> > hPackedCand;
iEvent.getByToken(pCandToken_,hPackedCand);
Handle< View< pat::Photon> > hPhotonProduct;
iEvent.getByToken(phoToken_,hPhotonProduct);
Handle< View< pat::Photon> > hPhotonFixProduct;
iEvent.getByToken(phoFixToken_,hPhotonFixProduct);

Handle< View< pat::Electron> > hElectronProduct;
iEvent.getByToken(eleToken_,hElectronProduct);
Handle< View< pat::Electron> > hElectronFixProduct;
iEvent.getByToken(eleFixToken_,hElectronFixProduct);
Handle<edm::View<pat::Jet> >PFJets;
iEvent.getByToken(JetTok_, PFJets);

std::auto_ptr< std::vector<pat::Jet> > patJets ( new std::vector<pat::Jet>() );
std::vector<unsigned int>PhoPFParticles;
std::vector<unsigned int >PhoCorrection;
auto phoFix=hPhotonFixProduct->begin();
auto pho=hPhotonProduct->begin();
for(pho=hPhotonProduct->begin(); pho!=hPhotonProduct->end(); ++pho, ++phoFix){
	unsigned int phoindex=pho-hPhotonProduct->begin();
        if(pho->energy()/phoFix->energy()>1.000001 || pho->energy()/phoFix->energy()<0.999999){

        for( const edm::Ref<pat::PackedCandidateCollection> & ref : pho->associatedPackedPFCandidates() ) {
        PhoPFParticles.push_back(ref.key());
        PhoCorrection.push_back(phoindex);
        }
     }
  }
for(auto Jet = PFJets->begin(); Jet != PFJets->end(); ++Jet){
	bool EGMatch=false;
	unsigned int idx=Jet-PFJets->begin();
	edm::RefToBase<pat::Jet> jetRef = PFJets->refAt(idx);
	pat::Jet ajet(jetRef);
	unsigned int PhoIndx=-1;
	for( const edm::Ptr<reco::Candidate> & ref : Jet->getJetConstituents() ) {
	    std::vector<unsigned int>::iterator it = find (PhoPFParticles.begin(), PhoPFParticles.end(), ref.key());	
	    if(it!=PhoPFParticles.end()){
		EGMatch=true;
		PhoIndx=PhoCorrection[it-PhoPFParticles.begin()];
		break; 
		}
      }
	if(EGMatch){
		math::XYZTLorentzVector pVecShift=ajet.p4();
		//math::PtEtaPhiMLorentzVector PolarLorentzVector pVecShift=ajet.p4();
		phoFix=hPhotonFixProduct->begin();
		pho=hPhotonProduct->begin();
		std::advance(pho,PhoIndx);
		std::advance(phoFix,PhoIndx);	
		float dpx=phoFix->px()-pho->px();
		float dpy=phoFix->py()-pho->py();
		float dpz=phoFix->pz()-pho->pz();
		//float dpt=phoFix->pt()-pho->pt();
		float dE=phoFix->energy()-pho->energy();
		pVecShift.SetXYZT(pVecShift.Px()+dpx, pVecShift.Py()+dpy, pVecShift.Pz()+dpz,pVecShift.E()+dE);
		ajet.setP4(pVecShift);
         }
		patJets->push_back(ajet);
	
  }
iEvent.put(patJets,"CorrPatJets");
 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
PatJetFix::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
PatJetFix::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
PatJetFix::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
PatJetFix::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
PatJetFix::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
PatJetFix::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PatJetFix::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatJetFix);
