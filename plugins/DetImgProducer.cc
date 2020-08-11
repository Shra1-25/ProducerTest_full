// -*- C++ -*-
//
// Package:    ProdTutorial/DetImgProducer
// Class:      DetImgProducer
// 
/**\class ProducerTest DetImgProducer.cc ProdTutorial/ProducerTest/plugins/DetImgProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shravan Chaudhari
//         Created:  Fri, 05 Jun 2020 16:11:58 GMT
//
//


// system include files
/*#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"*/
#include "ProdTutorial/ProducerTest/plugins/DetImgProducer.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

using namespace tensorflow;
using namespace std;


//
// class declaration
//

/*class DetImgProducer : public edm::stream::EDProducer<> {
   public:
      explicit DetImgProducer(const edm::ParameterSet&);
      ~ProducerTest();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
       
};*/

DetImgProducer::DetImgProducer(const edm::ParameterSet& iConfig)
{
 EBRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
 photonCollectionT_ = consumes<PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonCollection"));
 HBHERecHitCollectionT_  = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));
 EERecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEERecHitCollection"));
 //TRKRecHitCollectionT_   = consumes<TrackingRecHitCollection>(iConfig.getParameter<edm::InputTag>("trackRecHitCollection"));
 trackCollectionT_       = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackCollection"));
 vertexCollectionT_       = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
 jetCollectionT_         = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("ak4PFJetCollection"));
 TRKRecHitCollectionT_   = consumes<TrackingRecHitCollection>(iConfig.getParameter<edm::InputTag>("trackRecHitCollection"));
 genParticleCollectionT_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
 genJetCollectionT_      = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetCollection"));
 pfCollectionT_          = consumes<PFCollection>(iConfig.getParameter<edm::InputTag>("pfCollection"));
 recoJetsT_              = consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("recoJetsForBTagging"));
 jetTagCollectionT_      = consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("jetTagCollection"));
 ipTagInfoCollectionT_   = consumes<std::vector<reco::CandIPTagInfo> > (iConfig.getParameter<edm::InputTag>("ipTagInfoCollection"));
 
 /*mode_      = iConfig.getParameter<std::string>("mode");
 minJetPt_  = iConfig.getParameter<double>("minJetPt");
 maxJetEta_ = iConfig.getParameter<double>("maxJetEta");
 z0PVCut_   = iConfig.getParameter<double>("z0PVCut");*/
 
 //std::cout << " >> Mode set to " << mode_ << std::endl;
 /*if ( mode_ == "JetLevel" ) {
   doJets_ = true;
   nJets_ = iConfig.getParameter<int>("nJets");
   std::cout << "\t>> nJets set to " << nJets_ << std::endl;
 } else if ( mode_ == "EventLevel" ) {
   doJets_ = false;
 } else {
   std::cout << " >> Assuming EventLevel Config. " << std::endl;
   doJets_ = false;
 }
 edm::Service<TFileService> fs;
 RHTree = fs->make<TTree>("RHTree", "RecHit tree");
  if ( doJets_ ) {
    branchesEvtSel_jet( RHTree, fs );
  } else {
    branchesEvtSel( RHTree, fs );
  }*/

 //usesResource("TFileService");
 edm::Service<TFileService> fs;
 RHTree = fs->make<TTree>("RHTree", "RecHit tree");
 //RHTree->Branch("SC_iphi", &vIphi_Emax_);
 //RHTree->Branch("SC_ieta", &vIeta_Emax_);
 branchesEB           ( RHTree, fs );
 //branchesPhotonSel ( RHTree, fs );
 branchesHBHE (RHTree, fs );
 branchesECALstitched (RHTree, fs);
 branchesTracksAtECALstitched (RHTree, fs);
 std::cout<<"Branches done "<<std::endl;
 
 //produces<float>("photonClasses").setBranchAlias("PhotonClass");
 produces<std::vector<float>>("EBenergy");
 produces<std::vector<float>>("HBHEenergy");
 produces<std::vector<float>>("HBHEenergyEB");
 produces<std::vector<float>>("ECALstitchedenergy");
 produces<std::vector<float>>("TracksAtECALstitched");
 //produces<std::vector<int>>("JetSeedieta");
 //produces<std::vector<int>>("JetSeediphi");
 //if (!fw) { return; }
}


DetImgProducer::~DetImgProducer()
{
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DetImgProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //std::cout<<"New Event started"<<std::endl;
   using namespace edm;
   nTotal++;
   // ----- Apply event selection cuts ----- //
   /*bool passedSelection = false;
   if ( doJets_ ) {
     std::cout<<" >> doJets set"<<std::endl;
     passedSelection = runEvtSel_jet( iEvent, iSetup );
     std::cout<<" >> Size of JetSeed vector (JetSeed_eta_size, JetSeed_phi_size) is: ("<<vJetSeed_ieta_.size()<<", "<<vJetSeed_iphi_.size()<<")"<<std::endl;
     std::cout<<" >> The jet seeds are (ieta,iphi): ";
     if (vJetSeed_ieta_.size()==0){vJetSeed_ieta_.push_back(-1); vJetSeed_iphi_.push_back(-1); std::cout<<"(-1, -1)"<<std::endl;}
     else{
     	for (int idx=0;idx<int(vJetSeed_ieta_.size());idx++){
     		std::cout<<"("<<vJetSeed_ieta_[idx]<<","<<vJetSeed_iphi_[idx]<<") ";
     	}
     	std::cout<<std::endl;
     }
     if (vJetSeed_ieta_.size()==vJetSeed_iphi_.size()){
     	for (int idx=0;idx<int(vJetSeed_ieta_.size());idx++){
     		if(vJetSeed_ieta_[idx]>=0){vJetSeed_ieta_[idx]=int(vJetSeed_ieta_[idx]*5+2);}  //5 EB xtals per HB tower
		if(vJetSeed_iphi_[idx]>=0){vJetSeed_iphi_[idx]=int(vJetSeed_iphi_[idx]*5+2);}  //5 EB xtals per HB tower
		//std::cout<<vJetSeed_ieta_[idx]<<" "<<vJetSeed_iphi_[idx];
     	}
     }
     std::unique_ptr<std::vector<int>> JetSeedieta_edm (new std::vector<int>(vJetSeed_ieta_));
     std::unique_ptr<std::vector<int>> JetSeediphi_edm (new std::vector<int>(vJetSeed_iphi_));
     iEvent.put(std::move(JetSeedieta_edm),"JetSeedieta");
     iEvent.put(std::move(JetSeediphi_edm),"JetSeediphi");
     vJetSeed_ieta_.clear(); vJetSeed_iphi_.clear();
   } else {
     std::cout<<" >> doJets not set"<<std::endl;
     passedSelection = runEvtSel( iEvent, iSetup );
     std::cout<<" >> Size of JetSeed vector (JetSeed_eta_size, JetSeed_phi_size) is: ("<<vJetSeed_ieta_.size()<<", "<<vJetSeed_iphi_.size()<<")"<<std::endl;
     std::cout<<" The jet seeds are (ieta,iphi): ";
     if (vJetSeed_ieta_.size()==0){vJetSeed_ieta_.push_back(-1); vJetSeed_iphi_.push_back(-1); std::cout<<"(-1, -1)"<<std::endl;}
     else{
	   for (int idx=0;idx<int(vJetSeed_ieta_.size());idx++){
     		std::cout<<" The jet seeds are (ieta,iphi): "<<"("<<vJetSeed_ieta_[idx]<<","<<vJetSeed_iphi_[idx]<<") ";
     	}
     	std::cout<<std::endl;
     }
     if (vJetSeed_ieta_.size()==vJetSeed_iphi_.size()){
     	for (int idx=0;idx<int(vJetSeed_ieta_.size());idx++){
     		if(vJetSeed_ieta_[idx]>=0){vJetSeed_ieta_[idx]=int(vJetSeed_ieta_[idx]*5+2);}  //5 EB xtals per HB tower
		if(vJetSeed_iphi_[idx]>=0){vJetSeed_iphi_[idx]=int(vJetSeed_iphi_[idx]*5+2);}  //5 EB xtals per HB tower
		//std::cout<<vJetSeed_ieta_[idx]<<" "<<vJetSeed_iphi_[idx];
     	}
     }
     std::unique_ptr<std::vector<int>> JetSeedieta_edm (new std::vector<int>(vJetSeed_ieta_));
     std::unique_ptr<std::vector<int>> JetSeediphi_edm (new std::vector<int>(vJetSeed_iphi_));
     iEvent.put(std::move(JetSeedieta_edm),"JetSeedieta");
     iEvent.put(std::move(JetSeediphi_edm),"JetSeediphi");
     vJetSeed_ieta_.clear(); vJetSeed_iphi_.clear();
   }

   if ( !passedSelection ) {
     h_sel->Fill( 0. );;
     return;
   }*/
   //auto photon_classes = std::make_unique<float>(10.0);
   fillEB( iEvent, iSetup );
   std::unique_ptr<std::vector<float>> EBenergy_edm (new std::vector<float>(vEB_energy_));
   /*for (unsigned int i=0;i<vEB_energy_.size();i++){
    std::cout<<"( "<<i/vEB_energy_width<<", "<<i%vEB_energy_width<<" ) = "<<vEB_energy_[i]<<" ";
   }*/
   std::cout<<" >> Adding EB done "<<std::endl;
   //EBEnergy_edm->push_back(vEB_energy_);
   //std::cout<<"Size1 is: "<<vEB_energy_.size()<<std::endl;
   std::cout<<" >> Size of EB Energy vector is: "<<std::move(EBenergy_edm).get()->size()<<std::endl;
   // PhotonCollection 
   //*photon_classes=get_photons(iEvent, iSetup );
   //iEvent.put(std::move(photon_classes),"photonClasses");
   iEvent.put(std::move(EBenergy_edm),"EBenergy");
 
   fillHBHE (iEvent, iSetup );
   std::unique_ptr<std::vector<float>> HBHEenergy_edm (new std::vector<float>(vHBHE_energy_));
   std::unique_ptr<std::vector<float>> HBHEenergyEB_edm (new std::vector<float>(vHBHE_energy_EB_));
   std::cout<<" >> Size of HBHE Energy vector is: "<<std::move(HBHEenergy_edm).get()->size()<<std::endl;
   std::cout<<" >> Size of EB HBHE Energy vector is: "<<std::move(HBHEenergyEB_edm).get()->size()<<std::endl;
   iEvent.put(std::move(HBHEenergy_edm),"HBHEenergy");
   iEvent.put(std::move(HBHEenergyEB_edm),"HBHEenergyEB");
   
   fillECALstitched (iEvent, iSetup);
   std::unique_ptr<std::vector<float>> ECALstitched_energy_edm (new std::vector<float>(vECAL_energy_));
   std::cout<<" >> Size of Stitched ECAL Energy vector is: "<<std::move(ECALstitched_energy_edm).get()->size()<<std::endl;
   iEvent.put(std::move(ECALstitched_energy_edm), "ECALstitchedenergy");

   
   fillTracksAtECALstitched (iEvent, iSetup );
   std::unique_ptr<std::vector<float>> TracksECALstitchedPt_edm (new std::vector<float>(vECAL_tracksPt_));
   std::cout<<" >> Size of Pt Tracks vector at Stitched ECAL is: "<<std::move(TracksECALstitchedPt_edm).get()->size()<<std::endl;
   iEvent.put(std::move(TracksECALstitchedPt_edm), "TracksAtECALstitched");
 
   std::cout<<" >> Added EB, HBHE, HBHE_EB, ECALstitched, Tracks_at_ECALstitched to edm root file"<<std::endl;
   //EBEnergy_edm->clear();
   //iEvent.put(photon_classes,"photon_classes");
   // Fill RHTree
   RHTree->Fill();
   //vEB_photon_frames.clear();
   //h_sel->Fill( 1. );
   nPassed++;
   /*for (int frame_x=0; frame_x<vEB_frame_height;frame_x++){
    for (int frame_y=0;frame_y<vEB_frame_width;frame_y++){
     std::cout<<"yes "<<
    }
   }*/
   
   std::cout<<std::endl;
   //predict_tf();
   //std::cout<<"TF_predict done "<<std::endl;
   //int vec_size=61200;
   //std::vector<float> vEB_energy;
   //vEB_energy=read_vEB_energy(vec_size);
   return;
}
/*void ProducerTest::predict_tf(){
 tensorflow::Session* session;
 tensorflow::GraphDef graph_def;
 tensorflow::SessionOptions opts;
 std::vector<tensorflow::Tensor> outputs; // Store outputs
 std::string graph_definition="\\home\\cmsusr\\CMSSW_10_6_8\\src\\ProdTutorial\\ProducerTest\\plugins\\graph3.pb";
 
 tensorflow::Tensor x(tensorflow::DT_FLOAT, tensorflow::TensorShape({100, 32}));
 tensorflow::Tensor y(tensorflow::DT_FLOAT, tensorflow::TensorShape({100, 8}));
 auto _XTensor = x.matrix<float>();
 auto _YTensor = y.matrix<float>();
 _XTensor.setRandom();
 _YTensor.setRandom();
 
 TF_CHECK_OK(ReadBinaryProto(Env::Default(), graph_definition, &graph_def));
 // load the graph definition, i.e. an object that contains the computational graph
 //tensorflow::GraphDef* graphDef = tensorflow::loadGraphDef("graph3.pb");
 // Set GPU options
 //graph::SetDefaultDevice("/gpu:0", &graph_def);
 //opts.config.mutable_gpu_options()->set_per_process_gpu_memory_fraction(0.5);
 //opts.config.mutable_gpu_options()->set_allow_growth(true);
 
 // create a new session
 //TF_CHECK_OK(NewSession(opts, &session));
 
 // Load graph into session
 TF_CHECK_OK(session->Create(graph_def));
 
 // create a session
 //session = tensorflow::createSession(graphDef);
 
 // Initialize our variables
 TF_CHECK_OK(session->Run({}, {}, {"init_all_vars_op"}, nullptr));
 //tensorflow::run(session, {}, {"init_all_vars_op"}, nullptr);
 
 //for (int i = 0; i < 10; ++i) {
        
 TF_CHECK_OK(session->Run({{"x", x}, {"y", y}}, {"cost"}, {}, &outputs)); // Get cost
 //tensorflow::run(session, { { "x", x }, {"y", y} }, { "cost" }, &outputs);
 float cost = outputs[0].scalar<float>()(0);
 std::cout << "Cost: " <<  cost << std::endl;
 //TF_CHECK_OK(session->Run({{"x", x}, {"y", y}}, {}, {"train"}, nullptr)); // Train
 //tensorflow::run(session, { { "x", x }, {"y", y} }, {}, {"train"}, &outputs);
 outputs.clear();
  
 session->Close();
 delete session;
 //std::cout<<_YTensor(0,0)<<" "<<_YTensor(0,1)<<" "<<_YTensor(0,2)<<" "<<_YTensor(0,3)<<" "<<_YTensor(0,4)<<" "<<_YTensor(0,5)<<" "<<_YTensor(0,6)<<" "<<_YTensor(0,7)<<" "<<_YTensor(0,8)<<" "<<_YTensor(0,9)<<endl;
 //std::cout<<_YTensor(1,0)<<" "<<_YTensor(1,1)<<" "<<_YTensor(1,2)<<" "<<_YTensor(1,3)<<" "<<_YTensor(1,4)<<" "<<_YTensor(1,5)<<" "<<_YTensor(1,6)<<" "<<_YTensor(1,7)<<" "<<_YTensor(1,8)<<" "<<_YTensor(1,9)<<endl;
 std::cout<<"All done"<<endl;
 // cleanup
 //tensorflow::closeSession(session);
 //delete graphDef;
}*/

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
DetImgProducer::beginStream(edm::StreamID)
{
 nTotal = 0;
 nPassed = 0;
 std::cout<<"'DetImgProducer' Stream began"<<std::endl;
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
DetImgProducer::endStream() {
 std::cout << "'DetImgProducer' selected: " << nPassed << "/" << nTotal << std::endl;
 //fw->Write();
 //fw->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void
ProducerTest::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ProducerTest::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ProducerTest::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ProducerTest::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DetImgProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

/*const reco::PFCandidate*
DetImgProducer::getPFCand(edm::Handle<PFCollection> pfCands, float eta, float phi, float& minDr, bool debug ) {

  minDr = 10;
  const reco::PFCandidate* minDRCand = nullptr;
  
  for ( PFCollection::const_iterator iPFC = pfCands->begin();
        iPFC != pfCands->end(); ++iPFC ) {

    const reco::Track* thisTrk = iPFC->bestTrack();
    if ( !thisTrk ) continue;

    float thisdR = reco::deltaR( eta, phi, thisTrk->eta(), thisTrk->phi() );
    if (debug) std::cout << "\tthisdR: " << thisdR << " " << thisTrk->pt() << " " << iPFC->particleId() << std::endl;

    const reco::PFCandidate& thisPFCand = (*iPFC);
      
    if ( (thisdR < 0.01) && (thisdR <minDr) ) {
      minDr    = thisdR; 
      minDRCand = &thisPFCand;
    }
  }

  return minDRCand;  
}

const reco::Track*
DetImgProducer::getTrackCand(edm::Handle<reco::TrackCollection> trackCands, float eta, float phi, float& minDr, bool debug ) {

  minDr = 10;
  const reco::Track* minDRCand = nullptr;
  reco::Track::TrackQuality tkQt_ = reco::Track::qualityByName("highPurity");

  for ( reco::TrackCollection::const_iterator iTk = trackCands->begin();
        iTk != trackCands->end(); ++iTk ) {
    if ( !(iTk->quality(tkQt_)) ) continue;  

    float thisdR = reco::deltaR( eta, phi, iTk->eta(),iTk->phi() );
    if (debug) std::cout << "\tthisdR: " << thisdR << " " << iTk->pt() << std::endl;

    const reco::Track& thisTrackCand = (*iTk);
      
    if ( (thisdR < 0.01) && (thisdR <minDr) ) {
      minDr    = thisdR; 
      minDRCand = &thisTrackCand;
    }
  }

  return minDRCand;  
}




int DetImgProducer::getTruthLabel(const reco::PFJetRef& recJet, edm::Handle<reco::GenParticleCollection> genParticles, float dRMatch , bool debug ){
  if ( debug ) {
    std::cout << " Mathcing reco jetPt:" << recJet->pt() << " jetEta:" << recJet->eta() << " jetPhi:" << recJet->phi() << std::endl;
  }

  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
       iGen != genParticles->end();
       ++iGen) {

    // From: (page 7/ Table 1.5.2)
    //https://indico.desy.de/indico/event/7142/session/9/contribution/31/material/slides/6.pdf
    //code range explanation:
    // 11 - 19 beam particles
    // 21 - 29 particles of the hardest subprocess
    // 31 - 39 particles of subsequent subprocesses in multiparton interactions
    // 41 - 49 particles produced by initial-state-showers
    // 51 - 59 particles produced by final-state-showers
    // 61 - 69 particles produced by beam-remnant treatment
    // 71 - 79 partons in preparation of hadronization process
    // 81 - 89 primary hadrons produced by hadronization process
    // 91 - 99 particles produced in decay process, or by Bose-Einstein effects

    // Do not want to match to the final particles in the shower
    if ( iGen->status() > 99 ) continue;
    
    // Only want to match to partons/leptons/bosons
    if ( iGen->pdgId() > 25 ) continue;

    float dR = reco::deltaR( recJet->eta(),recJet->phi(), iGen->eta(),iGen->phi() );

    if ( debug ) std::cout << " \t >> dR " << dR << " id:" << iGen->pdgId() << " status:" << iGen->status() << " nDaught:" << iGen->numberOfDaughters() << " pt:"<< iGen->pt() << " eta:" <<iGen->eta() << " phi:" <<iGen->phi() << " nMoms:" <<iGen->numberOfMothers()<< std::endl;

    if ( dR > dRMatch ) continue; 
    if ( debug ) std::cout << " Matched pdgID " << iGen->pdgId() << std::endl;

    return iGen->pdgId();

  } // gen particles 





  return -99;
}


float DetImgProducer::getBTaggingValue(const reco::PFJetRef& recJet, edm::Handle<edm::View<reco::Jet> >& recoJetCollection, edm::Handle<reco::JetTagCollection>& btagCollection, float dRMatch, bool debug ){

  // loop over jets
  for( edm::View<reco::Jet>::const_iterator jetToMatch = recoJetCollection->begin(); jetToMatch != recoJetCollection->end(); ++jetToMatch )
    {
      reco::Jet thisJet = *jetToMatch;
      float dR = reco::deltaR( recJet->eta(),recJet->phi(), thisJet.eta(),thisJet.phi() );
      if(dR > 0.1) continue;

      size_t idx = (jetToMatch - recoJetCollection->begin());
      edm::RefToBase<reco::Jet> jetRef = recoJetCollection->refAt(idx);

      if(debug) std::cout << "btag discriminator value = " << (*btagCollection)[jetRef] << std::endl;
      return (*btagCollection)[jetRef];
  
    }

  if(debug){
    std::cout << "ERROR  No btag match: " << std::endl;
    
    // loop over jets
    for( edm::View<reco::Jet>::const_iterator jetToMatch = recoJetCollection->begin(); jetToMatch != recoJetCollection->end(); ++jetToMatch )
      {
	const reco::Jet thisJet = *jetToMatch;
	std::cout << "\t Match attempt pt: " <<  thisJet.pt() << " vs " <<  recJet->pt()
		  << " eta: " << thisJet.eta() << " vs " << recJet->eta()
		  << "phi: "<< thisJet.phi() << " vs " << recJet->phi()
		  << std::endl;
	float dR = reco::deltaR( recJet->eta(),recJet->phi(), thisJet.eta(),thisJet.phi() );
	std::cout << "dR " << dR << std::endl;
      }
  }    

  return -99;
}*/



/*
//____ Fill FC diphoton variables _____//
void RecHitAnalyzer::fillFC ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  vFC_inputs_.clear();
  int ptOrder[2] = {0, 1};
  if ( vPho_[1].Pt() > vPho_[0].Pt() ) {
      ptOrder[0] = 1;
      ptOrder[1] = 0;
  }
  for ( int i = 0; i < 2; i++ ) {
    vFC_inputs_.push_back( vPho_[ptOrder[i]].Pt()/m0_ );
    vFC_inputs_.push_back( vPho_[ptOrder[i]].Eta() );
  }
  vFC_inputs_.push_back( TMath::Cos(vPho_[0].Phi()-vPho_[1].Phi()) );
} // fillFC() 
*/

//define this as a plug-in

//define this as a plug-in
DEFINE_FWK_MODULE(DetImgProducer);
