#include "ProdTutorial/ProducerTest/plugins/DetImgProducer.h"
#include "ProdTutorial/ProducerTest/plugins/EGProducer.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
//#include <fstream>

EGProducer::EGProducer(const edm::ParameterSet& iConfig)
{
 //vEB_photon_frames = consumes<std::vector<std::vector<float>>>(iConfig.getParameter<edm::InputTag>("frames_"));
 EBRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
 photonCollectionT_ = consumes<PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonCollection"));
 vEB_energy_token = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("EBEnergy"));
 std::cout<<"Reading data collection done "<<nTotal<<std::endl;
 
 edm::Service<TFileService> fs;
 EGTree = fs->make<TTree>("EGTree", "RecHit tree");
 branchesPhotonSel ( EGTree, fs );
 produces<std::vector<float>>("EBenergyClass");
}
  
 EGProducer::~EGProducer()
{
 
}

void
EGProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   nTotal++;

   edm::Handle<std::vector<float>> vEB_energy_handle;
   iEvent.getByToken(vEB_energy_token, vEB_energy_handle);
  
   vEB_energy_=*vEB_energy_handle;
   vEB_flat_frame.clear();
   vEB_frame.clear();
   vclasses.clear();
   vEB_photon_frames.clear();
 
   get_photons(iEvent, iSetup );//stored in vEB_frames vectors
   std::unique_ptr<std::vector<float>> vclasses_edm (new std::vector<float>(vclasses));
   iEvent.put(std::move(vclasses_edm),"EBenergyClass");
   std::cout<<std::endl;
   nPassed++;
  // ----- Apply event selection cuts ----- //
   //std::cout<<"Event "<<nPassed-1<<"finished"<<"Prceeding to event "<<nPassed<<std::endl;
   return;
}

void
EGProducer::beginStream(edm::StreamID)
{
 nTotal = 0;
 nPassed = 0;
 std::cout<<"'EGProducer' Stream began"<<std::endl;
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
EGProducer::endStream() {
 std::cout << "'ProducerInference' selected: " << nPassed << "/" << nTotal << std::endl;
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
EGProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EGProducer);
