#include "ProdTutorial/ProducerTest/plugins/ProducerInference.h"
#include "ProdTutorial/ProducerTest/plugins/DetImgProducer.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include <fstream>

using namespace tensorflow;
using namespace std;

ProducerInference::ProducerInference(const edm::ParameterSet& iConfig)
{
 //vEB_photon_frames = consumes<std::vector<std::vector<float>>>(iConfig.getParameter<edm::InputTag>("frames_"));
 EBRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
 photonCollectionT_ = consumes<PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonCollection"));
 vEB_energy_token = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("EBEnergy"));
 ECALstitched_energy_token=consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("ECALstitchedenergy"));
 TracksAtECALstitched_token=consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("TracksAtECALstitched"));
 JetSeed_ieta_token=consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("JetSeedieta"));
 JetSeed_iphi_token=consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("JetSeediphi"));
 HBHEenergy_token = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("HBHEenergy"));
 std::cout<<"Reading data collection done "<<nTotal<<std::endl;
 
 produces<std::vector<float>>("EBenergyClass");
 produces<std::vector<int>>("ECALstitchedClass");
 produces<std::vector<int>>("TracksAtECALstitchedClass");
 produces<std::vector<int>>("HBHEenergyClass");
}

ProducerInference::~ProducerInference()
{
 
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
ProducerInference::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   nTotal++;

   edm::Handle<std::vector<float>> vEB_energy_handle;
   iEvent.getByToken(vEB_energy_token, vEB_energy_handle);
   edm::Handle<std::vector<float>> ECALstitched_energy_handle;
   iEvent.getByToken(ECALstitched_energy_token, ECALstitched_energy_handle);
   edm::Handle<std::vector<float>> TracksAtECALstitched_handle;
   iEvent.getByToken(TracksAtECALstitched_token, TracksAtECALstitched_handle);
   edm::Handle<std::vector<int>> JetSeed_ieta_handle;
   iEvent.getByToken(JetSeed_ieta_token, JetSeed_ieta_handle);
   edm::Handle<std::vector<int>> JetSeed_iphi_handle;
   iEvent.getByToken(JetSeed_iphi_token, JetSeed_iphi_handle);
   edm::Handle<std::vector<float>> HBHEenergy_handle;
   iEvent.getByToken(HBHEenergy_token, HBHEenergy_handle);
   
   std::vector<float>vECALstitched=*ECALstitched_energy_handle;
   std::vector<float>vTracksAtECALstitched=*TracksAtECALstitched_handle;
   std::vector<int>vJetSeed_ieta=*JetSeed_ieta_handle;
   std::vector<int>vJetSeed_iphi=*JetSeed_iphi_handle;
   std::vector<int>vECALstitchedClass;
   std::vector<int>vTracksAtECALstitchedClass;
   std::vector<int>vHBHEenergyClass;
   vECALstitchedClass.clear();
   vTracksAtECALstitchedClass.clear();
   vHBHEenergyClass.clear();
   for (int idx=0;idx<int(vJetSeed_ieta.size());idx++){
    std::cout<<" >> Generating Stitched ECAL frames and their track frames from the jet seed "<<idx+1<<"/"<<vJetSeed_ieta.size()<<" with seed value: ("<<vJetSeed_ieta[idx]<<","<<vJetSeed_iphi[idx]<<")"<<std::endl;
    if(vJetSeed_ieta[idx]>=0) {vECALstitched_frame=croppingFrames(vECALstitched, vJetSeed_ieta[idx], vJetSeed_iphi[idx], 280, 360, 125, 125); 
                               vTracksAtECALstitched_frame=croppingFrames(vTracksAtECALstitched, vJetSeed_ieta[idx], vJetSeed_iphi[idx], 280, 360, 125, 125);
    /*string filename="ECALstitched_"+std::to_string(nPassed+1)+"_"+std::to_string(idx+1)+".csv";
    std::ofstream file1(filename);
    for (int i=0;i<int(vECALstitched_frame.size());i++){
     for(int j=0;j<int(vECALstitched_frame[0].size());j++){
      file1<<vECALstitched_frame[i][j];
       if (j!=int(vECALstitched_frame[0].size())-1){file1<<",";}
      }
     file1<<"\n";
    }*/
    /*filename="TracksAtECALstitched_"+std::to_string(nPassed+1)+"_"+std::to_string(idx+1)+".csv";
    std::ofstream file2(filename);
    for (int i=0;i<int(vTracksAtECALstitched_frame.size());i++){
     for(int j=0;j<int(vTracksAtECALstitched_frame[0].size());j++){
       file2<<vTracksAtECALstitched_frame[i][j];
       if (j!=int(vTracksAtECALstitched_frame[0].size())-1){file2<<",";}
      }
     file2<<"\n";
    }*/
    vECALstitchedClass.push_back(predict_tf(vECALstitched_frame, "qg_model.pb", "inputs","softmax_1/Sigmoid"));
    vTracksAtECALstitchedClass.push_back(predict_tf(vTracksAtECALstitched_frame, "qg_model.pb", "inputs", "softmax_1/Sigmoid"));
    }
    else {
     vECALstitchedClass.push_back(-1);
     vTracksAtECALstitchedClass.push_back(-1);
    }
    std::cout<<" >> Predicted Class of Stitched ECAL: "<<vECALstitchedClass[idx]<<std::endl;
    std::cout<<" >> Predicted Class of Tracks at Stitched ECAL: "<<vTracksAtECALstitchedClass[idx]<<std::endl;
   }
   //std::cout<<std::endl; //Stitched ECAL and their track frames created.
   
  
   std::vector<float> vHBHEenergy=*HBHEenergy_handle;
   std::cout<<"Size of HBHE energy vector read: "<<vHBHEenergy.size()<<std::endl;
   std::vector<std::vector<float>> vHBHEenergy_strided = frameStriding(vHBHEenergy,56,72,5,5);
   std::cout<<"Size of Strided HBHE energy vector: ("<<vHBHEenergy_strided.size()<<","<<vHBHEenergy_strided[0].size()<<")"<<std::endl; //HBHE energy vector upsampled.
   std::vector<float> vHBHE_strided_flat (vHBHEenergy_strided.size()*vHBHEenergy_strided[0].size(),0);
   for (int x=0;x<int(vHBHEenergy_strided.size());x++){
    for (int y=0;y<int(vHBHEenergy_strided[0].size());y++){
     vHBHE_strided_flat[x*vHBHEenergy_strided[0].size()+y]=vHBHEenergy_strided[x][y];
    }
   }
   std::vector<std::vector<float>> vHBHEenergy_frame;
   for (int idx=0;idx<int(vJetSeed_ieta.size());idx++){
    std::cout<<" >> Generating HBHE energy frames from the jet seed "<<idx+1<<"/"<<vJetSeed_ieta.size()<<" with seed value: ("<<vJetSeed_ieta[idx]<<","<<vJetSeed_iphi[idx]<<")"<<std::endl;
    if(vJetSeed_ieta[idx]>=0) {vHBHEenergy_frame=croppingFrames(vHBHE_strided_flat, vJetSeed_ieta[idx], vJetSeed_iphi[idx], 280, 360, 125, 125); 
   /*for(int i=140;i<141;i++){
    for (int ki=0; ki<5;ki++){
     for (int j=0;j<360;j++){
      for (int kj=0;kj<5;kj++){
      std::cout<<"("<<i<<","<<j<<"): "<<vHBHEenergy_strided[5*i+ki][5*j+kj]<<" "<<vHBHEenergy[i*360+j]/25<<" ";
     }
     }
    }
   }
   std::cout<<std::endl;*/
    //std::cout<<std::endl;
    /*string filename="HBHEenergy"+std::to_string(nPassed+1)+".csv";
    std::ofstream file3(filename);
    for (int i=0;i<int(vHBHEenergy_frame.size());i++){
     for(int j=0;j<int(vHBHEenergy_frame[0].size());j++){
      file3<<vHBHEenergy_frame[i][j];
      if (j!=int(vHBHEenergy_frame[0].size())-1){file3<<",";}
     }
    file3<<"\n";
   }*/
   vHBHEenergyClass.push_back(predict_tf(vHBHEenergy_frame, "qg_model.pb", "inputs", "softmax_1/Sigmoid"));
   }
   else {vHBHEenergyClass.push_back(-1);}
   std::cout<<" >> Predicted Class of HBHE energy: "<<vHBHEenergyClass[idx]<<std::endl;
   }
   std::unique_ptr<std::vector<int>> vECALstitchedClass_edm (new std::vector<int>(vECALstitchedClass));
   iEvent.put(std::move(vECALstitchedClass_edm),"ECALstitchedClass");
   std::unique_ptr<std::vector<int>> vTracksAtECALstitchedClass_edm (new std::vector<int>(vTracksAtECALstitchedClass));
   iEvent.put(std::move(vTracksAtECALstitchedClass_edm),"TracksAtECALstitchedClass");
   std::unique_ptr<std::vector<int>> vHBHEenergyClass_edm (new std::vector<int>(vHBHEenergyClass));
   iEvent.put(std::move(vHBHEenergyClass_edm),"HBHEenergyClass");
 
   //std::cout<<"Size1: "<<vEB_energy_handle->size()<<std::endl;
   vEB_energy_=*vEB_energy_handle;
   //std::cout<<"Size2: "<<vEB_energy_.size();
   get_photons(iEvent, iSetup );//stored in vEB_frames vectors
   std::unique_ptr<std::vector<float>> vclasses_edm (new std::vector<float>(vclasses));
   iEvent.put(std::move(vclasses_edm),"EBenergyClass");
   std::cout<<std::endl;
   nPassed++;
   // ----- Apply event selection cuts ----- //
   //std::cout<<"Event "<<nPassed-1<<"finished"<<"Prceeding to event "<<nPassed<<std::endl;
   return;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ProducerInference::beginStream(edm::StreamID)
{
 nTotal = 0;
 nPassed = 0;
 std::cout<<"'ProducerInference' Stream began"<<std::endl;
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ProducerInference::endStream() {
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
ProducerInference::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ProducerInference);
