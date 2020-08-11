#ifndef QGProducer_h
#define QGProducer_h

#include <memory>
//#include <iostream>
//#include <fstream>
//#include <sstream>
#include <vector>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
//#include "FWCore/Framework/interface/one/EDAnalyzer.h"
//#include "FWCore/Framework/src/one/implementorsMethods.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"
#include "DQM/HcalCommon/interface/Constants.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" // reco::PhotonCollection defined here
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TFrame.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
/*#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"*/
#include <iostream>
using namespace std;
/*using pat::PhotonCollection;
using pat::PhotonRef;*/
using reco::PhotonCollection;
using reco::PhotonRef;

class QGProducer : public edm::stream::EDProducer<> {
   public:
      
      explicit QGProducer(const edm::ParameterSet&);
      ~QGProducer();
      
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
      // Tokens 
      //edm::EDGetTokenT<PhotonCollection> photonCollectionT_;
      //edm::EDGetTokenT<std::vector<std::vector<float>>> frames_;
      edm::EDGetTokenT<PhotonCollection> photonCollectionT_;
      edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_;
      edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitCollectionT_;
      edm::EDGetTokenT<std::vector<float>> ECALstitched_energy_token;
      edm::EDGetTokenT<std::vector<float>> TracksAtECALstitched_token;
      //edm::EDGetTokenT<std::vector<int>> JetSeed_ieta_token;
      //edm::EDGetTokenT<std::vector<int>> JetSeed_iphi_token;
      edm::EDGetTokenT<std::vector<float>> HBHEenergy_token;
      edm::EDGetTokenT<reco::VertexCollection> vertexCollectionT_;
      edm::EDGetTokenT<reco::PFJetCollection> jetCollectionT_;
      edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionT_;
      edm::EDGetTokenT<reco::JetTagCollection> jetTagCollectionT_;
      edm::EDGetTokenT<std::vector<reco::CandIPTagInfo> >    ipTagInfoCollectionT_;
      edm::EDGetTokenT<reco::GenJetCollection> genJetCollectionT_;
      edm::EDGetTokenT<TrackingRecHitCollection> TRKRecHitCollectionT_;
      edm::EDGetTokenT<edm::View<reco::Jet> > recoJetsT_;
      edm::EDGetTokenT<EBDigiCollection>     EBDigiCollectionT_;
      static const int nPhotons = 2;
   
      TH1F *h_sel;
      TTree* QGTree;
      
      void branchesEvtSel         ( TTree*, edm::Service<TFileService>& );
      void branchesEvtSel_jet     ( TTree*, edm::Service<TFileService>& );
      void branchesEvtSel_jet_dijet      ( TTree*, edm::Service<TFileService>& );
      void branchesEvtSel_jet_dijet_gg_qq( TTree*, edm::Service<TFileService>& );
      bool runEvtSel_jet_dijet      ( const edm::Event&, const edm::EventSetup& );
      bool runEvtSel_jet_dijet_gg_qq( const edm::Event&, const edm::EventSetup& );
      void fillEvtSel_jet_dijet      ( const edm::Event&, const edm::EventSetup& );
      void fillEvtSel_jet_dijet_gg_qq( const edm::Event&, const edm::EventSetup& );
      bool runEvtSel          ( const edm::Event&, const edm::EventSetup& );
      bool runEventSel_jet       ( const edm::Event&, const edm::EventSetup& );
      bool runEvtSel_jet      ( const edm::Event&, const edm::EventSetup& );
      void branchesJetInfoAtECALstitched   ( TTree*, edm::Service<TFileService>& );
      
      typedef std::vector<reco::PFCandidate>  PFCollection;
      edm::EDGetTokenT<PFCollection> pfCollectionT_;
      //std::vector<std::vector<float>> vEB_frame; //= std::vector<std::vector<float>> (vEB_frame_height,std::vector<float> (vEB_frame_width, 0.0));
      //std::vector<float> vEB_flat_frame = std::vector<float> (vEB_frame_height*vEB_frame_width,0.0);
      std::vector<std::vector<float>> vHBHEenergy_frame;
      std::vector<std::vector<float>> vECALstitched_frame;
      std::vector<std::vector<float>> vTracksAtECALstitched_frame;
      std::vector<float> vclasses;
      vector<int> vJetSeed_iphi_;
      vector<int> vJetSeed_ieta_;
      std::vector<float> vSC_eta_;
      std::vector<float> vSC_phi_;
   
      unsigned int nPho;
     
      //std::vector<vector<float>> croppingFrames   (std::vector<float>&, int ,int, int, int, int, int);
      //std::vector<vector<float>> frameStriding    (std::vector<float>&, int, int, int, int);
      //int predict_tf                              (std::vector<std::vector<float>>&, string, string, string);
      std::string mode_;  // EventLevel / JetLevel
      std::vector<int> vJetIdxs;
      bool doJets_;
      int  nJets_;
      int iphi_Emax, ieta_Emax;
      double minJetPt_;
      double maxJetEta_;
      double z0PVCut_;
   
      const reco::PFCandidate* getPFCand(edm::Handle<PFCollection> pfCands, float eta, float phi, float& minDr, bool debug = false);
      const reco::Track* getTrackCand(edm::Handle<reco::TrackCollection> trackCands, float eta, float phi, float& minDr, bool debug = false);
      int   getTruthLabel(const reco::PFJetRef& recJet, edm::Handle<reco::GenParticleCollection> genParticles, float dRMatch = 0.4, bool debug = false);
      float getBTaggingValue(const reco::PFJetRef& recJet, edm::Handle<edm::View<reco::Jet> >& recoJetCollection, edm::Handle<reco::JetTagCollection>& btagCollection, float dRMatch = 0.1, bool debug= false );
      /*TTree* RHTree;
      unsigned int nPho;
      
      void branchesEB             ( TTree*, edm::Service<TFileService>& );
      void branchesPhotonSel      ( TTree*, edm::Service<TFileService>& );*/
      //void fill_photons             ( const edm::Event&, const edm::EventSetup& );
      
      
      /*std::vector<float> vIphi_Emax_;
      std::vector<float> vIeta_Emax_;
      std::vector<float> vSC_eta_;
      std::vector<float> vSC_phi_;
      std::vector<int> vPreselPhoIdxs_;*/
      int nTotal, nPassed;
};
//static const float zs = 0.;
/*static const int EB_IPHI_MIN = EBDetId::MIN_IPHI;//1;
static const int EB_IPHI_MAX = EBDetId::MAX_IPHI;//360;
static const int EB_IETA_MIN = EBDetId::MIN_IETA;//1;
static const int EB_IETA_MAX = EBDetId::MAX_IETA;//85;
static const float zs = 0.;*/
#endif
//DEFINE_FWK_MODULE(QGProducer);
