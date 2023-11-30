#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/Math/interface/libminifloat.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "CommonTools/Egamma/interface/EffectiveAreas.h"
#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "RecoSimStudies/Dumpers/plugins/RecoSimDumper.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHadTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsMCRcd.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAlphas.h"
#include "CondFormats/DataRecord/interface/EcalLaserAlphasRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatiosRef.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRefRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatios.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRcd.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"

#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TTree.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include <set>
#include <assert.h>
#include <time.h>

#include <TMath.h>
#include <Math/VectorUtil.h>
//#include <boost/tokenizer.hpp>

using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;

//
// constructors and destructor
//
RecoSimDumper::RecoSimDumper(const edm::ParameterSet& iConfig):
   caloTopologyToken_(esConsumes()),
   caloGeometryToken_(esConsumes()),
   ADCtoGeVToken_(esConsumes()),
   alphaToken_(esConsumes()),
   APDPNRatiosToken_(esConsumes()),
   icalToken_(esConsumes()),
   icalMCToken_(esConsumes()),
   channelStatusToken_(esConsumes()),
   pedsToken_(esConsumes()),
   ratioToken_(esConsumes())
{
   usesResource(TFileService::kSharedResource);

   pileupSummaryToken_            = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummary"));
   vtxToken_                      = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
   rhoToken_                      = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"));
   genToken_                      = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   caloPartToken_                 = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
   puCaloPartToken_               = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("puCaloParticleCollection"));
   ootpuCaloPartToken_            = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("ootpuCaloParticleCollection"));
   ebRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
   pfRecHitToken_                 = consumes<std::vector<reco::PFRecHit> >(iConfig.getParameter<edm::InputTag>("pfRechitCollection")); 
   pfClusterToken_                = consumes<std::vector<reco::PFCluster> >(iConfig.getParameter<edm::InputTag>("pfClusterCollection")); 
   ebSuperClusterToken_           = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebSuperClusterCollection"));
   eeSuperClusterToken_           = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeSuperClusterCollection"));
   ebRetunedSuperClusterToken_    = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebRetunedSuperClusterCollection"));
   eeRetunedSuperClusterToken_    = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeRetunedSuperClusterCollection"));
   ebDeepSuperClusterToken_       = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebDeepSuperClusterCollection"));
   eeDeepSuperClusterToken_       = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeDeepSuperClusterCollection"));
   gsfElectronToken_              = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("gsfElectronCollection"));
   gedPhotonToken_                = consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("gedPhotonCollection"));
   patElectronToken_              = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("patElectronCollection"));
   patPhotonToken_                = consumes<std::vector<pat::Photon> >(iConfig.getParameter<edm::InputTag>("patPhotonCollection")); 
   patJetToken_                   = consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("patJetCollection")); 
   patMETToken_                   = consumes<std::vector<pat::MET> >(iConfig.getParameter<edm::InputTag>("patMETCollection"));
   
   nBits_                         = iConfig.getParameter<int>("nBits");
   doCompression_                 = iConfig.getParameter<bool>("doCompression");
   isMC_                          = iConfig.getParameter<bool>("isMC"); 
   saveGenParticles_              = iConfig.getParameter<bool>("saveGenParticles");
   saveCaloParticles_             = iConfig.getParameter<bool>("saveCaloParticles");
   saveCaloParticlesPU_           = iConfig.getParameter<bool>("saveCaloParticlesPU");
   saveCaloParticlesOOTPU_        = iConfig.getParameter<bool>("saveCaloParticlesOOTPU");   
   subtractSignalCalo_            = iConfig.getParameter<bool>("subtractSignalCalo"); 
   saveSimhits_             	  = iConfig.getParameter<bool>("saveSimhits");
   saveSimhitsPU_             	  = iConfig.getParameter<bool>("saveSimhitsPU");
   saveRechits_                   = iConfig.getParameter<bool>("saveRechits");
   savePFRechits_                 = iConfig.getParameter<bool>("savePFRechits"); 
   savePFCluster_                 = iConfig.getParameter<bool>("savePFCluster");
   savePFClusterhits_             = iConfig.getParameter<bool>("savePFClusterhits");
   saveShowerShapes_              = iConfig.getParameter<bool>("saveShowerShapes"); 
   saveSuperCluster_              = iConfig.getParameter<bool>("saveSuperCluster"); 
   saveRetunedSC_                 = iConfig.getParameter<bool>("saveRetunedSC");  
   saveDeepSC_                    = iConfig.getParameter<bool>("saveDeepSC");  
   saveGsfElectrons_              = iConfig.getParameter<bool>("saveGsfElectrons");
   saveGedPhotons_                = iConfig.getParameter<bool>("saveGedPhotons");
   savePatPhotons_                = iConfig.getParameter<bool>("savePatPhotons"); 
   savePatElectrons_              = iConfig.getParameter<bool>("savePatElectrons");
   savePatJets_                   = iConfig.getParameter<bool>("savePatJets");   

   egmCutBasedElectronIDVeto_     = iConfig.getParameter<std::string>("egmCutBasedElectronIDVeto"); 
   egmCutBasedElectronIDloose_    = iConfig.getParameter<std::string>("egmCutBasedElectronIDloose");  
   egmCutBasedElectronIDmedium_   = iConfig.getParameter<std::string>("egmCutBasedElectronIDmedium"); 
   egmCutBasedElectronIDtight_    = iConfig.getParameter<std::string>("egmCutBasedElectronIDtight");   
   egmMVAElectronIDloose_         = iConfig.getParameter<std::string>("egmMVAElectronIDloose");  
   egmMVAElectronIDmedium_        = iConfig.getParameter<std::string>("egmMVAElectronIDmedium"); 
   egmMVAElectronIDtight_         = iConfig.getParameter<std::string>("egmMVAElectronIDtight");   
   egmMVAElectronIDlooseNoIso_    = iConfig.getParameter<std::string>("egmMVAElectronIDlooseNoIso");  
   egmMVAElectronIDmediumNoIso_   = iConfig.getParameter<std::string>("egmMVAElectronIDmediumNoIso"); 
   egmMVAElectronIDtightNoIso_    = iConfig.getParameter<std::string>("egmMVAElectronIDtightNoIso");
   heepElectronID_                = iConfig.getParameter<std::string>("heepElectronID");   
   egmCutBasedPhotonIDloose_      = iConfig.getParameter<std::string>("egmCutBasedPhotonIDloose");  
   egmCutBasedPhotonIDmedium_     = iConfig.getParameter<std::string>("egmCutBasedPhotonIDmedium"); 
   egmCutBasedPhotonIDtight_      = iConfig.getParameter<std::string>("egmCutBasedPhotonIDtight"); 
   egmMVAPhotonIDmedium_          = iConfig.getParameter<std::string>("egmMVAPhotonIDmedium"); 
   egmMVAPhotonIDtight_           = iConfig.getParameter<std::string>("egmMVAPhotonIDtight");
   jmeCutBasedPFJetIDloose_         = iConfig.getParameter<std::vector<std::string> >("jmeCutBasedPFJetIDloose");
   jmeCutBasedPFJetIDtight_         = iConfig.getParameter<std::vector<std::string> >("jmeCutBasedPFJetIDtight");
   jmeCutBasedPFJetIDtightLepVeto_  = iConfig.getParameter<std::vector<std::string> >("jmeCutBasedPFJetIDtightLepVeto");
  
   if(nBits_>23 && doCompression_){
      cout << "WARNING: float compression bits > 23 ---> Using 23 (i.e. no compression) instead!" << endl;
      nBits_=23;
   }

   //output file, historgrams and trees
   tree = iFile->make<TTree>("caloTree","caloTree"); 
   setTree(tree);
}

RecoSimDumper::~RecoSimDumper()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void RecoSimDumper::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
   //calo geometry
   edm::ESHandle<CaloGeometry> pCaloGeometry = iSetup.getHandle(caloGeometryToken_);
   geometry = pCaloGeometry.product();
   _ebGeom = pCaloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
   _eeGeom = pCaloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
   _esGeom = pCaloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
   if (_esGeom) {
    for (uint32_t ic = 0; ic < _esGeom->getValidDetIds().size() && (!_esPlus || !_esMinus); ++ic) {
      const double z = _esGeom->getGeometry(_esGeom->getValidDetIds()[ic])->getPosition().z();
      _esPlus = _esPlus || (0 < z);
      _esMinus = _esMinus || (0 > z);
    }
   }

   edm::ESHandle<CaloTopology> pCaloTopology = iSetup.getHandle(caloTopologyToken_); 
   topology = pCaloTopology.product();

   edm::ESHandle<EcalADCToGeVConstant> pADCtoGeV = iSetup.getHandle(ADCtoGeVToken_);
   adcToGeV = pADCtoGeV.product();

   edm::ESHandle<EcalLaserAlphas> pAlpha = iSetup.getHandle(alphaToken_); 
   laserAlpha = pAlpha.product();

   edm::ESHandle<EcalLaserAPDPNRatios> pAPDPNRatios = iSetup.getHandle(APDPNRatiosToken_); 
   laserRatio = pAPDPNRatios.product(); 

   edm::ESHandle<EcalIntercalibConstants> pIcal = iSetup.getHandle(icalToken_); 
   ical = pIcal.product();

   if(isMC_){ 
    edm::ESHandle<EcalIntercalibConstantsMC> pIcalMC = iSetup.getHandle(icalMCToken_); 
    icalMC = pIcalMC.product();
   }

   edm::ESHandle<EcalChannelStatus> pChannelStatus = iSetup.getHandle(channelStatusToken_); 
   chStatus = pChannelStatus.product();

   edm::ESHandle<EcalPedestals> pPeds = iSetup.getHandle(pedsToken_); 
   ped = pPeds.product();
 
   edm::ESHandle<EcalGainRatios> pRatio = iSetup.getHandle(ratioToken_);  
   gr = pRatio.product();

   //MC-only info and collections
   truePU=-1.;
   obsPU=-1.;
   edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
   edm::Handle<std::vector<reco::GenParticle> > genParticlesTot;
   edm::Handle<std::vector<reco::GenParticle> > genParticles; 
   edm::Handle<std::vector<CaloParticle> > caloParticles;
   edm::Handle<std::vector<CaloParticle> > puCaloParticle;
   edm::Handle<std::vector<CaloParticle> > ootpuCaloParticle;

   if(isMC_){    
    ev.getByToken(pileupSummaryToken_, PupInfo);
    if (PupInfo.isValid()) 
    {
       for(auto &pu : *PupInfo){
           if(pu.getBunchCrossing() == 0 ){
              truePU = pu.getTrueNumInteractions();
              obsPU = pu.getPU_NumInteractions();
              break;
           } 
       } 
    }else{
       std::cerr << "Analyze --> PupInfo not found" << std::endl;
    }

    ev.getByToken(genToken_,genParticles);
    if (!genParticles.isValid()) {
        std::cerr << "Analyze --> genParticles not found" << std::endl; 
        return;
    }
    
    if(saveCaloParticles_) {  
       ev.getByToken(caloPartToken_,caloParticles);
       if (!caloParticles.isValid() && isMC_) {
           std::cerr << "Analyze --> caloParticles not found" << std::endl; 
           return;
       }
    }

    if(saveCaloParticlesPU_) {  
      ev.getByToken(puCaloPartToken_,puCaloParticle);
      if (!puCaloParticle.isValid()) {
          std::cerr << "Analyze --> puCaloParticle not found" << std::endl; 
          return;
      }
    }

    if(saveCaloParticlesOOTPU_) {  
      ev.getByToken(ootpuCaloPartToken_,ootpuCaloParticle);
      if (!ootpuCaloParticle.isValid()) {
          std::cerr << "Analyze --> ootpuCaloParticle not found" << std::endl; 
          return;
      }
    }
   }

   //Other collections
   edm::Handle<double> rhos;
   ev.getByToken(rhoToken_,rhos);
   if (!rhos.isValid()) {
       std::cerr << "Analyze --> rhos not found" << std::endl; 
       return;
   }

   edm::Handle<reco::VertexCollection> vertices;
   ev.getByToken(vtxToken_,vertices);
   if (!vertices.isValid()) {
       std::cerr << "Analyze --> vertices not found" << std::endl; 
       return;
   }
    
   edm::Handle<EcalRecHitCollection> recHitsEB;
   ev.getByToken(ebRechitToken_, recHitsEB);
   if (!recHitsEB.isValid()) {
       std::cerr << "Analyze --> recHitsEB not found" << std::endl; 
       return;
   }

   edm::Handle<EcalRecHitCollection> recHitsEE;
   ev.getByToken(eeRechitToken_, recHitsEE);
   if (!recHitsEE.isValid()) {
       std::cerr << "Analyze --> recHitsEE not found" << std::endl; 
       return;
   } 

   edm::Handle<std::vector<reco::PFRecHit> > pfRecHits;
   if(savePFRechits_) {
      ev.getByToken(pfRecHitToken_, pfRecHits);
      if (!pfRecHits.isValid()) {
          std::cerr << "Analyze --> pfRecHits not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::PFCluster> > pfClusters;
   if(savePFCluster_) {
      ev.getByToken(pfClusterToken_, pfClusters);
      if (!pfClusters.isValid()) {
          std::cerr << "Analyze --> pfClusters not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEB;
   edm::Handle<std::vector<reco::SuperCluster> > superClusterEE;
   if(saveSuperCluster_) {
      ev.getByToken(ebSuperClusterToken_, superClusterEB);
      if (!superClusterEB.isValid()) {
          std::cerr << "Analyze --> superClusterEB not found" << std::endl; 
          return;
      }
      ev.getByToken(eeSuperClusterToken_, superClusterEE);
      if (!superClusterEE.isValid()) {
          std::cerr << "Analyze --> superClusterEE not found" << std::endl; 
          return;
      }
   }

   edm::Handle<std::vector<reco::SuperCluster> > retunedSuperClusterEB;
   edm::Handle<std::vector<reco::SuperCluster> > retunedSuperClusterEE;
   if(saveRetunedSC_) {
      ev.getByToken(ebRetunedSuperClusterToken_, retunedSuperClusterEB);
      if (!retunedSuperClusterEB.isValid()) {
          std::cerr << "Analyze --> retunedSuperClusterEB not found" << std::endl; 
          return;
      }
      ev.getByToken(eeRetunedSuperClusterToken_, retunedSuperClusterEE);
      if (!retunedSuperClusterEE.isValid()) {
          std::cerr << "Analyze --> retunedSuperClusterEE not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > deepSuperClusterEB;
   edm::Handle<std::vector<reco::SuperCluster> > deepSuperClusterEE;
   if(saveDeepSC_) {
      ev.getByToken(ebDeepSuperClusterToken_, deepSuperClusterEB);
      if (!deepSuperClusterEB.isValid()) {
          std::cerr << "Analyze --> deepSuperClusterEB not found" << std::endl; 
          return;
      }
      ev.getByToken(eeDeepSuperClusterToken_, deepSuperClusterEE);
      if (!deepSuperClusterEE.isValid()) {
          std::cerr << "Analyze --> deepSuperClusterEE not found" << std::endl; 
          return;
      }
   }

   edm::Handle<std::vector<reco::GsfElectron> > gsfElectron;
   edm::Handle<std::vector<reco::Photon> > gedPhoton;
   if(saveGsfElectrons_){
      ev.getByToken(gsfElectronToken_, gsfElectron);
      if (!gsfElectron.isValid()) {
          std::cerr << "Analyze --> gsfElectron not found" << std::endl; 
          return;
      }
   }
   if(saveGedPhotons_){
      ev.getByToken(gedPhotonToken_, gedPhoton);
      if (!gedPhoton.isValid()) {
          std::cerr << "Analyze --> gedPhoton not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<pat::Electron> > patElectron;
   edm::Handle<std::vector<pat::Photon> > patPhoton;
   edm::Handle<std::vector<pat::Jet> > patJet;
   edm::Handle< std::vector<pat::MET> > patMET;
   if(savePatPhotons_) {
      ev.getByToken(patPhotonToken_,patPhoton);
      if (!patPhoton.isValid()) {
          std::cerr << "Analyze --> patPhotons not found" << std::endl; 
          return;
      }
   }
   if(savePatElectrons_) {
      ev.getByToken(patElectronToken_,patElectron);
      if (!patElectron.isValid()) {
          std::cerr << "Analyze --> patElectrons not found" << std::endl; 
          return;
      }
   }
   if(savePatJets_) {
      ev.getByToken(patJetToken_,patJet);
      if (!patJet.isValid()) {
          std::cerr << "Analyze --> patJets not found" << std::endl; 
          return;
      }
   }
   if(savePatPhotons_ || savePatElectrons_ || savePatJets_) {
      ev.getByToken(patMETToken_,patMET);
      if (!patMET.isValid()) {
          std::cerr << "Analyze --> patMET not found" << std::endl; 
          return;
      }
   }

   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();
   nVtx = vertices->size();
   rho = *(rhos.product());

   genParticle_genMotherIndex.clear();
   genParticle_pdgId.clear();
   genParticle_status.clear();
   genParticle_statusFlag.clear();
   genParticle_energy.clear();
   genParticle_pt.clear();
   genParticle_eta.clear();
   genParticle_phi.clear();

   genParts.clear();
   int nGenParticles = 0;
   
   if(isMC_){ 
     
    const std::vector<reco::GenParticle>& genParts_total = *(genParticles.product());
    genParts = genParts_total;
    nGenParticles = genParts.size(); 
  
    genParticle_genMotherIndex.resize(nGenParticles);
    genParticle_pdgId.resize(nGenParticles);
    genParticle_status.resize(nGenParticles);
    genParticle_statusFlag.resize(nGenParticles);
    genParticle_energy.resize(nGenParticles);
    genParticle_pt.resize(nGenParticles);
    genParticle_eta.resize(nGenParticles);
    genParticle_phi.resize(nGenParticles);

    // loop on all the selected genParticles 
    for(unsigned int iGen=0; iGen<genParts.size(); iGen++)
    {
       genParticle_genMotherIndex[iGen] = getGenMother(&genParts.at(iGen)); 
       genParticle_pdgId[iGen] = genParts.at(iGen).pdgId(); 
       genParticle_status[iGen] = genParts.at(iGen).status(); 
       genParticle_statusFlag[iGen] = getGenStatusFlag(&genParts.at(iGen)); 
       genParticle_energy[iGen] = genParts.at(iGen).energy(); 
       genParticle_pt[iGen] = genParts.at(iGen).pt();
       genParticle_eta[iGen] = genParts.at(iGen).eta();
       genParticle_phi[iGen] = genParts.at(iGen).phi();
    } 
    
    genParticle_size = nGenParticles; 
    //std::cout << "GenParticles size  : " << nGenParticles << std::endl;
   }

   GlobalPoint cell;
   hitsAndEnergies_CaloPart.clear();
   hitsAndEnergies_CaloPartPU.clear();
   hitsAndEnergies_CaloPartOOTPU.clear(); 
   
   std::vector<CaloParticle> caloParts;
   std::vector<GlobalPoint> caloParts_position;
   caloParticle_size = 0;
   caloParts.clear();

   if(saveCaloParticles_ && isMC_)
   { 
    for(const auto& iCalo : *(caloParticles.product()))
    {
       caloParticle_size++;
       std::vector<std::pair<DetId, float> > caloParticle_hitsAndEnergies = *getHitsAndEnergiesCaloPart(&iCalo,-1.);
       GlobalPoint caloParticle_position = calculateAndSetPositionActual(&caloParticle_hitsAndEnergies, 7.4, 3.1, 1.2, 4.2, 0.89, 0.,true);
       if(caloParticle_position == GlobalPoint(-999999., -999999., -999999.)){
          std::cout << "Invalid position for caloParticle, skipping caloParticle!" << std::endl;
          continue;
       }    

       hitsAndEnergies_CaloPart.push_back(caloParticle_hitsAndEnergies);
       caloParts_position.push_back(caloParticle_position);
       caloParts.push_back(iCalo); 
    }
   }
 
   int nCaloParticles = caloParts.size(); 
   //std::cout << "CaloParticles size  : " << nCaloParticles << " - " << caloParticle_size << " - " << caloParticlePU_size << " - " << caloParticleOOTPU_size << " - " << nVtx << std::endl;

   if(saveCaloParticlesPU_ && isMC_)
   {
      caloParticlePU_xtalEnergy.clear();
      caloParticlePU_xtalEta.clear();
      caloParticlePU_xtalPhi.clear();
      caloParticlePU_xtalIeta.clear();   
      caloParticlePU_xtalIphi.clear();   
      caloParticlePU_xtalIz.clear();   
      caloParticlePU_xtalIplane.clear(); 
      for(const auto& iCalo : *(puCaloParticle.product()))
      {
          const auto& simClusters = iCalo.simClusters();
          auto hits_and_fractions = simClusters[0]->hits_and_fractions();
          caloParticlePU_nHitsWithES = hits_and_fractions.size(); 
          caloParticlePU_nHits = 0; 
          std::vector<std::pair<DetId, float>> hitsAndEnergies;   
          caloParticlePU_totEnergyWithES = 0.;   
          caloParticlePU_totEnergy = 0.; 
          for(unsigned int i = 0; i < hits_and_fractions.size(); i++){
              float energyWithES=0.;    
              float energy = 0.;
  
              DetId id(hits_and_fractions[i].first);
              if(id.subdetId()!=EcalPreshower){
                 energy = hits_and_fractions[i].second;
                 caloParticlePU_nHits += 1;
              }  
              energyWithES = hits_and_fractions[i].second;  

              if(subtractSignalCalo_){
                 for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                     for(const std::pair<DetId, float>& hit_CaloPart : hitsAndEnergies_CaloPart.at(iCalo)){  
                         if(hit_CaloPart.first.rawId() == DetId(hits_and_fractions[i].first).rawId()){ 
                            energyWithES -= hit_CaloPart.second; 
                            if(id.subdetId()!=EcalPreshower) energy -= hit_CaloPart.second; 
                         }
                     }
                 }         
              }
              caloParticlePU_totEnergyWithES += energyWithES;
              caloParticlePU_totEnergy += energy;

              hitsAndEnergies.push_back(std::make_pair(DetId(hits_and_fractions[i].first),energy));
              if(saveSimhitsPU_){               
                 cell = geometry->getPosition(DetId(hits_and_fractions[i].first)); 
                 caloParticlePU_xtalEnergy.push_back(reduceFloat(energy,nBits_));
                 caloParticlePU_xtalEta.push_back(reduceFloat(cell.eta(),nBits_));
                 caloParticlePU_xtalPhi.push_back(reduceFloat(cell.phi(),nBits_)); 
                     
                 int ieta = -99; 
                 int iphi = -99; 
                 int iz = -99;  
                 int iplane = -99;
                 if(id.subdetId()==EcalBarrel){
                    EBDetId eb_id(id);
                    ieta = eb_id.ieta(); 
                    iphi = eb_id.iphi();  
                    iz = 0;   
                    iplane = 0;
                 }else if(id.subdetId()==EcalEndcap){
                    EEDetId ee_id(id);
                    ieta = ee_id.ix(); 
                    iphi = ee_id.iy();
                    if(ee_id.zside()<0) iz=-1;
                    if(ee_id.zside()>0) iz=1;           
                    iplane = 0;
                 }else if(id.subdetId()==EcalPreshower){
                    ESDetId es_id(id);
                    ieta = es_id.six(); 
                    iphi = es_id.siy();
                    if(es_id.zside()<0) iz=-1;
                    if(es_id.zside()>0) iz=1;           
                    iplane = es_id.plane();  
                 } 
                 caloParticlePU_xtalIeta.push_back(ieta); 
                 caloParticlePU_xtalIphi.push_back(iphi);   
                 caloParticlePU_xtalIz.push_back(iz);
                 caloParticlePU_xtalIplane.push_back(iplane); 
              }  
          } 
          hitsAndEnergies_CaloPartPU.push_back(hitsAndEnergies);
      }   
   }

   if(saveCaloParticlesOOTPU_ && isMC_)
   {
      caloParticleOOTPU_xtalEnergy.clear();
      caloParticleOOTPU_xtalEta.clear();
      caloParticleOOTPU_xtalPhi.clear(); 
      caloParticleOOTPU_xtalIeta.clear();   
      caloParticleOOTPU_xtalIphi.clear();   
      caloParticleOOTPU_xtalIz.clear();   
      caloParticleOOTPU_xtalIplane.clear(); 
      for(const auto& iCalo : *(ootpuCaloParticle.product()))
      {
          const auto& simClusters = iCalo.simClusters();
          auto hits_and_fractions = simClusters[0]->hits_and_fractions();
          caloParticleOOTPU_nHitsWithES = hits_and_fractions.size(); 
          caloParticleOOTPU_nHits = 0; 
          std::vector<std::pair<DetId, float>> hitsAndEnergies;   
          caloParticleOOTPU_totEnergyWithES = 0.;   
          caloParticleOOTPU_totEnergy = 0.; 
          for(unsigned int i = 0; i < hits_and_fractions.size(); i++){
              float energyWithES=0.;    
              float energy = 0.;
  
              DetId id(hits_and_fractions[i].first);
              if(id.subdetId()!=EcalPreshower){
                 energy = hits_and_fractions[i].second;
                 caloParticlePU_nHits += 1;
              }  
              energyWithES = hits_and_fractions[i].second;  

              if(subtractSignalCalo_){
                 for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                     for(const std::pair<DetId, float>& hit_CaloPart : hitsAndEnergies_CaloPart.at(iCalo)){  
                         if(hit_CaloPart.first.rawId() == DetId(hits_and_fractions[i].first).rawId()){ 
                            energyWithES -= hit_CaloPart.second; 
                            if(id.subdetId()!=EcalPreshower) energy -= hit_CaloPart.second; 
                         }
                     }
                 }         
              }
              caloParticleOOTPU_totEnergyWithES += energyWithES;
              caloParticleOOTPU_totEnergy += energy;

              hitsAndEnergies.push_back(std::make_pair(DetId(hits_and_fractions[i].first),energy));
              if(saveSimhitsPU_){               
                 cell = geometry->getPosition(DetId(hits_and_fractions[i].first)); 
                 caloParticleOOTPU_xtalEnergy.push_back(reduceFloat(energy,nBits_));
                 caloParticleOOTPU_xtalEta.push_back(reduceFloat(cell.eta(),nBits_));
                 caloParticleOOTPU_xtalPhi.push_back(reduceFloat(cell.phi(),nBits_)); 
                     
                 int ieta = -99; 
                 int iphi = -99; 
                 int iz = -99;  
                 int iplane = -99;
                 if(id.subdetId()==EcalBarrel){
                    EBDetId eb_id(id);
                    ieta = eb_id.ieta(); 
                    iphi = eb_id.iphi();  
                    iz = 0;   
                    iplane = 0;
                 }else if(id.subdetId()==EcalEndcap){
                    EEDetId ee_id(id);
                    ieta = ee_id.ix(); 
                    iphi = ee_id.iy();
                    if(ee_id.zside()<0) iz=-1;
                    if(ee_id.zside()>0) iz=1;           
                    iplane = 0;
                 }else if(id.subdetId()==EcalPreshower){
                    ESDetId es_id(id);
                    ieta = es_id.six(); 
                    iphi = es_id.siy();
                    if(es_id.zside()<0) iz=-1;
                    if(es_id.zside()>0) iz=1;           
                    iplane = es_id.plane();  
                 } 
                 caloParticleOOTPU_xtalIeta.push_back(ieta); 
                 caloParticleOOTPU_xtalIphi.push_back(iphi);   
                 caloParticleOOTPU_xtalIz.push_back(iz);
                 caloParticleOOTPU_xtalIplane.push_back(iplane); 
              }  
          } 
          hitsAndEnergies_CaloPartOOTPU.push_back(hitsAndEnergies);
      }   
      
   } 

   int nSuperClustersEB = 0;
   int nSuperClustersEE = 0;
   int nRetunedSuperClustersEB = 0;
   int nRetunedSuperClustersEE = 0;
   int nDeepSuperClustersEB = 0;
   int nDeepSuperClustersEE = 0;
   int nPFClusters = 0;
   if(savePFCluster_){
      nPFClusters = (pfClusters.product())->size();
   }
   if(saveSuperCluster_){
      nSuperClustersEB = (superClusterEB.product())->size();
      nSuperClustersEE = (superClusterEE.product())->size();
   }
   if(saveRetunedSC_){
      nRetunedSuperClustersEB = (retunedSuperClusterEB.product())->size();
      nRetunedSuperClustersEE = (retunedSuperClusterEE.product())->size();
   }
   if(saveDeepSC_){ 
      nDeepSuperClustersEB = (deepSuperClusterEB.product())->size();
      nDeepSuperClustersEE = (deepSuperClusterEE.product())->size();
   }

   setVectors(nGenParticles, nCaloParticles, nPFClusters, nSuperClustersEB, nSuperClustersEE, nRetunedSuperClustersEB, nRetunedSuperClustersEE, nDeepSuperClustersEB, nDeepSuperClustersEE); 

   hitsAndEnergies_PFCluster.clear();
   hitsAndEnergies_SuperClusterEB.clear();
   hitsAndEnergies_SuperClusterEE.clear();
   hitsAndEnergies_RetunedSuperClusterEB.clear();
   hitsAndEnergies_RetunedSuperClusterEE.clear();
   hitsAndEnergies_DeepSuperClusterEB.clear();
   hitsAndEnergies_DeepSuperClusterEE.clear();

   std::vector<DetId> hits;
   std::vector<float> energies; 
   hits_CaloParticle.clear();
   energies_CaloParticle.clear();
 
   if(isMC_ && saveCaloParticles_){
   
    for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
   
       caloParticle_index.push_back(iCalo); 
       caloParticle_nXtals.push_back(hitsAndEnergies_CaloPart.at(iCalo).size()); 
       int genIndex = caloParts.at(iCalo).g4Tracks()[0].genpartIndex()-1;     
       const auto& genParts_tmp = *(genParticles.product());    
       auto genParticle = genParts_tmp[genIndex]; 
       int partonIndex = -1;
       if(genParticle.numberOfMothers()!=0) partonIndex = getGenParton(&genParts_tmp,genIndex); 
       if(partonIndex>=0){
          auto genParton = genParts_tmp[partonIndex]; 
          caloParticle_partonIndex.push_back(partonIndex);
          caloParticle_partonPdgId.push_back(genParton.pdgId());
          caloParticle_partonCharge.push_back(genParton.charge());
          caloParticle_partonEnergy.push_back(reduceFloat(genParton.energy(),nBits_));
          caloParticle_partonPt.push_back(reduceFloat(genParton.pt(),nBits_));
          caloParticle_partonEta.push_back(reduceFloat(genParton.eta(),nBits_));
          caloParticle_partonPhi.push_back(reduceFloat(genParton.phi(),nBits_)); 
       }else{
          caloParticle_partonIndex.push_back(-99);
          caloParticle_partonPdgId.push_back(0);
          caloParticle_partonCharge.push_back(-99);
          caloParticle_partonEnergy.push_back(-999.);
          caloParticle_partonPt.push_back(-999.);
          caloParticle_partonEta.push_back(-999.);
          caloParticle_partonPhi.push_back(-999.);
       } 
       if(genParticle.numberOfMothers()!=0){
          caloParticle_genMotherPdgId.push_back(genParticle.mother()->pdgId());
          caloParticle_genMotherStatus.push_back(genParticle.mother()->status());
          caloParticle_genMotherCharge.push_back(genParticle.mother()->charge());
          caloParticle_genMotherEnergy.push_back(reduceFloat(genParticle.mother()->energy(),nBits_));
          caloParticle_genMotherPt.push_back(reduceFloat(genParticle.mother()->pt(),nBits_));
          caloParticle_genMotherEta.push_back(reduceFloat(genParticle.mother()->eta(),nBits_));
          caloParticle_genMotherPhi.push_back(reduceFloat(genParticle.mother()->phi(),nBits_)); 
       }else{
          caloParticle_genMotherPdgId.push_back(0);
          caloParticle_genMotherStatus.push_back(-999);
          caloParticle_genMotherCharge.push_back(-99);
          caloParticle_genMotherEnergy.push_back(-999.);
          caloParticle_genMotherPt.push_back(-999.);
          caloParticle_genMotherEta.push_back(-999.);
          caloParticle_genMotherPhi.push_back(-999.);  
       }
       caloParticle_pdgId.push_back(genParticle.pdgId());
       caloParticle_status.push_back(genParticle.status());
       caloParticle_charge.push_back(genParticle.charge());
       caloParticle_genEnergy.push_back(reduceFloat(genParticle.energy(),nBits_));
       caloParticle_genPt.push_back(reduceFloat(genParticle.pt(),nBits_));
       caloParticle_genEta.push_back(reduceFloat(genParticle.eta(),nBits_));
       caloParticle_genPhi.push_back(reduceFloat(genParticle.phi(),nBits_));
       
       GlobalPoint caloParticle_position = caloParts_position.at(iCalo);
       caloParticle_simEta.push_back(reduceFloat(caloParticle_position.eta(),nBits_));
       caloParticle_simPhi.push_back(reduceFloat(caloParticle_position.phi(),nBits_));
       caloParticle_simPt.push_back(reduceFloat(sqrt(caloParts.at(iCalo).px()*caloParts.at(iCalo).px() + caloParts.at(iCalo).py()*caloParts.at(iCalo).py() + caloParts.at(iCalo).pz()*caloParts.at(iCalo).pz())/TMath::CosH(caloParticle_position.eta()),nBits_));   
       if(std::abs(caloParticle_position.eta()) < 1.479){  
          EBDetId eb_id(_ebGeom->getClosestCell(caloParticle_position));  
          caloParticle_simIeta.push_back(eb_id.ieta());
          caloParticle_simIphi.push_back(eb_id.iphi());
          caloParticle_simIz.push_back(0); 
       }else{            
          int iz=-99;
          EEDetId ee_id(_eeGeom->getClosestCell(caloParticle_position));   
          caloParticle_simIeta.push_back(ee_id.ix());
          caloParticle_simIphi.push_back(ee_id.iy());
          if(ee_id.zside()<0) iz=-1;
          if(ee_id.zside()>0) iz=1;  
          caloParticle_simIz.push_back(iz); 
       } 

       float calo_simEnergy=0.; 
       float calo_simEnergyGoodStatus=0.;  
       float calo_simEnergyWithES=0.;  
       for(auto const& hit: hitsAndEnergies_CaloPart[iCalo])
       {
           DetId id(hit.first);
           int status = (*chStatus->getMap().find(id.rawId())).getStatusCode();
           if(id.subdetId()==EcalBarrel || id.subdetId()==EcalEndcap) calo_simEnergy += hit.second; 
           if((id.subdetId()==EcalBarrel || id.subdetId()==EcalEndcap) && status<3) calo_simEnergyGoodStatus += hit.second; 
           if(id.subdetId()!=EcalBarrel && id.subdetId()!=EcalEndcap && id.subdetId() != EcalPreshower) continue;
               
           calo_simEnergyWithES += hit.second; 

           cell = geometry->getPosition(id);
           float eta = cell.eta();  
           float phi = cell.phi();  
           int ieta = -99; 
           int iphi = -99;
           int iz = -99;  
           int iplane = -99;
           if(id.subdetId()==EcalBarrel){
              EBDetId eb_id(id);
              ieta = eb_id.ieta(); 
              iphi = eb_id.iphi();  
              iz = 0;   
              iplane = 0;
           }else if(id.subdetId()==EcalEndcap){
              EEDetId ee_id(id);
              ieta = ee_id.ix(); 
              iphi = ee_id.iy();
              if(ee_id.zside()<0) iz=-1;
              if(ee_id.zside()>0) iz=1;           
              iplane = 0;
           }else if(id.subdetId()==EcalPreshower){
              ESDetId es_id(id);
              ieta = es_id.six(); 
              iphi = es_id.siy();
              if(es_id.zside()<0) iz=-1;
              if(es_id.zside()>0) iz=1;           
              iplane = es_id.plane(); 
           } 
           if(saveSimhits_ && saveCaloParticles_){
              simHit_energy[iCalo].push_back(reduceFloat(hit.second,nBits_));
              simHit_eta[iCalo].push_back(reduceFloat(eta,nBits_));
              simHit_phi[iCalo].push_back(reduceFloat(phi,nBits_));
              simHit_ieta[iCalo].push_back(ieta);
              simHit_iphi[iCalo].push_back(iphi);
              simHit_iz[iCalo].push_back(iz); 
              simHit_iplane[iCalo].push_back(iplane); 
              simHit_chStatus[iCalo].push_back(status); 
           }
       } 
       caloParticle_simEnergy.push_back(reduceFloat(calo_simEnergy,nBits_));
       caloParticle_simEnergyGoodStatus.push_back(reduceFloat(calo_simEnergyGoodStatus,nBits_));
       caloParticle_simEnergyWithES.push_back(reduceFloat(calo_simEnergyWithES,nBits_));
    }
   }

   //check shared crystals among caloParticles 
   if(isMC_ && saveCaloParticles_){ 

    for(unsigned int i=0; i<hitsAndEnergies_CaloPart.size(); i++){
       for(unsigned int j=0; j<hitsAndEnergies_CaloPart.size(); j++)
       {
           if(i>j || i==j) continue;
           
           std::vector<std::pair<DetId, std::pair<float,float> > >* shareHits =  getSharedHitsAndEnergies(&hitsAndEnergies_CaloPart.at(i), &hitsAndEnergies_CaloPart.at(j));
           
           if((int)shareHits->size()>0){
              double sharedEnergy1 = 0.; 
              double sharedEnergy2 = 0.; 
              for(unsigned int d=0; d<shareHits->size(); d++)
              {
                  sharedEnergy1 += shareHits->at(d).second.first;  
                  sharedEnergy2 += shareHits->at(d).second.second;  
              }  
              caloParticle_sharedIndex1.push_back((int)i);  
              caloParticle_sharedIndex2.push_back((int)j);  
              caloParticle_nSharedXtals.push_back((int)shareHits->size()); 
              caloParticle_sharedEnergyFrac1.push_back(reduceFloat(sharedEnergy1/caloParticle_simEnergyWithES.at(i),nBits_));   
              caloParticle_sharedEnergyFrac2.push_back(reduceFloat(sharedEnergy2/caloParticle_simEnergyWithES.at(j),nBits_));  
           }   
       } 
    } 
   }

   //save hitsAndEnergies for each PFcluster and SuperCluster
   if(savePFCluster_){
      for(const auto& iPFCluster : *(pfClusters.product())){  
          reco::CaloCluster caloBC(iPFCluster);
          hitsAndEnergies_PFCluster.push_back(*getHitsAndEnergiesBC(&caloBC,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
      }
   }
        
   if(saveSuperCluster_){
     for(const auto& iSuperCluster : *(superClusterEB.product())) 
         hitsAndEnergies_SuperClusterEB.push_back(*getHitsAndEnergiesSC(&iSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
     for(const auto& iSuperCluster : *(superClusterEE.product())) 
         hitsAndEnergies_SuperClusterEE.push_back(*getHitsAndEnergiesSC(&iSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
   }

   if(saveRetunedSC_){
     for(const auto& iRetunedSuperCluster : *(retunedSuperClusterEB.product())) 
         hitsAndEnergies_RetunedSuperClusterEB.push_back(*getHitsAndEnergiesSC(&iRetunedSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
     for(const auto& iRetunedSuperCluster : *(retunedSuperClusterEE.product())) 
         hitsAndEnergies_RetunedSuperClusterEE.push_back(*getHitsAndEnergiesSC(&iRetunedSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
   }

   if(saveDeepSC_){
     for(const auto& iSuperCluster : *(deepSuperClusterEB.product())) 
         hitsAndEnergies_DeepSuperClusterEB.push_back(*getHitsAndEnergiesSC(&iSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
     for(const auto& iSuperCluster : *(deepSuperClusterEE.product())) 
         hitsAndEnergies_DeepSuperClusterEE.push_back(*getHitsAndEnergiesSC(&iSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
   }
   
   //Save PFClusters 
   if(savePFCluster_){
     
      int iPFCl=0;
      //std::cout << "PFClusters size     : " << (pfClusters.product())->size() << std::endl;
      for(const auto& iPFCluster : *(pfClusters.product())){  

          pfCluster_rawEnergy.push_back(reduceFloat(iPFCluster.energy(),nBits_));
          pfCluster_energy.push_back(reduceFloat(iPFCluster.correctedEnergy(),nBits_));
          pfCluster_rawPt.push_back(reduceFloat(iPFCluster.energy()/TMath::CosH(iPFCluster.eta()),nBits_));
          pfCluster_pt.push_back(reduceFloat(iPFCluster.correctedEnergy()/TMath::CosH(iPFCluster.eta()),nBits_));
          pfCluster_eta.push_back(reduceFloat(iPFCluster.eta(),nBits_));
          pfCluster_phi.push_back(reduceFloat(iPFCluster.phi(),nBits_));

          reco::CaloCluster caloBC(iPFCluster);

          math::XYZPoint caloPos = caloBC.position();
          if(iPFCluster.layer() == PFLayer::ECAL_BARREL){  
             EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
             pfCluster_ieta.push_back(eb_id.ieta());
             pfCluster_iphi.push_back(eb_id.iphi());
             pfCluster_iz.push_back(0); 
          }else if(iPFCluster.layer() == PFLayer::ECAL_ENDCAP){ 
             int iz=-99;
             EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
             if(ee_id.zside()<0) iz=-1;
             if(ee_id.zside()>0) iz=1;     
             pfCluster_ieta.push_back(ee_id.ix());
             pfCluster_iphi.push_back(ee_id.iy());
             pfCluster_iz.push_back(iz); 
          } 
           
          if(isMC_ && saveCaloParticles_ && saveCaloParticlesPU_){
             std::vector<double> noise = getNoise(&iPFCluster,&hitsAndEnergies_CaloPart,&hitsAndEnergies_CaloPartPU.at(0), &(*(recHitsEB.product())), &(*(recHitsEE.product())), laserAlpha, laserRatio, ical, icalMC, ped, adcToGeV, gr, true);       
             pfCluster_noise.push_back(reduceFloat(noise[1],nBits_));
             pfCluster_noiseUncalib.push_back(reduceFloat(noise[2],nBits_));
             pfCluster_noiseDB.push_back(reduceFloat(noise[3],nBits_));
             pfCluster_noiseDBUncalib.push_back(reduceFloat(noise[4],nBits_));
             pfCluster_rawEnergyUncalib.push_back(reduceFloat(noise[0],nBits_));
       
             std::vector<double> noiseNoFractions = getNoise(&iPFCluster,&hitsAndEnergies_CaloPart,&hitsAndEnergies_CaloPartPU.at(0), &(*(recHitsEB.product())), &(*(recHitsEE.product())), laserAlpha, laserRatio, ical, icalMC, ped, adcToGeV, gr, false);     
             pfCluster_noiseNoFractions.push_back(reduceFloat(noiseNoFractions[1],nBits_));
             pfCluster_noiseUncalibNoFractions.push_back(reduceFloat(noiseNoFractions[2],nBits_)); 
             pfCluster_noiseDBNoFractions.push_back(reduceFloat(noiseNoFractions[3],nBits_));
             pfCluster_noiseDBUncalibNoFractions.push_back(reduceFloat(noiseNoFractions[4],nBits_)); 
          }

          if(saveShowerShapes_ && iPFCluster.layer() == PFLayer::ECAL_BARREL){
             widths_ = calculateCovariances(&iPFCluster, &(*(recHitsEB.product())), &(*_ebGeom));
             pfCluster_etaWidth.push_back(reduceFloat(widths_.first,nBits_));
             pfCluster_phiWidth.push_back(reduceFloat(widths_.second,nBits_));
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), &(*topology));  
             pfCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             pfCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             pfCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     pfCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             pfCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             pfCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             pfCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             pfCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             pfCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             pfCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             pfCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             pfCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             pfCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             pfCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             pfCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             pfCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             pfCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             pfCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             pfCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             pfCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             pfCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             pfCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             pfCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             pfCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             pfCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             pfCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             pfCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             pfCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             pfCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             pfCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             pfCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             pfCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             pfCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             pfCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             pfCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             pfCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             pfCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             pfCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 
          }else if(saveShowerShapes_ && iPFCluster.layer() == PFLayer::ECAL_ENDCAP){ 
             widths_ = calculateCovariances(&iPFCluster, &(*(recHitsEE.product())), &(*_eeGeom));
             pfCluster_etaWidth.push_back(reduceFloat(widths_.first,nBits_));
             pfCluster_phiWidth.push_back(reduceFloat(widths_.second,nBits_)); 
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), &(*topology));  
             pfCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             pfCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             pfCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     pfCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             pfCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             pfCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             pfCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             pfCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             pfCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             pfCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             pfCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             pfCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             pfCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             pfCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             pfCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             pfCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             pfCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             pfCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             pfCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             pfCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             pfCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             pfCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             pfCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             pfCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             pfCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             pfCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             pfCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             pfCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             pfCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             pfCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             pfCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             pfCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             pfCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             pfCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             pfCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             pfCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             pfCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             pfCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 
          }  
          
          if(savePFClusterhits_ && savePFCluster_){
             //for save PFClusterHit    
             const std::vector<std::pair<DetId,float> > &hitsAndFractions = iPFCluster.hitsAndFractions();  
             for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++){      
                 cell = geometry->getPosition(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first);
                 for(unsigned int hits=0; hits<hitsAndFractions.size(); hits++){
                     if(hitsAndFractions.at(hits).first.rawId() == hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId())
                         pfClusterHit_fraction[iPFCl].push_back(hitsAndFractions.at(hits).second);
                 } 
                 float agv = 1.;
                 float ic = *ical->getMap().find(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId()); 
                 float icMC = 1.;
                 if(isMC_) icMC = *icalMC->getMap().find(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId());     
                 double alpha = *laserAlpha->getMap().find(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId());
                 double apdpn = (*laserRatio->getLaserMap().find(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId())).p2;
                 float laserCorr = pow(apdpn,-alpha);
                 int status = (*chStatus->getMap().find(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId())).getStatusCode();
                 pfClusterHit_eta[iPFCl].push_back(reduceFloat(cell.eta(),nBits_));
                 pfClusterHit_phi[iPFCl].push_back(reduceFloat(cell.phi(),nBits_));
                 if(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.subdetId()==EcalBarrel){ 
                    EBDetId eb_id(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first); 
                    pfClusterHit_rechitEnergy[iPFCl].push_back(reduceFloat((*(recHitsEB.product())->find(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first)).energy(),nBits_)); 
                    pfClusterHit_ieta[iPFCl].push_back(eb_id.ieta());
                    pfClusterHit_iphi[iPFCl].push_back(eb_id.iphi());
                    pfClusterHit_iz[iPFCl].push_back(0); 
                    agv = adcToGeV->getEBValue();
                 }else if(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.subdetId()==EcalEndcap){  
                    int iz=-99;
                    EEDetId ee_id(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first);  
                    pfClusterHit_rechitEnergy[iPFCl].push_back(reduceFloat((*(recHitsEE.product())->find(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first)).energy(),nBits_)); 
                    pfClusterHit_ieta[iPFCl].push_back(ee_id.ix());
                    pfClusterHit_iphi[iPFCl].push_back(ee_id.iy());
                    if(ee_id.zside()<0) iz=-1;
                    if(ee_id.zside()>0) iz=1;   
                    pfClusterHit_iz[iPFCl].push_back(iz); 
                    agv = adcToGeV->getEEValue();
                 }    
                 pfClusterHit_adcToGeV[iPFCl].push_back(reduceFloat(agv,nBits_)); 
                 pfClusterHit_laserCorr[iPFCl].push_back(reduceFloat(laserCorr,nBits_)); 
                 pfClusterHit_ic[iPFCl].push_back(reduceFloat(ic,nBits_)); 
                 pfClusterHit_icMC[iPFCl].push_back(reduceFloat(icMC,nBits_)); 
                 pfClusterHit_chStatus[iPFCl].push_back(status); 
             }
          }
          
          //compute scores     
          if(saveGenParticles_ && isMC_){
             dR_genScore.clear();
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             pfCluster_dR_genScore[iPFCl] = dR_genScore;        
             pfCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_dR_genScore, 999., false, 0., iPFCl));
          } 
        
          if(saveCaloParticlesPU_ && isMC_){ 
             simPU_nSharedXtals.clear();
             simEnergy_sharedXtalsPU.clear(); 
             recoEnergy_sharedXtalsPU.clear();
             simEnergy_noHitsFraction_sharedXtalsPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iPFCluster,&hitsAndEnergies_CaloPartPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[8]);   
             } 
             pfCluster_simPU_nSharedXtals[iPFCl] = simPU_nSharedXtals[0];   
             pfCluster_simEnergy_sharedXtalsPU[iPFCl] = simEnergy_sharedXtalsPU[0];   
             pfCluster_recoEnergy_sharedXtalsPU[iPFCl] = recoEnergy_sharedXtalsPU[0];     
             pfCluster_simEnergy_noHitsFraction_sharedXtalsPU[iPFCl] = simEnergy_noHitsFraction_sharedXtalsPU[0];   
             pfCluster_recoEnergy_noHitsFraction_sharedXtalsPU[iPFCl] = recoEnergy_noHitsFraction_sharedXtalsPU[0];        
          }

          if(saveCaloParticlesOOTPU_ && isMC_){ 
             simOOTPU_nSharedXtals.clear();
             simEnergy_sharedXtalsOOTPU.clear(); 
             recoEnergy_sharedXtalsOOTPU.clear();
             simEnergy_noHitsFraction_sharedXtalsOOTPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsOOTPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartOOTPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iPFCluster,&hitsAndEnergies_CaloPartOOTPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simOOTPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsOOTPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsOOTPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[8]);   
             } 
             pfCluster_simOOTPU_nSharedXtals[iPFCl] = simOOTPU_nSharedXtals[0];   
             pfCluster_simEnergy_sharedXtalsOOTPU[iPFCl] = simEnergy_sharedXtalsOOTPU[0];   
             pfCluster_recoEnergy_sharedXtalsOOTPU[iPFCl] = recoEnergy_sharedXtalsOOTPU[0];     
             pfCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU[iPFCl] = simEnergy_noHitsFraction_sharedXtalsOOTPU[0];   
             pfCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU[iPFCl] = recoEnergy_noHitsFraction_sharedXtalsOOTPU[0];        
          }

          pfCluster_nXtals.push_back((double)(iPFCluster.hitsAndFractions()).size());   
    
          if(saveCaloParticles_ && isMC_){ 
             dR_simScore.clear();
             sim_nSharedXtals.clear();
             sim_fraction_noHitsFraction.clear();
             sim_fraction.clear();
             recoToSim_fraction.clear();
             recoToSim_fraction_sharedXtals.clear();  
             simEnergy_sharedXtals.clear(); 
             recoEnergy_sharedXtals.clear(); 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 GlobalPoint caloParticle_position = caloParts_position.at(iCalo);
                 std::vector<double> scores = getScores(&iPFCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             pfCluster_dR_simScore[iPFCl] = dR_simScore;  
             pfCluster_sim_nSharedXtals[iPFCl] = sim_nSharedXtals;   
             pfCluster_sim_fraction_noHitsFraction[iPFCl] = sim_fraction_noHitsFraction;    
             pfCluster_sim_fraction[iPFCl] = sim_fraction;   
             pfCluster_recoToSim_fraction[iPFCl] = recoToSim_fraction;   
             pfCluster_recoToSim_fraction_sharedXtals[iPFCl] = recoToSim_fraction_sharedXtals;        
             pfCluster_simEnergy_sharedXtals[iPFCl] = simEnergy_sharedXtals;   
             pfCluster_recoEnergy_sharedXtals[iPFCl] = recoEnergy_sharedXtals;        

             pfCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_dR_simScore, 999., false, 0., iPFCl));
             pfCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_nSharedXtals, -999., true, 0., iPFCl));
             pfCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_noHitsFraction, -999., true, 0., iPFCl));
             pfCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction, -999., true, 0., iPFCl));
             pfCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_recoToSim_fraction, 999., false, 1., iPFCl)); 
             pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&pfCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iPFCl)); 
             pfCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&pfCluster_simEnergy_sharedXtals, -999., true, 0., iPFCl)); 
             pfCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&pfCluster_recoEnergy_sharedXtals, -999., true, 0., iPFCl));
             
          }    
          iPFCl++;        
      }

      //save inverse of matchings
      if(saveGenParticles_ && isMC_){ 
         fillParticleMatchedIndex(&genParticle_pfCluster_dR_genScore_MatchedIndex,&pfCluster_dR_genScore_MatchedIndex);
      } 
      if(saveCaloParticles_ && isMC_){ 
         fillParticleMatchedIndex(&caloParticle_pfCluster_dR_simScore_MatchedIndex,&pfCluster_dR_simScore_MatchedIndex);
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_nSharedXtals_MatchedIndex,&pfCluster_sim_nSharedXtals_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_noHitsFraction_MatchedIndex,&pfCluster_sim_fraction_noHitsFraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_MatchedIndex,&pfCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_recoToSim_fraction_MatchedIndex,&pfCluster_recoToSim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex,&pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_simEnergy_sharedXtals_MatchedIndex,&pfCluster_simEnergy_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_recoEnergy_sharedXtals_MatchedIndex,&pfCluster_recoEnergy_sharedXtals_MatchedIndex);  
      }
   } 
   
   if(saveGsfElectrons_)
   {
      //save gsfElectron infos
      int iEle=0;
      for(const auto& iElectron : *(gsfElectron.product())){ 
       
          reco::SuperClusterRef scRef = iElectron.superCluster();
          double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
          double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y()); 

          reco::SuperCluster matchedEcalSC; 
          bool isMatchedInEB = false;  
          bool isMatchedInEE = false;  
          for(const auto& iSuperCluster : *(superClusterEB.product())){
              if(scRef->seed()->seed().rawId()==iSuperCluster.seed()->seed().rawId()){ 
                 matchedEcalSC = iSuperCluster;   
                 isMatchedInEB = true;
                 break;
              }     
          }    
          if(!isMatchedInEB){
             for(const auto& iSuperCluster : *(superClusterEE.product())){
                 if(scRef->seed()->seed().rawId()==iSuperCluster.seed()->seed().rawId()){ 
                    matchedEcalSC = iSuperCluster;  
                    isMatchedInEE = true; 
                    break;
                 }     
             }
          }    
          double ecalR  = TMath::Sqrt(matchedEcalSC.x()*matchedEcalSC.x() + matchedEcalSC.y()*matchedEcalSC.y() + matchedEcalSC.z()*matchedEcalSC.z());
          double ecalRt = TMath::Sqrt(matchedEcalSC.x()*matchedEcalSC.x() + matchedEcalSC.y()*matchedEcalSC.y());  

          double swissCross = -999.;
          double eMax = -999.;
          double e2x2 = -999.;
          double e3x3 = -999.;
          double e5x5 = -999.; 
          double full5x5_e2x2 = -999.; 
          double full5x5_e3x3 = -999.; 
          double full5x5_e5x5 = -999.;  
          double full5x5_eMax = -999.; 
 
          reco::GsfElectron::ShowerShape eleSS = iElectron.showerShape(); 
          reco::GsfElectron::ShowerShape full5x5_eleSS = iElectron.full5x5_showerShape();
          const std::vector<std::pair<DetId,float> > &hits= iElectron.superCluster()->hitsAndFractions();
          if(iElectron.isEB())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEB.product())));  
             swissCross = EcalTools::swissCross(id.first,*(recHitsEB.product()),0.);
             e2x2 = EcalClusterTools::e2x2( *scRef, &(*(recHitsEB.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *scRef, &(*(recHitsEB.product())), topology);
             e5x5 = EcalClusterTools::e5x5( *scRef, &(*(recHitsEB.product())), topology);
             eMax = EcalClusterTools::eMax( *scRef, &(*(recHitsEB.product())));      
             full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *scRef, &(*(recHitsEB.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *scRef, &(*(recHitsEB.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *scRef, &(*(recHitsEB.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *scRef, &(*(recHitsEB.product())));     
          }
          if(iElectron.isEE())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEE.product())));  
             swissCross = EcalTools::swissCross(id.first,*(recHitsEE.product()),0.);
             e2x2 = EcalClusterTools::e2x2( *scRef, &(*(recHitsEE.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *scRef, &(*(recHitsEE.product())), topology);
             e5x5 = EcalClusterTools::e5x5( *scRef, &(*(recHitsEE.product())), topology);
             eMax = EcalClusterTools::eMax( *scRef, &(*(recHitsEE.product())));      
             full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *scRef, &(*(recHitsEE.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *scRef, &(*(recHitsEE.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *scRef, &(*(recHitsEE.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *scRef, &(*(recHitsEE.product())));     
          }

          gsfElectron_index.push_back(iEle);
          gsfElectron_seedRawId.push_back(iElectron.superCluster()->seed()->seed().rawId());
          gsfElectron_isEB.push_back(iElectron.isEB()); 
          gsfElectron_isEE.push_back(iElectron.isEE());  
          gsfElectron_isEBEEGap.push_back(iElectron.isEBEEGap());  
          gsfElectron_isEBEtaGap.push_back(iElectron.isEBEtaGap());  
          gsfElectron_isEBPhiGap.push_back(iElectron.isEBPhiGap());  
          gsfElectron_isEEDeeGap.push_back(iElectron.isEEDeeGap());  
          gsfElectron_isEERingGap.push_back(iElectron.isEERingGap());   
          gsfElectron_isEcalDriven.push_back(iElectron.ecalDrivenSeed()); 
          gsfElectron_isTrackerDriven.push_back(iElectron.trackerDrivenSeed()); 
          gsfElectron_classification.push_back(iElectron.classification()); 
          gsfElectron_scNPFClusters.push_back(scRef->clusters().size());
          if(isMatchedInEB || isMatchedInEE) gsfElectron_ecalSCNPFClusters.push_back(matchedEcalSC.clusters().size());
          else gsfElectron_ecalSCNPFClusters.push_back(-1);  
          gsfElectron_p.push_back(reduceFloat(iElectron.trackMomentumAtVtx().R(),nBits_));
          gsfElectron_pt.push_back(reduceFloat(TMath::Sqrt(iElectron.trackMomentumAtVtx().Perp2()),nBits_));
          gsfElectron_et.push_back(reduceFloat(iElectron.et(),nBits_));
          gsfElectron_energy.push_back(reduceFloat(iElectron.energy(),nBits_));
          gsfElectron_energyErr.push_back(reduceFloat(iElectron.p4Error(reco::GsfElectron::P4_COMBINATION),nBits_));
          gsfElectron_ecalEnergy.push_back(reduceFloat(iElectron.ecalEnergy(),nBits_));
          gsfElectron_ecalEnergyErr.push_back(reduceFloat(iElectron.ecalEnergyError(),nBits_));
          gsfElectron_eta.push_back(reduceFloat(iElectron.eta(),nBits_));
          gsfElectron_phi.push_back(reduceFloat(iElectron.phi(),nBits_));
          gsfElectron_trkEtaMode.push_back(reduceFloat(iElectron.gsfTrack()->etaMode(),nBits_));
          gsfElectron_trkPhiMode.push_back(reduceFloat(iElectron.gsfTrack()->phiMode(),nBits_));
          gsfElectron_trkPMode.push_back(reduceFloat(iElectron.gsfTrack()->pMode(),nBits_));
          gsfElectron_trkPModeErr.push_back(reduceFloat(std::abs(iElectron.gsfTrack()->qoverpModeError())*iElectron.gsfTrack()->pMode()*iElectron.gsfTrack()->pMode(),nBits_));
          gsfElectron_trkPInn.push_back(reduceFloat(iElectron.gsfTrack()->p(),nBits_));
          gsfElectron_trkPtInn.push_back(reduceFloat(iElectron.gsfTrack()->pt(),nBits_));
          gsfElectron_trkPVtx.push_back(reduceFloat(std::sqrt(iElectron.trackMomentumAtVtx().Mag2()),nBits_));
          gsfElectron_trkPOut.push_back(reduceFloat(std::sqrt(iElectron.trackMomentumOut().Mag2()),nBits_));
          gsfElectron_trkChi2.push_back(reduceFloat(iElectron.gsfTrack()->chi2(),nBits_));
          gsfElectron_trkNDof.push_back(reduceFloat(iElectron.gsfTrack()->ndof(),nBits_));
          gsfElectron_trackFbrem.push_back(reduceFloat(iElectron.trackFbrem(),nBits_));
          gsfElectron_superClusterFbrem.push_back(reduceFloat(iElectron.superClusterFbrem(),nBits_));
          gsfElectron_hademTow.push_back(reduceFloat(iElectron.hcalOverEcalBc(),nBits_));
          gsfElectron_hademCone.push_back(reduceFloat(iElectron.hcalOverEcal(),nBits_));
          gsfElectron_ecalDrivenSeed.push_back(reduceFloat(iElectron.ecalDrivenSeed(),nBits_));
          gsfElectron_nrSatCrys.push_back(reduceFloat(iElectron.nSaturatedXtals(),nBits_));
          gsfElectron_scEnergy.push_back(reduceFloat(scRef->energy(),nBits_)); 
          gsfElectron_scRawEnergy.push_back(reduceFloat(scRef->rawEnergy(),nBits_)); 
          gsfElectron_scRawESEnergy.push_back(reduceFloat(scRef->preshowerEnergy(),nBits_)); 
          gsfElectron_scEt.push_back(reduceFloat(scRef->energy()*(Rt/R),nBits_));
          gsfElectron_scPhiWidth.push_back(reduceFloat(scRef->phiWidth(),nBits_));
          gsfElectron_scEtaWidth.push_back(reduceFloat(scRef->etaWidth(),nBits_)); 
          gsfElectron_scEoP.push_back(reduceFloat(scRef->energy()/iElectron.trackMomentumAtVtx().R(),nBits_)); 
          gsfElectron_scEta.push_back(reduceFloat(scRef->eta(),nBits_));
          gsfElectron_scPhi.push_back(reduceFloat(scRef->phi(),nBits_));  
          if(isMatchedInEB || isMatchedInEE){
             gsfElectron_ecalSCEta.push_back(reduceFloat(matchedEcalSC.eta(),nBits_)); 
             gsfElectron_ecalSCPhi.push_back(reduceFloat(matchedEcalSC.phi(),nBits_)); 
             gsfElectron_ecalSCEnergy.push_back(reduceFloat(matchedEcalSC.energy(),nBits_)); 
             gsfElectron_ecalSCRawEnergy.push_back(reduceFloat(matchedEcalSC.rawEnergy(),nBits_)); 
             gsfElectron_ecalSCRawESEnergy.push_back(reduceFloat(matchedEcalSC.preshowerEnergy(),nBits_)); 
             gsfElectron_ecalSCEt.push_back(reduceFloat(matchedEcalSC.energy()*(ecalRt/ecalR),nBits_));
             gsfElectron_ecalSCPhiWidth.push_back(reduceFloat(matchedEcalSC.phiWidth(),nBits_));
             gsfElectron_ecalSCEtaWidth.push_back(reduceFloat(matchedEcalSC.etaWidth(),nBits_)); 
             gsfElectron_ecalSCEoP.push_back(reduceFloat(matchedEcalSC.energy()/iElectron.trackMomentumAtVtx().R(),nBits_));  
          }else{
             gsfElectron_ecalSCEta.push_back(-999.); 
             gsfElectron_ecalSCPhi.push_back(-999.); 
             gsfElectron_ecalSCEnergy.push_back(-999.); 
             gsfElectron_ecalSCRawEnergy.push_back(-999.); 
             gsfElectron_ecalSCRawESEnergy.push_back(-999.); 
             gsfElectron_ecalSCEt.push_back(-999.);
             gsfElectron_ecalSCPhiWidth.push_back(-999.);
             gsfElectron_ecalSCEtaWidth.push_back(-999.); 
             gsfElectron_ecalSCEoP.push_back(-999.);  
          } 
          gsfElectron_scSwissCross.push_back(reduceFloat(swissCross,nBits_));  
          gsfElectron_scEMax.push_back(reduceFloat(eMax,nBits_));  
          gsfElectron_scE2x2.push_back(reduceFloat(e2x2,nBits_));
          gsfElectron_scE3x3.push_back(reduceFloat(e3x3,nBits_));
          gsfElectron_scE5x5.push_back(reduceFloat(e5x5,nBits_));
          gsfElectron_scR9.push_back(reduceFloat(e3x3/scRef->rawEnergy(),nBits_));
          gsfElectron_scSigmaIEtaIEta.push_back(reduceFloat(eleSS.sigmaIetaIeta,nBits_));
          gsfElectron_scSigmaIEtaIPhi.push_back(reduceFloat(eleSS.sigmaIetaIphi,nBits_));
          gsfElectron_scSigmaIPhiIPhi.push_back(reduceFloat(eleSS.sigmaIphiIphi,nBits_));   
          gsfElectron_full5x5_scEMax.push_back(reduceFloat(full5x5_eMax,nBits_));  
          gsfElectron_full5x5_scE2x2.push_back(reduceFloat(full5x5_e2x2,nBits_));
          gsfElectron_full5x5_scE3x3.push_back(reduceFloat(full5x5_e3x3,nBits_));
          gsfElectron_full5x5_scE5x5.push_back(reduceFloat(full5x5_e5x5,nBits_));
          gsfElectron_full5x5_scR9.push_back(reduceFloat(full5x5_e3x3/scRef->rawEnergy(),nBits_));
          gsfElectron_full5x5_scSigmaIEtaIEta.push_back(reduceFloat(full5x5_eleSS.sigmaIetaIeta,nBits_));
          gsfElectron_full5x5_scSigmaIEtaIPhi.push_back(reduceFloat(full5x5_eleSS.sigmaIetaIphi,nBits_));
          gsfElectron_full5x5_scSigmaIPhiIPhi.push_back(reduceFloat(full5x5_eleSS.sigmaIphiIphi,nBits_)); 
          gsfElectron_HoE.push_back(reduceFloat(iElectron.hcalOverEcal(),nBits_));   
          gsfElectron_trkIso03.push_back(reduceFloat(iElectron.dr03TkSumPt(),nBits_));    
          gsfElectron_ecalIso03.push_back(reduceFloat(iElectron.dr03EcalRecHitSumEt(),nBits_));   
          gsfElectron_hcalIso03.push_back(reduceFloat(iElectron.dr03HcalTowerSumEt(),nBits_));  
          gsfElectron_trkIso04.push_back(reduceFloat(iElectron.dr04TkSumPt(),nBits_));    
          gsfElectron_ecalIso04.push_back(reduceFloat(iElectron.dr04EcalRecHitSumEt(),nBits_));   
          gsfElectron_hcalIso04.push_back(reduceFloat(iElectron.dr04HcalTowerSumEt(),nBits_));  
          gsfElectron_pfPhotonIso.push_back(reduceFloat(iElectron.pfIsolationVariables().sumPhotonEt,nBits_));      
          gsfElectron_pfChargedHadronIso.push_back(reduceFloat(iElectron.pfIsolationVariables().sumChargedHadronPt,nBits_));      
          gsfElectron_pfNeutralHadronIso.push_back(reduceFloat(iElectron.pfIsolationVariables().sumNeutralHadronEt,nBits_));   
          gsfElectron_mva_Isolated.push_back(reduceFloat(iElectron.mva_Isolated(),nBits_));   
          gsfElectron_mva_e_pi.push_back(reduceFloat(iElectron.mva_e_pi(),nBits_));   
          gsfElectron_dnn_signal_Isolated.push_back(reduceFloat(iElectron.dnn_signal_Isolated(),nBits_));   
          gsfElectron_dnn_signal_nonIsolated.push_back(reduceFloat(iElectron.dnn_signal_nonIsolated(),nBits_));   
          gsfElectron_dnn_bkg_nonIsolated.push_back(reduceFloat(iElectron.dnn_bkg_nonIsolated(),nBits_));   
          gsfElectron_dnn_bkg_Tau.push_back(reduceFloat(iElectron.dnn_bkg_Tau(),nBits_));   
          gsfElectron_dnn_bkg_Photon.push_back(reduceFloat(iElectron.dnn_bkg_Photon(),nBits_));   
    
          iEle++; 
      }
   }
 
   if(saveGedPhotons_)
   {
      //save gedPhoton infos
      int iPho=0;
      for(const auto& iPhoton : *(gedPhoton.product())){ 
       
          reco::SuperClusterRef scRef = iPhoton.superCluster();
          reco::PhotonCoreRef phoCoreRef = iPhoton.photonCore();
          double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
          double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());

          reco::SuperCluster matchedEcalSC; 
          bool isMatchedInEB = false;  
          bool isMatchedInEE = false;  
          for(const auto& iSuperCluster : *(superClusterEB.product())){
              if(scRef->seed()->seed().rawId()==iSuperCluster.seed()->seed().rawId()){ 
                 matchedEcalSC = iSuperCluster;   
                 isMatchedInEB = true;
                 break;
              }     
          }    
          if(!isMatchedInEB){
             for(const auto& iSuperCluster : *(superClusterEE.product())){
                 if(scRef->seed()->seed().rawId()==iSuperCluster.seed()->seed().rawId()){ 
                    matchedEcalSC = iSuperCluster;  
                    isMatchedInEE = true; 
                    break;
                 }     
             }
          }  
          double ecalR  = TMath::Sqrt(matchedEcalSC.x()*matchedEcalSC.x() + matchedEcalSC.y()*matchedEcalSC.y() + matchedEcalSC.z()*matchedEcalSC.z());
          double ecalRt = TMath::Sqrt(matchedEcalSC.x()*matchedEcalSC.x() + matchedEcalSC.y()*matchedEcalSC.y());  

          double swissCross = -999.;
          double eMax = -999.;
          double e2x2 = -999.;
          double e3x3 = -999.;
          double e5x5 = -999.; 
          double full5x5_e2x2 = -999.; 
          double full5x5_e3x3 = -999.; 
          double full5x5_e5x5 = -999.;  
          double full5x5_eMax = -999.; 

          reco::Photon::ShowerShape phoSS = iPhoton.showerShapeVariables();   
          reco::Photon::ShowerShape full5x5_phoSS = iPhoton.full5x5_showerShapeVariables(); 
          const std::vector<std::pair<DetId,float> > &hits= iPhoton.superCluster()->hitsAndFractions();
          if(iPhoton.isEB())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEB.product())));  
             swissCross = EcalTools::swissCross(id.first,*(recHitsEB.product()),0.);
             e2x2 = EcalClusterTools::e2x2( *scRef, &(*(recHitsEB.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *scRef, &(*(recHitsEB.product())), topology);
             e5x5 = EcalClusterTools::e5x5( *scRef, &(*(recHitsEB.product())), topology);
             eMax = EcalClusterTools::eMax( *scRef, &(*(recHitsEB.product())));      
             full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *scRef, &(*(recHitsEB.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *scRef, &(*(recHitsEB.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *scRef, &(*(recHitsEB.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *scRef, &(*(recHitsEB.product())));      
          }
          if(iPhoton.isEE())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEE.product())));  
             swissCross = EcalTools::swissCross(id.first,*(recHitsEE.product()),0.);
             e2x2 = EcalClusterTools::e2x2( *scRef, &(*(recHitsEE.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *scRef, &(*(recHitsEE.product())), topology);
             e5x5 = EcalClusterTools::e5x5( *scRef, &(*(recHitsEE.product())), topology);
             eMax = EcalClusterTools::eMax( *scRef, &(*(recHitsEE.product())));      
             full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *scRef, &(*(recHitsEE.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *scRef, &(*(recHitsEE.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *scRef, &(*(recHitsEE.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *scRef, &(*(recHitsEE.product())));         
          }

          gedPhoton_index.push_back(iPho);  
          gedPhoton_seedRawId.push_back(iPhoton.superCluster()->seed()->seed().rawId());
          gedPhoton_isEB.push_back(iPhoton.isEB());
          gedPhoton_isEE.push_back(iPhoton.isEE());  
          gedPhoton_isEBEEGap.push_back(iPhoton.isEBEEGap());  
          gedPhoton_isEBEtaGap.push_back(iPhoton.isEBEtaGap());  
          gedPhoton_isEBPhiGap.push_back(iPhoton.isEBPhiGap());  
          gedPhoton_isEEDeeGap.push_back(iPhoton.isEEDeeGap());  
          gedPhoton_isEERingGap.push_back(iPhoton.isEERingGap());   
          gedPhoton_scNPFClusters.push_back(scRef->clusters().size());
          if(isMatchedInEB || isMatchedInEE) gedPhoton_ecalSCNPFClusters.push_back(matchedEcalSC.clusters().size());  
          else gedPhoton_ecalSCNPFClusters.push_back(-1);  
          gedPhoton_hasConversionTracks.push_back(iPhoton.hasConversionTracks());
          gedPhoton_nConversions.push_back(phoCoreRef->conversions().size());
          gedPhoton_nConversionsOneLeg.push_back(phoCoreRef->conversionsOneLeg().size());
          gedPhoton_et.push_back(reduceFloat(iPhoton.et(),nBits_));
          gedPhoton_energy.push_back(reduceFloat(iPhoton.energy(),nBits_));
          gedPhoton_energyErr.push_back(reduceFloat(iPhoton.getCorrectedEnergyError(reco::Photon::regression2),nBits_));
          gedPhoton_ecalEnergy.push_back(reduceFloat(iPhoton.energyCorrections().phoEcalEnergy,nBits_));
          gedPhoton_ecalEnergyErr.push_back(reduceFloat(iPhoton.energyCorrections().phoEcalEnergyError,nBits_));
          gedPhoton_eta.push_back(reduceFloat(iPhoton.eta(),nBits_));
          gedPhoton_phi.push_back(reduceFloat(iPhoton.phi(),nBits_));
          gedPhoton_hademTow.push_back(reduceFloat(iPhoton.hadTowOverEm(),nBits_));
          gedPhoton_hademCone.push_back(reduceFloat(iPhoton.hadronicOverEm(),nBits_));
          gedPhoton_nrSatCrys.push_back(reduceFloat(iPhoton.nSaturatedXtals(),nBits_));
          gedPhoton_scEta.push_back(reduceFloat(scRef->eta(),nBits_));
          gedPhoton_scPhi.push_back(reduceFloat(scRef->phi(),nBits_));
          gedPhoton_scEnergy.push_back(reduceFloat(scRef->energy(),nBits_));
          gedPhoton_scRawEnergy.push_back(reduceFloat(scRef->rawEnergy(),nBits_));
          gedPhoton_scRawESEnergy.push_back(reduceFloat(scRef->preshowerEnergy(),nBits_));
          gedPhoton_scEt.push_back(reduceFloat(scRef->energy()*(Rt/R),nBits_));
          gedPhoton_scEtaWidth.push_back(reduceFloat(scRef->etaWidth(),nBits_));  
          gedPhoton_scPhiWidth.push_back(reduceFloat(scRef->phiWidth(),nBits_));
          if(isMatchedInEB || isMatchedInEE){
             gedPhoton_ecalSCEta.push_back(reduceFloat(matchedEcalSC.eta(),nBits_));
             gedPhoton_ecalSCPhi.push_back(reduceFloat(matchedEcalSC.phi(),nBits_)); 
             gedPhoton_ecalSCEnergy.push_back(reduceFloat(matchedEcalSC.energy(),nBits_)); 
             gedPhoton_ecalSCRawEnergy.push_back(reduceFloat(matchedEcalSC.rawEnergy(),nBits_)); 
             gedPhoton_ecalSCRawESEnergy.push_back(reduceFloat(matchedEcalSC.preshowerEnergy(),nBits_)); 
             gedPhoton_ecalSCEt.push_back(reduceFloat(matchedEcalSC.energy()*(ecalRt/ecalR),nBits_));
             gedPhoton_ecalSCPhiWidth.push_back(reduceFloat(matchedEcalSC.phiWidth(),nBits_));
             gedPhoton_ecalSCEtaWidth.push_back(reduceFloat(matchedEcalSC.etaWidth(),nBits_)); 
          }else{
             gedPhoton_ecalSCEta.push_back(-999.);
             gedPhoton_ecalSCPhi.push_back(-999.); 
             gedPhoton_ecalSCEnergy.push_back(-999.); 
             gedPhoton_ecalSCRawEnergy.push_back(-999.); 
             gedPhoton_ecalSCRawESEnergy.push_back(-999.); 
             gedPhoton_ecalSCEt.push_back(-999.);
             gedPhoton_ecalSCPhiWidth.push_back(-999.);
             gedPhoton_ecalSCEtaWidth.push_back(-999.); 
          }
          gedPhoton_scSwissCross.push_back(reduceFloat(swissCross,nBits_));  
          gedPhoton_scEMax.push_back(reduceFloat(eMax,nBits_));  
          gedPhoton_scE2x2.push_back(reduceFloat(e2x2,nBits_));
          gedPhoton_scE3x3.push_back(reduceFloat(e3x3,nBits_));
          gedPhoton_scE5x5.push_back(reduceFloat(e5x5,nBits_));
          gedPhoton_scR9.push_back(reduceFloat(e3x3/scRef->rawEnergy(),nBits_));
          gedPhoton_scSigmaIEtaIEta.push_back(reduceFloat(phoSS.sigmaIetaIeta,nBits_));
          gedPhoton_scSigmaIEtaIPhi.push_back(reduceFloat(phoSS.sigmaIetaIphi,nBits_));
          gedPhoton_scSigmaIPhiIPhi.push_back(reduceFloat(phoSS.sigmaIphiIphi,nBits_));    
          gedPhoton_full5x5_scEMax.push_back(reduceFloat(full5x5_eMax,nBits_));  
          gedPhoton_full5x5_scE2x2.push_back(reduceFloat(full5x5_e2x2,nBits_));
          gedPhoton_full5x5_scE3x3.push_back(reduceFloat(full5x5_e3x3,nBits_));
          gedPhoton_full5x5_scE5x5.push_back(reduceFloat(full5x5_e5x5,nBits_));
          gedPhoton_full5x5_scR9.push_back(reduceFloat(full5x5_e3x3/scRef->rawEnergy(),nBits_));        
          gedPhoton_full5x5_scSigmaIEtaIEta.push_back(reduceFloat(full5x5_phoSS.sigmaIetaIeta,nBits_));
          gedPhoton_full5x5_scSigmaIEtaIPhi.push_back(reduceFloat(full5x5_phoSS.sigmaIetaIphi,nBits_));
          gedPhoton_full5x5_scSigmaIPhiIPhi.push_back(reduceFloat(full5x5_phoSS.sigmaIphiIphi,nBits_));    
          gedPhoton_HoE.push_back(reduceFloat(iPhoton.hcalOverEcal(),nBits_));    
          gedPhoton_trkIso03.push_back(reduceFloat(iPhoton.trkSumPtSolidConeDR03(),nBits_));    
          gedPhoton_ecalIso03.push_back(reduceFloat(iPhoton.ecalRecHitSumEtConeDR03(),nBits_));   
          gedPhoton_hcalIso03.push_back(reduceFloat(iPhoton.hcalTowerSumEtConeDR03(),nBits_)); 
          gedPhoton_trkIso04.push_back(reduceFloat(iPhoton.trkSumPtSolidConeDR04(),nBits_));    
          gedPhoton_ecalIso04.push_back(reduceFloat(iPhoton.ecalRecHitSumEtConeDR04(),nBits_));   
          gedPhoton_hcalIso04.push_back(reduceFloat(iPhoton.hcalTowerSumEtConeDR04(),nBits_));  
          gedPhoton_pfPhotonIso.push_back(reduceFloat(iPhoton.photonIso(),nBits_));      
          gedPhoton_pfChargedHadronIso.push_back(reduceFloat(iPhoton.chargedHadronIso(),nBits_));      
          gedPhoton_pfNeutralHadronIso.push_back(reduceFloat(iPhoton.neutralHadronIso(),nBits_));  
          gedPhoton_nClusterOutsideMustache.push_back(reduceFloat(iPhoton.nClusterOutsideMustache(),nBits_)); 
          gedPhoton_etOutsideMustache.push_back(reduceFloat(iPhoton.etOutsideMustache(),nBits_)); 
          gedPhoton_pfMVA.push_back(reduceFloat(iPhoton.pfMVA(),nBits_));  
          gedPhoton_pfDNN.push_back(reduceFloat(iPhoton.pfDNN(),nBits_));      

          iPho++; 
      }
   }

   //save MET info
   pat::METRef metP;
   if(savePatElectrons_ || savePatPhotons_ || savePatJets_)
   { 
      metP = pat::METRef(patMET,0);
      patMET_sumEt = reduceFloat(metP->sumEt(),nBits_);
      patMET_et   = reduceFloat(metP->corPt(),nBits_); 
   }

   //save patElectron info
   if(savePatElectrons_)
   { 
      //save Mll info
      reco::Candidate::LorentzVector L1;
      reco::Candidate::LorentzVector L2;
      bool found1 = false;
      bool found2 = false;
 
      mll = -999.;
      for(const auto& iElectron : *(patElectron.product())){  
          if(found1 == false){
             L1 = iElectron.p4();
             found1 = true;
          }else if(found2 == false){
             L2 = iElectron.p4();
             found2 = true;
             mll = reduceFloat((L1+L2).mass(),nBits_);
          }
      } 

      int iEle=0;
      for(const auto& iElectron : *(patElectron.product())){ 
 
          reco::SuperClusterRef scRef = iElectron.superCluster();
          double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
          double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());

          reco::SuperCluster matchedEcalSC; 
          bool isMatchedInEB = false;  
          bool isMatchedInEE = false;  
          for(const auto& iSuperCluster : *(superClusterEB.product())){
              if(scRef->seed()->seed().rawId()==iSuperCluster.seed()->seed().rawId()){ 
                 matchedEcalSC = iSuperCluster;   
                 isMatchedInEB = true;
                 break;
              }     
          }    
          if(!isMatchedInEB){
             for(const auto& iSuperCluster : *(superClusterEE.product())){
                 if(scRef->seed()->seed().rawId()==iSuperCluster.seed()->seed().rawId()){ 
                    matchedEcalSC = iSuperCluster;  
                    isMatchedInEE = true; 
                    break;
                 }     
             }
          }   
          double ecalR  = TMath::Sqrt(matchedEcalSC.x()*matchedEcalSC.x() + matchedEcalSC.y()*matchedEcalSC.y() + matchedEcalSC.z()*matchedEcalSC.z());
          double ecalRt = TMath::Sqrt(matchedEcalSC.x()*matchedEcalSC.x() + matchedEcalSC.y()*matchedEcalSC.y());  

          double swissCross = -999.;
          double e3x3 = -999.;
          double e2x2 = -999.;
          double e5x5 = -999.;  
          double eMax = -999.;
          double full5x5_e3x3 = -999.;
          double full5x5_e2x2 = -999.;
          double full5x5_e5x5 = -999.;  
          double full5x5_eMax = -999.; 
         
          reco::GsfElectron::ShowerShape eleSS = iElectron.showerShape();
          reco::GsfElectron::ShowerShape full5x5_eleSS = iElectron.full5x5_showerShape();
          const std::vector<std::pair<DetId,float> > &hits= iElectron.superCluster()->hitsAndFractions();
          if(iElectron.isEB())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEB.product())));   
             swissCross = EcalTools::swissCross(id.first,*(recHitsEB.product()),0.);
             e2x2 = EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             e5x5 = EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));             
             full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));
          }
          if(iElectron.isEE())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEE.product())));   
             swissCross = EcalTools::swissCross(id.first,*(recHitsEE.product()),0.);
             e2x2 = EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             e5x5 = EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
             full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
          }
              
          double cphi = (iElectron.p4().x() * metP->corPx() + iElectron.p4().y() * metP->corPy()) / (metP->corPt()*iElectron.p4().Pt());
          double Mt = sqrt(2 * iElectron.p4().Et() * metP->corPt() * (1-cphi));
          double dphiMET = deltaPhi(metP->corPhi(),iElectron.p4().phi());

          unsigned int nPhotons = 0;
          if(iElectron.hasOverlaps("photons")) nPhotons = iElectron.overlaps("photons").size();  
          std::vector<int> photonIndices;
          photonIndices.resize(nPhotons);
          for(unsigned int iPho=0; iPho<nPhotons; iPho++)
              photonIndices[iPho] = iElectron.overlaps("photons")[iPho].key();

          unsigned int nJets = 0;
          int jetIndex = -1;
          if(iElectron.hasUserCand("jet")){ 
             nJets = 1;
             jetIndex = iElectron.userCand("jet").key();
          }  
          
          patElectron_index.push_back(iEle);
          patElectron_seedRawId.push_back(iElectron.superCluster()->seed()->seed().rawId());
          patElectron_classification.push_back(iElectron.classification());
          patElectron_scNPFClusters.push_back(scRef->clusters().size());
          if(isMatchedInEB || isMatchedInEE) patElectron_ecalSCNPFClusters.push_back(matchedEcalSC.clusters().size());  
          else patElectron_ecalSCNPFClusters.push_back(-1);   
          patElectron_charge.push_back(iElectron.charge()); 
          patElectron_isEB.push_back(iElectron.isEB()); 
          patElectron_isEE.push_back(iElectron.isEE());  
          patElectron_isEBEEGap.push_back(iElectron.isEBEEGap());  
          patElectron_isEBEtaGap.push_back(iElectron.isEBEtaGap());  
          patElectron_isEBPhiGap.push_back(iElectron.isEBPhiGap());  
          patElectron_isEEDeeGap.push_back(iElectron.isEEDeeGap());  
          patElectron_isEERingGap.push_back(iElectron.isEERingGap());   
          patElectron_isEcalDriven.push_back(iElectron.ecalDrivenSeed()); 
          patElectron_isTrackerDriven.push_back(iElectron.trackerDrivenSeed()); 
          patElectron_passConversionVeto.push_back(iElectron.passConversionVeto()); 
          patElectron_nOverlapPhotons.push_back(nPhotons);  
          patElectron_overlapPhotonIndices.push_back(photonIndices);   
          patElectron_hasOverlapJet.push_back(nJets);  
          patElectron_overlapJetIndex.push_back(jetIndex);  
          patElectron_eta.push_back(reduceFloat(iElectron.p4().eta(),nBits_));
          patElectron_phi.push_back(reduceFloat(iElectron.p4().phi(),nBits_));
          patElectron_p.push_back(reduceFloat(iElectron.trackMomentumAtVtx().R(),nBits_));
          patElectron_pt.push_back(reduceFloat(TMath::Sqrt(iElectron.trackMomentumAtVtx().Perp2()),nBits_));
          patElectron_pIn.push_back(reduceFloat(iElectron.trackMomentumAtVtx().R(),nBits_));
          patElectron_pOut.push_back(reduceFloat(iElectron.trackMomentumOut().R(),nBits_));
          patElectron_pAtCalo.push_back(reduceFloat(iElectron.trackMomentumAtCalo().R(),nBits_));
          patElectron_deltaEtaIn.push_back(reduceFloat(iElectron.deltaEtaSuperClusterTrackAtVtx(),nBits_));
          patElectron_deltaPhiIn.push_back(reduceFloat(iElectron.deltaPhiSuperClusterTrackAtVtx(),nBits_));
          patElectron_deltaEtaSeedClusterAtCalo.push_back(reduceFloat(iElectron.deltaEtaSeedClusterTrackAtCalo(),nBits_));
          patElectron_deltaEtaEleClusterAtCalo.push_back(reduceFloat(iElectron.deltaEtaEleClusterTrackAtCalo(),nBits_));
          patElectron_deltaPhiEleClusterAtCalo.push_back(reduceFloat(iElectron.deltaPhiEleClusterTrackAtCalo(),nBits_));
          patElectron_deltaPhiSeedClusterAtCalo.push_back(reduceFloat(iElectron.deltaPhiSeedClusterTrackAtCalo(),nBits_));
          patElectron_misHits.push_back(iElectron.gsfTrack()->hitPattern().numberOfAllTrackerHits(reco::HitPattern::TRACK_HITS));
          patElectron_nAmbiguousGsfTracks.push_back(iElectron.ambiguousGsfTracksSize());
          patElectron_trackFbrem.push_back(reduceFloat(iElectron.trackFbrem(),nBits_));
          patElectron_superClusterFbrem.push_back(reduceFloat(iElectron.superClusterFbrem(),nBits_));
          patElectron_dz.push_back(reduceFloat(iElectron.dB(pat::Electron::PVDZ),nBits_)); 
          patElectron_dzError.push_back(reduceFloat(abs(iElectron.edB(pat::Electron::PVDZ)),nBits_)); 
          patElectron_dxy.push_back(reduceFloat(iElectron.dB(pat::Electron::PV2D),nBits_)); 
          patElectron_dxyError.push_back(reduceFloat(abs(iElectron.edB(pat::Electron::PV2D)),nBits_)); 
          patElectron_energy.push_back(reduceFloat(iElectron.energy(),nBits_));
          patElectron_energyErr.push_back(reduceFloat(iElectron.p4Error(reco::GsfElectron::P4_COMBINATION),nBits_));
          patElectron_ecalEnergy.push_back(reduceFloat(iElectron.ecalEnergy(),nBits_));
          patElectron_ecalEnergyErr.push_back(reduceFloat(iElectron.ecalEnergyError(),nBits_));
          patElectron_et.push_back(reduceFloat(iElectron.p4().Et(),nBits_));
          patElectron_mt.push_back(reduceFloat(Mt,nBits_));
          patElectron_dphiMET.push_back(reduceFloat(dphiMET,nBits_));
          patElectron_scEta.push_back(reduceFloat(scRef->eta(),nBits_));
          patElectron_scPhi.push_back(reduceFloat(scRef->phi(),nBits_));
          patElectron_scEnergy.push_back(reduceFloat(scRef->energy(),nBits_)); 
          patElectron_scRawEnergy.push_back(reduceFloat(scRef->rawEnergy(),nBits_)); 
          patElectron_scRawESEnergy.push_back(reduceFloat(scRef->preshowerEnergy(),nBits_)); 
          patElectron_scEt.push_back(reduceFloat(scRef->energy()*(Rt/R),nBits_));
          patElectron_scPhiWidth.push_back(reduceFloat(scRef->phiWidth(),nBits_));
          patElectron_scEtaWidth.push_back(reduceFloat(scRef->etaWidth(),nBits_)); 
          patElectron_scEoP.push_back(reduceFloat(scRef->energy()/iElectron.trackMomentumAtVtx().R(),nBits_)); 
          if(isMatchedInEB || isMatchedInEE){
             patElectron_ecalSCEta.push_back(reduceFloat(matchedEcalSC.eta(),nBits_)); 
             patElectron_ecalSCPhi.push_back(reduceFloat(matchedEcalSC.phi(),nBits_)); 
             patElectron_ecalSCEnergy.push_back(reduceFloat(matchedEcalSC.energy(),nBits_)); 
             patElectron_ecalSCRawEnergy.push_back(reduceFloat(matchedEcalSC.rawEnergy(),nBits_)); 
             patElectron_ecalSCRawESEnergy.push_back(reduceFloat(matchedEcalSC.preshowerEnergy(),nBits_)); 
             patElectron_ecalSCEt.push_back(reduceFloat(matchedEcalSC.energy()*(ecalRt/ecalR),nBits_));
             patElectron_ecalSCPhiWidth.push_back(reduceFloat(matchedEcalSC.phiWidth(),nBits_));
             patElectron_ecalSCEtaWidth.push_back(reduceFloat(matchedEcalSC.etaWidth(),nBits_)); 
             patElectron_ecalSCEoP.push_back(reduceFloat(matchedEcalSC.energy()/iElectron.trackMomentumAtVtx().R(),nBits_));   
          }else{
             patElectron_ecalSCEta.push_back(-999.); 
             patElectron_ecalSCPhi.push_back(-999.); 
             patElectron_ecalSCEnergy.push_back(-999.); 
             patElectron_ecalSCRawEnergy.push_back(-999.); 
             patElectron_ecalSCRawESEnergy.push_back(-999.); 
             patElectron_ecalSCEt.push_back(-999.);
             patElectron_ecalSCPhiWidth.push_back(-999.);
             patElectron_ecalSCEtaWidth.push_back(-999.); 
             patElectron_ecalSCEoP.push_back(-999.);  
          }
          patElectron_scSwissCross.push_back(reduceFloat(swissCross,nBits_)); 
          patElectron_scE2x2.push_back(reduceFloat(e2x2,nBits_)); 
          patElectron_scE3x3.push_back(reduceFloat(e3x3,nBits_)); 
          patElectron_scE5x5.push_back(reduceFloat(e5x5,nBits_)); 
          patElectron_scEMax.push_back(reduceFloat(eMax,nBits_)); 
          patElectron_scR9.push_back(reduceFloat(e3x3/scRef->rawEnergy(),nBits_)); 
          patElectron_scSigmaIEtaIEta.push_back(reduceFloat(eleSS.sigmaIetaIeta,nBits_));
          patElectron_scSigmaIEtaIPhi.push_back(reduceFloat(eleSS.sigmaIetaIphi,nBits_));
          patElectron_scSigmaIPhiIPhi.push_back(reduceFloat(eleSS.sigmaIphiIphi,nBits_));  
          patElectron_full5x5_scE2x2.push_back(reduceFloat(full5x5_e2x2,nBits_)); 
          patElectron_full5x5_scE3x3.push_back(reduceFloat(full5x5_e3x3,nBits_)); 
          patElectron_full5x5_scE5x5.push_back(reduceFloat(full5x5_e5x5,nBits_)); 
          patElectron_full5x5_scEMax.push_back(reduceFloat(full5x5_eMax,nBits_));
          patElectron_full5x5_scR9.push_back(reduceFloat(full5x5_e3x3/scRef->rawEnergy(),nBits_));
          patElectron_full5x5_scSigmaIEtaIEta.push_back(reduceFloat(full5x5_eleSS.sigmaIetaIeta,nBits_));
          patElectron_full5x5_scSigmaIEtaIPhi.push_back(reduceFloat(full5x5_eleSS.sigmaIetaIphi,nBits_));
          patElectron_full5x5_scSigmaIPhiIPhi.push_back(reduceFloat(full5x5_eleSS.sigmaIphiIphi,nBits_));  
          patElectron_HoE.push_back(reduceFloat(iElectron.hadronicOverEm(),nBits_)); 
          patElectron_trkIso03.push_back(reduceFloat(iElectron.dr03TkSumPt(),nBits_));    
          patElectron_ecalIso03.push_back(reduceFloat(iElectron.dr03EcalRecHitSumEt(),nBits_));   
          patElectron_hcalIso03.push_back(reduceFloat(iElectron.dr03HcalTowerSumEt(),nBits_));  
          patElectron_trkIso04.push_back(reduceFloat(iElectron.dr04TkSumPt(),nBits_));    
          patElectron_ecalIso04.push_back(reduceFloat(iElectron.dr04EcalRecHitSumEt(),nBits_));   
          patElectron_hcalIso04.push_back(reduceFloat(iElectron.dr04HcalTowerSumEt(),nBits_));  
          patElectron_pfPhotonIso.push_back(reduceFloat(iElectron.pfIsolationVariables().sumPhotonEt,nBits_));      
          patElectron_pfChargedHadronIso.push_back(reduceFloat(iElectron.pfIsolationVariables().sumChargedHadronPt,nBits_));      
          patElectron_pfNeutralHadronIso.push_back(reduceFloat(iElectron.pfIsolationVariables().sumNeutralHadronEt,nBits_));
          patElectron_mva_Isolated.push_back(reduceFloat(iElectron.mva_Isolated(),nBits_));   
          patElectron_mva_e_pi.push_back(reduceFloat(iElectron.mva_e_pi(),nBits_));   
          patElectron_dnn_signal_Isolated.push_back(reduceFloat(iElectron.dnn_signal_Isolated(),nBits_));   
          patElectron_dnn_signal_nonIsolated.push_back(reduceFloat(iElectron.dnn_signal_nonIsolated(),nBits_));   
          patElectron_dnn_bkg_nonIsolated.push_back(reduceFloat(iElectron.dnn_bkg_nonIsolated(),nBits_));   
          patElectron_dnn_bkg_Tau.push_back(reduceFloat(iElectron.dnn_bkg_Tau(),nBits_));   
          patElectron_dnn_bkg_Photon.push_back(reduceFloat(iElectron.dnn_bkg_Photon(),nBits_));      
          patElectron_egmCutBasedElectronIDVeto.push_back(iElectron.electronID(egmCutBasedElectronIDVeto_.c_str()));
          patElectron_egmCutBasedElectronIDloose.push_back(iElectron.electronID(egmCutBasedElectronIDloose_.c_str()));
          patElectron_egmCutBasedElectronIDmedium.push_back(iElectron.electronID(egmCutBasedElectronIDmedium_.c_str()));
          patElectron_egmCutBasedElectronIDtight.push_back(iElectron.electronID(egmCutBasedElectronIDtight_.c_str()));
          patElectron_egmMVAElectronIDloose.push_back(iElectron.electronID(egmMVAElectronIDloose_.c_str()));
          patElectron_egmMVAElectronIDmedium.push_back(iElectron.electronID(egmMVAElectronIDmedium_.c_str()));
          patElectron_egmMVAElectronIDtight.push_back(iElectron.electronID(egmMVAElectronIDtight_.c_str()));
          patElectron_egmMVAElectronIDlooseNoIso.push_back(iElectron.electronID(egmMVAElectronIDlooseNoIso_.c_str()));
          patElectron_egmMVAElectronIDmediumNoIso.push_back(iElectron.electronID(egmMVAElectronIDmediumNoIso_.c_str()));
          patElectron_egmMVAElectronIDtightNoIso.push_back(iElectron.electronID(egmMVAElectronIDtightNoIso_.c_str()));
          patElectron_heepElectronID.push_back(iElectron.electronID(heepElectronID_.c_str()));
          
          iEle++; 
      }
   }    
   
   //save patPhoton info
   if(savePatPhotons_)
   {    
      int iPho=0;
      for(const auto& iPhoton : *(patPhoton.product())){ 
 
          reco::SuperClusterRef scRef = iPhoton.superCluster();
          reco::PhotonCoreRef phoCoreRef = iPhoton.photonCore();
          double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
          double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());

          reco::SuperCluster matchedEcalSC; 
          bool isMatchedInEB = false;  
          bool isMatchedInEE = false;  
          for(const auto& iSuperCluster : *(superClusterEB.product())){
              if(scRef->seed()->seed().rawId()==iSuperCluster.seed()->seed().rawId()){ 
                 matchedEcalSC = iSuperCluster;   
                 isMatchedInEB = true;
                 break;
              }     
          }    
          if(!isMatchedInEB){
             for(const auto& iSuperCluster : *(superClusterEE.product())){
                 if(scRef->seed()->seed().rawId()==iSuperCluster.seed()->seed().rawId()){ 
                    matchedEcalSC = iSuperCluster;  
                    isMatchedInEE = true; 
                    break;
                 }     
             }
          }    
          double ecalR  = TMath::Sqrt(matchedEcalSC.x()*matchedEcalSC.x() + matchedEcalSC.y()*matchedEcalSC.y() + matchedEcalSC.z()*matchedEcalSC.z());
          double ecalRt = TMath::Sqrt(matchedEcalSC.x()*matchedEcalSC.x() + matchedEcalSC.y()*matchedEcalSC.y());  

          double swissCross = -999.;
          double e3x3 = -999.;
          double e2x2 = -999.;
          double e5x5 = -999.;  
          double eMax = -999.;
          double full5x5_e3x3 = -999.;
          double full5x5_e2x2 = -999.;
          double full5x5_e5x5 = -999.;  
          double full5x5_eMax = -999.; 
         
          reco::Photon::ShowerShape phoSS = iPhoton.showerShapeVariables(); 
          reco::Photon::ShowerShape full5x5_phoSS = iPhoton.full5x5_showerShapeVariables();    
          const std::vector<std::pair<DetId,float> > &hits= iPhoton.superCluster()->hitsAndFractions();
          if(iPhoton.isEB())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEB.product())));   

             swissCross = EcalTools::swissCross(id.first,*(recHitsEB.product()),0.);
             e2x2 = EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             e5x5 = EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));
             full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));
          }
          if(iPhoton.isEE())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEE.product())));   

             swissCross = EcalTools::swissCross(id.first,*(recHitsEE.product()),0.);
             e2x2 = EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             e5x5 = EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
             full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
          }

          double cphi = (iPhoton.p4().x() * metP->corPx() + iPhoton.p4().y() * metP->corPy()) / (metP->corPt()*iPhoton.p4().Pt());
          double Mt = sqrt(2 * iPhoton.p4().Et() * metP->corPt() * (1-cphi));
          double dphiMET = deltaPhi(metP->corPhi(),iPhoton.p4().phi());  

          unsigned int nElectrons = 0;
          int electronIndex = -1;
          if(iPhoton.hasUserCand("electron")){ 
             nElectrons = 1;
             electronIndex = iPhoton.userCand("electron").key();
          }  

          unsigned int nJets = 0;
          int jetIndex = -1;
          if(iPhoton.hasUserCand("jet")){ 
             nJets = 1;
             jetIndex = iPhoton.userCand("jet").key();
          }  

          patPhoton_index.push_back(iPho);
          patPhoton_seedRawId.push_back(iPhoton.superCluster()->seed()->seed().rawId());
          patPhoton_scNPFClusters.push_back(scRef->clusters().size());
          if(isMatchedInEB || isMatchedInEE) patPhoton_ecalSCNPFClusters.push_back(matchedEcalSC.clusters().size());  
          else patPhoton_ecalSCNPFClusters.push_back(-1);   
          patPhoton_isEB.push_back(iPhoton.isEB());
          patPhoton_isEE.push_back(iPhoton.isEE());  
          patPhoton_isEBEEGap.push_back(iPhoton.isEBEEGap());  
          patPhoton_isEBEtaGap.push_back(iPhoton.isEBEtaGap());  
          patPhoton_isEBPhiGap.push_back(iPhoton.isEBPhiGap());  
          patPhoton_isEEDeeGap.push_back(iPhoton.isEEDeeGap());  
          patPhoton_isEERingGap.push_back(iPhoton.isEERingGap());   
          patPhoton_passElectronVeto.push_back(iPhoton.passElectronVeto()); 
          patPhoton_hasPixelSeed.push_back(iPhoton.hasPixelSeed());  
          patPhoton_hasConversionTracks.push_back(iPhoton.hasConversionTracks());
          patPhoton_nConversions.push_back(phoCoreRef->conversions().size());
          patPhoton_nConversionsOneLeg.push_back(phoCoreRef->conversionsOneLeg().size());
          patPhoton_hasOverlapElectron.push_back(nElectrons);  
          patPhoton_overlapElectronIndex.push_back(electronIndex);   
          patPhoton_hasOverlapJet.push_back(nJets);  
          patPhoton_overlapJetIndex.push_back(jetIndex);  
          patPhoton_eta.push_back(reduceFloat(iPhoton.p4().eta(),nBits_));
          patPhoton_phi.push_back(reduceFloat(iPhoton.p4().phi(),nBits_));
          patPhoton_energy.push_back(reduceFloat(iPhoton.energy(),nBits_));
          patPhoton_energyErr.push_back(reduceFloat(iPhoton.getCorrectedEnergyError(reco::Photon::regression2),nBits_));
          patPhoton_ecalEnergy.push_back(reduceFloat(iPhoton.energyCorrections().phoEcalEnergy,nBits_));
          patPhoton_ecalEnergyErr.push_back(reduceFloat(iPhoton.energyCorrections().phoEcalEnergyError,nBits_)); 
          patPhoton_et.push_back(reduceFloat(iPhoton.p4().Et(),nBits_));
          patPhoton_mt.push_back(reduceFloat(Mt,nBits_));
          patPhoton_dphiMET.push_back(reduceFloat(dphiMET,nBits_));
          patPhoton_scEta.push_back(reduceFloat(scRef->eta(),nBits_));
          patPhoton_scPhi.push_back(reduceFloat(scRef->phi(),nBits_));
          patPhoton_scEnergy.push_back(reduceFloat(scRef->energy(),nBits_)); 
          patPhoton_scRawEnergy.push_back(reduceFloat(scRef->rawEnergy(),nBits_)); 
          patPhoton_scRawESEnergy.push_back(reduceFloat(scRef->preshowerEnergy(),nBits_)); 
          patPhoton_scEt.push_back(reduceFloat(scRef->energy()*(Rt/R),nBits_));
          patPhoton_scPhiWidth.push_back(reduceFloat(scRef->phiWidth(),nBits_));
          patPhoton_scEtaWidth.push_back(reduceFloat(scRef->etaWidth(),nBits_)); 
          if(isMatchedInEB || isMatchedInEE){
             patPhoton_ecalSCEta.push_back(reduceFloat(matchedEcalSC.eta(),nBits_)); 
             patPhoton_ecalSCPhi.push_back(reduceFloat(matchedEcalSC.phi(),nBits_));  
             patPhoton_ecalSCEnergy.push_back(reduceFloat(matchedEcalSC.energy(),nBits_)); 
             patPhoton_ecalSCRawEnergy.push_back(reduceFloat(matchedEcalSC.rawEnergy(),nBits_)); 
             patPhoton_ecalSCRawESEnergy.push_back(reduceFloat(matchedEcalSC.preshowerEnergy(),nBits_)); 
             patPhoton_ecalSCEt.push_back(reduceFloat(matchedEcalSC.energy()*(ecalRt/ecalR),nBits_));
             patPhoton_ecalSCPhiWidth.push_back(reduceFloat(matchedEcalSC.phiWidth(),nBits_));
             patPhoton_ecalSCEtaWidth.push_back(reduceFloat(matchedEcalSC.etaWidth(),nBits_)); 
          }else{
             patPhoton_ecalSCEta.push_back(-999.); 
             patPhoton_ecalSCPhi.push_back(-999.);  
             patPhoton_ecalSCEnergy.push_back(-999.); 
             patPhoton_ecalSCRawEnergy.push_back(-999.); 
             patPhoton_ecalSCRawESEnergy.push_back(-999.); 
             patPhoton_ecalSCEt.push_back(-999.);
             patPhoton_ecalSCPhiWidth.push_back(-999.);
             patPhoton_ecalSCEtaWidth.push_back(-999.);   
          }    
          patPhoton_scSwissCross.push_back(reduceFloat(swissCross,nBits_)); 
          patPhoton_scE2x2.push_back(reduceFloat(e2x2,nBits_)); 
          patPhoton_scE3x3.push_back(reduceFloat(e3x3,nBits_)); 
          patPhoton_scE5x5.push_back(reduceFloat(e5x5,nBits_)); 
          patPhoton_scEMax.push_back(reduceFloat(eMax,nBits_)); 
          patPhoton_scR9.push_back(reduceFloat(e3x3/scRef->rawEnergy(),nBits_)); 
          patPhoton_scSigmaIEtaIEta.push_back(reduceFloat(phoSS.sigmaIetaIeta,nBits_));
          patPhoton_scSigmaIEtaIPhi.push_back(reduceFloat(phoSS.sigmaIetaIphi,nBits_));
          patPhoton_scSigmaIPhiIPhi.push_back(reduceFloat(phoSS.sigmaIphiIphi,nBits_));    
          patPhoton_full5x5_scEMax.push_back(reduceFloat(full5x5_eMax,nBits_)); 
          patPhoton_full5x5_scE2x2.push_back(reduceFloat(full5x5_e2x2,nBits_)); 
          patPhoton_full5x5_scE3x3.push_back(reduceFloat(full5x5_e3x3,nBits_)); 
          patPhoton_full5x5_scE5x5.push_back(reduceFloat(full5x5_e5x5,nBits_)); 
          patPhoton_full5x5_scEMax.push_back(reduceFloat(full5x5_eMax,nBits_));
          patPhoton_full5x5_scR9.push_back(reduceFloat(full5x5_e3x3/scRef->rawEnergy(),nBits_)); 
          patPhoton_full5x5_scSigmaIEtaIEta.push_back(reduceFloat(full5x5_phoSS.sigmaIetaIeta,nBits_));
          patPhoton_full5x5_scSigmaIEtaIPhi.push_back(reduceFloat(full5x5_phoSS.sigmaIetaIphi,nBits_));
          patPhoton_full5x5_scSigmaIPhiIPhi.push_back(reduceFloat(full5x5_phoSS.sigmaIphiIphi,nBits_));    
          patPhoton_HoE.push_back(reduceFloat(iPhoton.hadronicOverEm(),nBits_)); 
          patPhoton_trkIso03.push_back(reduceFloat(iPhoton.trkSumPtSolidConeDR03(),nBits_));
          patPhoton_ecalIso03.push_back(reduceFloat(iPhoton.ecalRecHitSumEtConeDR03(),nBits_));
          patPhoton_hcalIso03.push_back(reduceFloat(iPhoton.hcalTowerSumEtConeDR03(),nBits_));
          patPhoton_trkIso04.push_back(reduceFloat(iPhoton.trkSumPtSolidConeDR04(),nBits_));
          patPhoton_ecalIso04.push_back(reduceFloat(iPhoton.ecalRecHitSumEtConeDR04(),nBits_));
          patPhoton_hcalIso04.push_back(reduceFloat(iPhoton.hcalTowerSumEtConeDR04(),nBits_));
          patPhoton_patParticleIso.push_back(reduceFloat(iPhoton.patParticleIso(),nBits_));
          patPhoton_pfChargedHadronIso.push_back(reduceFloat(iPhoton.chargedHadronIso(),nBits_));
          patPhoton_pfNeutralHadronIso.push_back(reduceFloat(iPhoton.neutralHadronIso(),nBits_));
          patPhoton_pfPhotonIso.push_back(reduceFloat(iPhoton.photonIso(),nBits_));
          patPhoton_pfPuChargedHadronIso.push_back(reduceFloat(iPhoton.puChargedHadronIso(),nBits_));
          patPhoton_nClusterOutsideMustache.push_back(reduceFloat(iPhoton.nClusterOutsideMustache(),nBits_)); 
          patPhoton_etOutsideMustache.push_back(reduceFloat(iPhoton.etOutsideMustache(),nBits_)); 
          patPhoton_pfMVA.push_back(reduceFloat(iPhoton.pfMVA(),nBits_));  
          patPhoton_pfDNN.push_back(reduceFloat(iPhoton.pfDNN(),nBits_));  
          patPhoton_egmCutBasedPhotonIDloose.push_back(iPhoton.photonID(egmCutBasedPhotonIDloose_.c_str()));
          patPhoton_egmCutBasedPhotonIDmedium.push_back(iPhoton.photonID(egmCutBasedPhotonIDmedium_.c_str()));
          patPhoton_egmCutBasedPhotonIDtight.push_back(iPhoton.photonID(egmCutBasedPhotonIDtight_.c_str()));
          patPhoton_egmMVAPhotonIDmedium.push_back(iPhoton.photonID(egmMVAPhotonIDmedium_.c_str()));
          patPhoton_egmMVAPhotonIDtight.push_back(iPhoton.photonID(egmMVAPhotonIDtight_.c_str()));
          
          iPho++; 
      }
   }
 
   //save patJet info
   if(savePatJets_)
   {    
      int iJ=0;
      for(const auto& iJet : *(patJet.product())){ 
          
          int nCandInEcal = 0; 
          int nCandInEcalWithCharge = 0;      
          std::vector<float> charge;
          std::vector<float> ecalEnergies; 
          std::vector<float> ecalEnergiesFraction; 
          std::vector<float> hcalEnergies; 
          std::vector<float> hcalEnergiesFraction;  
          std::vector<float> eta;  
          std::vector<float> phi;    
          
          int nTotal = iJet.numberOfDaughters();
          for(unsigned ic = 0; ic < iJet.numberOfDaughters(); ic++){
  
              auto ip = dynamic_cast<const pat::PackedCandidate*> (iJet.daughter(ic));   
              
              double ecalEnergy = ip->energy()*ip->caloFraction()*(1-ip->hcalFraction());
              double ecalEnergyFraction = ip->caloFraction()*(1-ip->hcalFraction()); 
              double hcalEnergy = ip->energy()*ip->caloFraction()*ip->hcalFraction(); 
              double hcalEnergyFraction = ip->caloFraction()*ip->hcalFraction(); 

              if(ecalEnergy>0.){
                 nCandInEcal++;
                 if(ip->charge()){ 
                    nCandInEcalWithCharge++;
                    charge.push_back(reduceFloat(ip->charge(),nBits_));
                 }else{ 
                    charge.push_back(0.);
                 }
                 ecalEnergies.push_back(reduceFloat(ecalEnergy,nBits_));
                 ecalEnergiesFraction.push_back(reduceFloat(ecalEnergyFraction,nBits_));
                 hcalEnergies.push_back(reduceFloat(hcalEnergy,nBits_));
                 hcalEnergiesFraction.push_back(reduceFloat(hcalEnergyFraction,nBits_));
                 eta.push_back(reduceFloat(ip->eta(),nBits_));
                 phi.push_back(reduceFloat(ip->phi(),nBits_));
              }     
          }
          
          unsigned int nMuons = 0;
          if(iJet.hasOverlaps("muons")) nMuons = iJet.overlaps("muons").size();   

          unsigned int nTaus = 0;
          if(iJet.hasOverlaps("taus")) nTaus = iJet.overlaps("taus").size();    

          unsigned int nElectrons = 0;
          if(iJet.hasOverlaps("electrons")) nElectrons = iJet.overlaps("electrons").size(); 
          std::vector<int> electronIndices;
          electronIndices.resize(nElectrons);
          for(unsigned int iEle=0; iEle<nElectrons; iEle++)
              electronIndices[iEle] = iJet.overlaps("electrons")[iEle].key();

          unsigned int nPhotons = 0;
          if(iJet.hasOverlaps("photons")) nPhotons = iJet.overlaps("photons").size();  
          std::vector<int> photonIndices;
          photonIndices.resize(nPhotons);
          for(unsigned int iPho=0; iPho<nPhotons; iPho++)
              photonIndices[iPho] = iJet.overlaps("photons")[iPho].key();
       
          double btag_score = iJet.bDiscriminator("pfDeepCSVJetTags:probb") + iJet.bDiscriminator("pfDeepCSVJetTags:probbb");
          double puID_score = iJet.userFloat("pileupJetId:fullDiscriminant"); 
          int puID = iJet.userInt("pileupJetId:fullId");

          bool pfJetID_loose = passesJedID(&iJet,jmeCutBasedPFJetIDloose_.at(0),jmeCutBasedPFJetIDloose_.at(1),jmeCutBasedPFJetIDloose_.at(2));
          bool pfJetID_tight = passesJedID(&iJet,jmeCutBasedPFJetIDtight_.at(0),jmeCutBasedPFJetIDtight_.at(1),jmeCutBasedPFJetIDtight_.at(2));
          bool pfJetID_tightLepVeto = passesJedID(&iJet,jmeCutBasedPFJetIDtightLepVeto_.at(0),jmeCutBasedPFJetIDtightLepVeto_.at(1),jmeCutBasedPFJetIDtightLepVeto_.at(2));

          double qgLikelihood = iJet.userFloat("QGTagger:qgLikelihood");  

          patJet_index.push_back(iJ);
          patJet_isCaloJet.push_back(iJet.isCaloJet());  
          patJet_isJPTJet.push_back(iJet.isJPTJet());  
          patJet_isPFJet.push_back(iJet.isPFJet());  
          patJet_isBasicJet.push_back(iJet.isBasicJet());  
          patJet_charge.push_back(reduceFloat(iJet.jetCharge(),nBits_)); 
          patJet_energy.push_back(reduceFloat(iJet.energy(),nBits_)); 
          patJet_uncorrectedEnergy.push_back(reduceFloat(iJet.correctedP4("Uncorrected").energy(),nBits_)); 
          patJet_pt.push_back(reduceFloat(iJet.pt(),nBits_)); 
          patJet_uncorrectedPt.push_back(reduceFloat(iJet.correctedP4("Uncorrected").pt(),nBits_));      
          patJet_eta.push_back(reduceFloat(iJet.eta(),nBits_));    
          patJet_phi.push_back(reduceFloat(iJet.phi(),nBits_));   
          patJet_area.push_back(reduceFloat(iJet.jetArea(),nBits_));     
          if(iJet.isCaloJet()){
             patJet_energyFractionHadronic.push_back(reduceFloat(iJet.energyFractionHadronic(),nBits_));
             patJet_hadEnergyInHB.push_back(reduceFloat(iJet.hadEnergyInHB(),nBits_));
             patJet_hadEnergyInHO.push_back(reduceFloat(iJet.hadEnergyInHO(),nBits_));
             patJet_hadEnergyInHE.push_back(reduceFloat(iJet.hadEnergyInHE(),nBits_));
             patJet_hadEnergyInHF.push_back(reduceFloat(iJet.hadEnergyInHF(),nBits_));
             patJet_emEnergyInEB.push_back(reduceFloat(iJet.emEnergyInEB(),nBits_));
             patJet_emEnergyInEE.push_back(reduceFloat(iJet.emEnergyInEE(),nBits_));
             patJet_emEnergyInHF.push_back(reduceFloat(iJet.emEnergyInHF(),nBits_));
             patJet_chargedHadronEnergyFraction.push_back(reduceFloat(-999.,nBits_));
             patJet_neutralHadronEnergyFraction.push_back(reduceFloat(-999.,nBits_)); 
             patJet_chargedEmEnergyFraction.push_back(reduceFloat(-999.,nBits_));
             patJet_neutralEmEnergyFraction.push_back(reduceFloat(-999.,nBits_)); 
          }else{
             patJet_energyFractionHadronic.push_back(reduceFloat(-999.,nBits_));
             patJet_hadEnergyInHB.push_back(reduceFloat(-999.,nBits_));
             patJet_hadEnergyInHO.push_back(reduceFloat(-999.,nBits_));
             patJet_hadEnergyInHE.push_back(reduceFloat(-999.,nBits_));
             patJet_hadEnergyInHF.push_back(reduceFloat(-999.,nBits_));
             patJet_emEnergyInEB.push_back(reduceFloat(-999.,nBits_));
             patJet_emEnergyInEE.push_back(reduceFloat(-999.,nBits_));
             patJet_emEnergyInHF.push_back(reduceFloat(-999.,nBits_));   
             patJet_chargedHadronEnergyFraction.push_back(reduceFloat(iJet.chargedHadronEnergyFraction(),nBits_));
             patJet_neutralHadronEnergyFraction.push_back(reduceFloat(iJet.neutralHadronEnergyFraction(),nBits_)); 
             patJet_chargedEmEnergyFraction.push_back(reduceFloat(iJet.chargedEmEnergyFraction(),nBits_));
             patJet_neutralEmEnergyFraction.push_back(reduceFloat(iJet.neutralEmEnergyFraction(),nBits_));
          }
          patJet_nOverlapMuons.push_back(nMuons);    
          patJet_nOverlapTaus.push_back(nTaus);    
          patJet_nOverlapElectrons.push_back(nElectrons);  
          patJet_overlapElectronIndices.push_back(electronIndices);   
          patJet_nOverlapPhotons.push_back(nPhotons);  
          patJet_overlapPhotonIndices.push_back(photonIndices);      
          patJet_photonEnergy.push_back(reduceFloat(iJet.photonEnergy(),nBits_));
          patJet_photonEnergyFraction.push_back(reduceFloat(iJet.photonEnergyFraction(),nBits_));
          patJet_electronEnergy.push_back(reduceFloat(iJet.electronEnergy(),nBits_));
          patJet_electronEnergyFraction.push_back(reduceFloat(iJet.electronEnergyFraction(),nBits_));
          patJet_muonEnergy.push_back(reduceFloat(iJet.muonEnergy(),nBits_));
          patJet_muonEnergyFraction.push_back(reduceFloat(iJet.muonEnergyFraction(),nBits_));
          patJet_HFHadronEnergy.push_back(reduceFloat(iJet.HFHadronEnergy(),nBits_));          
          patJet_HFHadronEnergyFraction.push_back(reduceFloat(iJet.HFHadronEnergyFraction(),nBits_));
          patJet_HFEMEnergy.push_back(reduceFloat(iJet.HFEMEnergy(),nBits_));   
          patJet_HFEMEnergyFraction.push_back(reduceFloat(iJet.HFEMEnergyFraction(),nBits_));
          patJet_chargedMuEnergy.push_back(reduceFloat(iJet.chargedMuEnergy(),nBits_));   
          patJet_chargedMuEnergyFraction.push_back(reduceFloat(iJet.chargedMuEnergyFraction(),nBits_));
          patJet_hoEnergy.push_back(reduceFloat(iJet.hoEnergy(),nBits_));   
          patJet_hoEnergyFraction.push_back(reduceFloat(iJet.hoEnergyFraction(),nBits_));
          patJet_nCandidates.push_back(nTotal); 
          patJet_nCandInEcal.push_back(nCandInEcal); 
          patJet_nCandInEcalWithCharge.push_back(nCandInEcalWithCharge); 
          patJet_candInEcal_charge.push_back(charge); 
          patJet_candInEcal_ecalEnergy.push_back(ecalEnergies); 
          patJet_candInEcal_ecalEnergyFraction.push_back(ecalEnergiesFraction); 
          patJet_candInEcal_hcalEnergy.push_back(hcalEnergies); 
          patJet_candInEcal_hcalEnergyFraction.push_back(hcalEnergiesFraction); 
          patJet_candInEcal_eta.push_back(eta);  
          patJet_candInEcal_phi.push_back(phi);     
          patJet_bTagScore_pfDeepCSV.push_back(reduceFloat(btag_score,nBits_));
          patJet_puIDScore.push_back(reduceFloat(puID_score,nBits_));
          patJet_puID.push_back(puID);  
          patJet_qgL.push_back(reduceFloat(qgLikelihood,nBits_));   
          patJet_jmeCutBasedPFJetIDloose.push_back(pfJetID_loose);
          patJet_jmeCutBasedPFJetIDtight.push_back(pfJetID_tight); 
          patJet_jmeCutBasedPFJetIDtightLepVeto.push_back(pfJetID_tightLepVeto); 
        
          iJ++; 
      }  
   }

   //Save SuperClusters 
   if(saveSuperCluster_){
      int iSC=0;
      //std::cout << "SuperClustersEB size: " << (superClusterEB.product())->size() << std::endl;
      for(const auto& iSuperCluster : *(superClusterEB.product())){  

          superCluster_seedRawId.push_back(iSuperCluster.seed()->seed().rawId());
          superCluster_rawEnergy.push_back(reduceFloat(iSuperCluster.rawEnergy(),nBits_));
          superCluster_rawESEnergy.push_back(reduceFloat(iSuperCluster.preshowerEnergy(),nBits_));
          superCluster_energy.push_back(reduceFloat(iSuperCluster.energy(),nBits_));
          superCluster_eta.push_back(reduceFloat(iSuperCluster.eta(),nBits_));
          superCluster_phi.push_back(reduceFloat(iSuperCluster.phi(),nBits_));
          superCluster_etaWidth.push_back(reduceFloat(iSuperCluster.etaWidth(),nBits_));
          superCluster_phiWidth.push_back(reduceFloat(iSuperCluster.phiWidth(),nBits_));
          superCluster_R.push_back(reduceFloat(iSuperCluster.position().R(),nBits_));
          superCluster_nPFClusters.push_back(iSuperCluster.clusters().size());
          math::XYZPoint caloPos = iSuperCluster.seed()->position();
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          superCluster_ieta.push_back(eb_id.ieta());
          superCluster_iphi.push_back(eb_id.iphi());
          superCluster_iz.push_back(0);   
 
          if(saveShowerShapes_){
             reco::CaloCluster caloBC(*iSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), topology);  
             superCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             superCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             superCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     superCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             superCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             superCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             superCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             superCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             superCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             superCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             superCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             superCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             superCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             superCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             superCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             superCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iSuperCluster.rawEnergy(),nBits_));
             superCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             superCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             superCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             superCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             superCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             superCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             superCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             superCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             superCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             superCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             superCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             superCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             superCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             superCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             superCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             superCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             superCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             superCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             superCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iSuperCluster.rawEnergy(),nBits_));
             superCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             superCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             superCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 
          } 
         
          //compute scores  
          if(saveGenParticles_ && isMC_){
             dR_genScore.clear();
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             superCluster_dR_genScore[iSC] = dR_genScore;        
             superCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_dR_genScore, 999., false, 0., iSC));
          } 

          if(saveCaloParticlesPU_ && isMC_){ 
             simPU_nSharedXtals.clear();
             simEnergy_sharedXtalsPU.clear(); 
             recoEnergy_sharedXtalsPU.clear();
             simEnergy_noHitsFraction_sharedXtalsPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iSuperCluster,&hitsAndEnergies_CaloPartPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[8]);   
             } 
             superCluster_simPU_nSharedXtals[iSC] = simPU_nSharedXtals[0];   
             superCluster_simEnergy_sharedXtalsPU[iSC] = simEnergy_sharedXtalsPU[0];   
             superCluster_recoEnergy_sharedXtalsPU[iSC] = recoEnergy_sharedXtalsPU[0];     
             superCluster_simEnergy_noHitsFraction_sharedXtalsPU[iSC] = simEnergy_noHitsFraction_sharedXtalsPU[0];   
             superCluster_recoEnergy_noHitsFraction_sharedXtalsPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsPU[0];        
          }

          if(saveCaloParticlesOOTPU_ && isMC_){ 
             simOOTPU_nSharedXtals.clear();
             simEnergy_sharedXtalsOOTPU.clear(); 
             recoEnergy_sharedXtalsOOTPU.clear();
             simEnergy_noHitsFraction_sharedXtalsOOTPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsOOTPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartOOTPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iSuperCluster,&hitsAndEnergies_CaloPartOOTPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simOOTPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsOOTPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsOOTPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[8]);   
             } 
             superCluster_simOOTPU_nSharedXtals[iSC] = simOOTPU_nSharedXtals[0];   
             superCluster_simEnergy_sharedXtalsOOTPU[iSC] = simEnergy_sharedXtalsOOTPU[0];   
             superCluster_recoEnergy_sharedXtalsOOTPU[iSC] = recoEnergy_sharedXtalsOOTPU[0];     
             superCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = simEnergy_noHitsFraction_sharedXtalsOOTPU[0];   
             superCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsOOTPU[0];        
          }

          if(saveCaloParticles_ && isMC_){ 
             dR_simScore.clear();
             sim_nSharedXtals.clear();
             sim_fraction_noHitsFraction.clear();
             sim_fraction.clear();
             recoToSim_fraction.clear();
             recoToSim_fraction_sharedXtals.clear();  
             simEnergy_sharedXtals.clear(); 
             recoEnergy_sharedXtals.clear(); 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 GlobalPoint caloParticle_position = caloParts_position.at(iCalo);
                 std::vector<double> scores = getScores(&iSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             superCluster_nXtals.push_back((iSuperCluster.hitsAndFractions()).size());   
             superCluster_dR_simScore[iSC] = dR_simScore;  
             superCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             superCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             superCluster_sim_fraction[iSC] = sim_fraction;   
             superCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             superCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             superCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             superCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             superCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_dR_simScore, 999., false, 0., iSC));
             superCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_nSharedXtals, -999., true, 0., iSC));
             superCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             superCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction, -999., true, 0., iSC));
             superCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             superCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             superCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             superCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoEnergy_sharedXtals, -999., true, 0., iSC)); 
          }
          
          if(savePFCluster_ && saveSuperCluster_){   
             //save clusters and superClusters mutual info
             reco::CaloCluster caloSeed(*iSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC)); 
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) superCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) superCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      } 

      // The global SuperCluster indexing for EE has an offset = nSuperClusterEB
      iSC = (superClusterEB.product())->size();
      int iSC_tmp=-1;
      //std::cout << "SuperClustersEE size: " << (superClusterEE.product())->size() << std::endl;

      for(const auto& iSuperCluster : *(superClusterEE.product())){    

          iSC_tmp++;
        
          superCluster_seedRawId.push_back(iSuperCluster.seed()->seed().rawId());
          superCluster_rawEnergy.push_back(reduceFloat(iSuperCluster.rawEnergy(),nBits_));
          superCluster_rawESEnergy.push_back(reduceFloat(iSuperCluster.preshowerEnergy(),nBits_));
          superCluster_energy.push_back(reduceFloat(iSuperCluster.energy(),nBits_));
          superCluster_eta.push_back(reduceFloat(iSuperCluster.eta(),nBits_));
          superCluster_phi.push_back(reduceFloat(iSuperCluster.phi(),nBits_));
          superCluster_etaWidth.push_back(reduceFloat(iSuperCluster.etaWidth(),nBits_));
          superCluster_phiWidth.push_back(reduceFloat(iSuperCluster.phiWidth(),nBits_));
          superCluster_R.push_back(reduceFloat(iSuperCluster.position().R(),nBits_)); 
          superCluster_nPFClusters.push_back(iSuperCluster.clusters().size());  
          math::XYZPoint caloPos = iSuperCluster.seed()->position(); 
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          superCluster_ieta.push_back(ee_id.ix());
          superCluster_iphi.push_back(ee_id.iy());
          superCluster_iz.push_back(ee_id.zside());   

          if(saveShowerShapes_){ 
             reco::CaloCluster caloBC(*iSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), topology);  
             superCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             superCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             superCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     superCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             superCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             superCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             superCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             superCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             superCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             superCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             superCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             superCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             superCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             superCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             superCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             superCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iSuperCluster.rawEnergy(),nBits_));
             superCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             superCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             superCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             superCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             superCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             superCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             superCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             superCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             superCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             superCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             superCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             superCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             superCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             superCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             superCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             superCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             superCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             superCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             superCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iSuperCluster.rawEnergy(),nBits_));
             superCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             superCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             superCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 
          }

          //compute scores  
          if(saveGenParticles_ && isMC_){
             dR_genScore.clear();
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             superCluster_dR_genScore[iSC] = dR_genScore;        
             superCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_dR_genScore, 999., false, 0., iSC));
          } 
          
          if(saveCaloParticlesPU_ && isMC_){ 
             simPU_nSharedXtals.clear();
             simEnergy_sharedXtalsPU.clear(); 
             recoEnergy_sharedXtalsPU.clear();
             simEnergy_noHitsFraction_sharedXtalsPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iSuperCluster,&hitsAndEnergies_CaloPartPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[8]);   
             } 
             superCluster_simPU_nSharedXtals[iSC] = simPU_nSharedXtals[0];   
             superCluster_simEnergy_sharedXtalsPU[iSC] = simEnergy_sharedXtalsPU[0];   
             superCluster_recoEnergy_sharedXtalsPU[iSC] = recoEnergy_sharedXtalsPU[0];     
             superCluster_simEnergy_noHitsFraction_sharedXtalsPU[iSC] = simEnergy_noHitsFraction_sharedXtalsPU[0];   
             superCluster_recoEnergy_noHitsFraction_sharedXtalsPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsPU[0];        
          }

          if(saveCaloParticlesOOTPU_ && isMC_){ 
             simOOTPU_nSharedXtals.clear();
             simEnergy_sharedXtalsOOTPU.clear(); 
             recoEnergy_sharedXtalsOOTPU.clear();
             simEnergy_noHitsFraction_sharedXtalsOOTPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsOOTPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartOOTPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iSuperCluster,&hitsAndEnergies_CaloPartOOTPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simOOTPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsOOTPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsOOTPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[8]);   
             } 
             superCluster_simOOTPU_nSharedXtals[iSC] = simOOTPU_nSharedXtals[0];   
             superCluster_simEnergy_sharedXtalsOOTPU[iSC] = simEnergy_sharedXtalsOOTPU[0];   
             superCluster_recoEnergy_sharedXtalsOOTPU[iSC] = recoEnergy_sharedXtalsOOTPU[0];     
             superCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = simEnergy_noHitsFraction_sharedXtalsOOTPU[0];   
             superCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsOOTPU[0];        
          }
 
          if(saveCaloParticles_ && isMC_){ 
             dR_simScore.clear();
             sim_nSharedXtals.clear();
             sim_fraction_noHitsFraction.clear();
             sim_fraction.clear();
             recoToSim_fraction.clear();
             recoToSim_fraction_sharedXtals.clear();  
             simEnergy_sharedXtals.clear(); 
             recoEnergy_sharedXtals.clear();  
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 GlobalPoint caloParticle_position = caloParts_position.at(iCalo);
                 std::vector<double> scores = getScores(&iSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             superCluster_nXtals.push_back((iSuperCluster.hitsAndFractions()).size());   
             superCluster_dR_simScore[iSC] = dR_simScore;  
             superCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             superCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             superCluster_sim_fraction[iSC] = sim_fraction;   
             superCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             superCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             superCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             superCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             superCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_dR_simScore, 999., false, 0., iSC));
             superCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_nSharedXtals, -999., true, 0., iSC));
             superCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             superCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction, -999., true, 0., iSC));
             superCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             superCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             superCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             superCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoEnergy_sharedXtals, -999., true, 0., iSC)); 
          }
          
          if(iSuperCluster.preshowerClusters().isAvailable()){
              for(unsigned int iPC=0; iPC<iSuperCluster.preshowerClusters().size(); iPC++){
                  if(!iSuperCluster.preshowerClusters()[iPC].isAvailable()) { continue; } 
                  superCluster_psCluster_energy[iSC_tmp].push_back(reduceFloat(iSuperCluster.preshowerClusters()[iPC]->energy(),nBits_));
                  superCluster_psCluster_eta[iSC_tmp].push_back(reduceFloat(iSuperCluster.preshowerClusters()[iPC]->eta(),nBits_));
                  superCluster_psCluster_phi[iSC_tmp].push_back(reduceFloat(iSuperCluster.preshowerClusters()[iPC]->phi(),nBits_));   
              }
          } 

          if(savePFCluster_ && saveSuperCluster_){   
             //save clusters and superClusters mutual info
             reco::CaloCluster caloSeed(*iSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC));  
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) superCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) superCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      }

      //save pfCluster_superClustersIndex
      if(savePFCluster_ && saveSuperCluster_){
         for(unsigned int iSC=0; iSC<superCluster_pfClustersIndex.size(); iSC++)
             for(unsigned int iPF=0; iPF<superCluster_pfClustersIndex.at(iSC).size(); iPF++)
                 if(superCluster_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_superClustersIndex[superCluster_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
      }

     //save inverse of matchings
     if(saveGenParticles_ && isMC_){ 
        fillParticleMatchedIndex(&genParticle_superCluster_dR_genScore_MatchedIndex,&superCluster_dR_genScore_MatchedIndex);
      } 
      if(saveCaloParticles_ && isMC_){ 
         fillParticleMatchedIndex(&caloParticle_superCluster_dR_simScore_MatchedIndex,&superCluster_dR_simScore_MatchedIndex);
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_nSharedXtals_MatchedIndex,&superCluster_sim_nSharedXtals_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_noHitsFraction_MatchedIndex,&superCluster_sim_fraction_noHitsFraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_MatchedIndex,&superCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_recoToSim_fraction_MatchedIndex,&superCluster_recoToSim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_recoToSim_fraction_sharedXtals_MatchedIndex,&superCluster_recoToSim_fraction_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_simEnergy_sharedXtals_MatchedIndex,&superCluster_simEnergy_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_recoEnergy_sharedXtals_MatchedIndex,&superCluster_recoEnergy_sharedXtals_MatchedIndex);  
      }
   }

   //Save retunedSuperClusters 
   //std::cout << "-----> retunedSuperClusters <-----" << std::endl;  
   if(saveRetunedSC_){
      int iSC=0;
      //std::cout << "retunedSuperClustersEB size: " << (retunedSuperClusterEB.product())->size() << std::endl;
      for(const auto& iRetunedSuperCluster : *(retunedSuperClusterEB.product())){  

          retunedSuperCluster_seedRawId.push_back(iRetunedSuperCluster.seed()->seed().rawId());   
          retunedSuperCluster_rawEnergy.push_back(reduceFloat(iRetunedSuperCluster.rawEnergy(),nBits_));
          retunedSuperCluster_rawESEnergy.push_back(reduceFloat(iRetunedSuperCluster.preshowerEnergy(),nBits_));
          retunedSuperCluster_energy.push_back(reduceFloat(iRetunedSuperCluster.energy(),nBits_)); 
          retunedSuperCluster_eta.push_back(reduceFloat(iRetunedSuperCluster.eta(),nBits_));
          retunedSuperCluster_phi.push_back(reduceFloat(iRetunedSuperCluster.phi(),nBits_));
          retunedSuperCluster_etaWidth.push_back(reduceFloat(iRetunedSuperCluster.etaWidth(),nBits_));
          retunedSuperCluster_phiWidth.push_back(reduceFloat(iRetunedSuperCluster.phiWidth(),nBits_));
          retunedSuperCluster_R.push_back(reduceFloat(iRetunedSuperCluster.position().R(),nBits_));
          retunedSuperCluster_nPFClusters.push_back(iRetunedSuperCluster.clusters().size());
          math::XYZPoint caloPos = iRetunedSuperCluster.seed()->position();
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          retunedSuperCluster_ieta.push_back(eb_id.ieta());
          retunedSuperCluster_iphi.push_back(eb_id.iphi());
          retunedSuperCluster_iz.push_back(0);   
 
          if(saveShowerShapes_){
             reco::CaloCluster caloBC(*iRetunedSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), topology);  
             retunedSuperCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             retunedSuperCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             retunedSuperCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     retunedSuperCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             retunedSuperCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             retunedSuperCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             retunedSuperCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             retunedSuperCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             retunedSuperCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             retunedSuperCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             retunedSuperCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             retunedSuperCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             retunedSuperCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             retunedSuperCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             retunedSuperCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             retunedSuperCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iRetunedSuperCluster.rawEnergy(),nBits_));
             retunedSuperCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             retunedSuperCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             retunedSuperCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             retunedSuperCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             retunedSuperCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             retunedSuperCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             retunedSuperCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             retunedSuperCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             retunedSuperCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             retunedSuperCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             retunedSuperCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             retunedSuperCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             retunedSuperCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             retunedSuperCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             retunedSuperCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             retunedSuperCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             retunedSuperCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             retunedSuperCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             retunedSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iRetunedSuperCluster.rawEnergy(),nBits_));
             retunedSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 
          } 
         
          //compute scores  
          if(saveGenParticles_ && isMC_){
             dR_genScore.clear();
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             retunedSuperCluster_dR_genScore[iSC] = dR_genScore;        
             retunedSuperCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_dR_genScore, 999., false, 0., iSC));
          } 

          if(saveCaloParticlesPU_ && isMC_){ 
             simPU_nSharedXtals.clear();
             simEnergy_sharedXtalsPU.clear(); 
             recoEnergy_sharedXtalsPU.clear();
             simEnergy_noHitsFraction_sharedXtalsPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iRetunedSuperCluster,&hitsAndEnergies_CaloPartPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[8]);   
             } 
             retunedSuperCluster_simPU_nSharedXtals[iSC] = simPU_nSharedXtals[0];   
             retunedSuperCluster_simEnergy_sharedXtalsPU[iSC] = simEnergy_sharedXtalsPU[0];   
             retunedSuperCluster_recoEnergy_sharedXtalsPU[iSC] = recoEnergy_sharedXtalsPU[0];     
             retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU[iSC] = simEnergy_noHitsFraction_sharedXtalsPU[0];   
             retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsPU[0];        
          }

          if(saveCaloParticlesOOTPU_ && isMC_){ 
             simOOTPU_nSharedXtals.clear();
             simEnergy_sharedXtalsOOTPU.clear(); 
             recoEnergy_sharedXtalsOOTPU.clear();
             simEnergy_noHitsFraction_sharedXtalsOOTPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsOOTPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartOOTPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iRetunedSuperCluster,&hitsAndEnergies_CaloPartOOTPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simOOTPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsOOTPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsOOTPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[8]);   
             } 
             retunedSuperCluster_simOOTPU_nSharedXtals[iSC] = simOOTPU_nSharedXtals[0];   
             retunedSuperCluster_simEnergy_sharedXtalsOOTPU[iSC] = simEnergy_sharedXtalsOOTPU[0];   
             retunedSuperCluster_recoEnergy_sharedXtalsOOTPU[iSC] = recoEnergy_sharedXtalsOOTPU[0];     
             retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = simEnergy_noHitsFraction_sharedXtalsOOTPU[0];   
             retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsOOTPU[0];        
          }

          if(saveCaloParticles_ && isMC_){ 
             dR_simScore.clear();
             sim_nSharedXtals.clear();
             sim_fraction_noHitsFraction.clear();
             sim_fraction.clear();
             recoToSim_fraction.clear();
             recoToSim_fraction_sharedXtals.clear();  
             simEnergy_sharedXtals.clear(); 
             recoEnergy_sharedXtals.clear(); 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 GlobalPoint caloParticle_position = caloParts_position.at(iCalo);
                 std::vector<double> scores = getScores(&iRetunedSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             retunedSuperCluster_nXtals.push_back((iRetunedSuperCluster.hitsAndFractions()).size());   
             retunedSuperCluster_dR_simScore[iSC] = dR_simScore;  
             retunedSuperCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             retunedSuperCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             retunedSuperCluster_sim_fraction[iSC] = sim_fraction;   
             retunedSuperCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             retunedSuperCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             retunedSuperCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             retunedSuperCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             retunedSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_dR_simScore, 999., false, 0., iSC));
             retunedSuperCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_nSharedXtals, -999., true, 0., iSC));
             retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             retunedSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction, -999., true, 0., iSC));
             retunedSuperCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoEnergy_sharedXtals, -999., true, 0., iSC));                
          }

          if(savePFCluster_){   
             //save clusters and retunedSuperClusters mutual info
             reco::CaloCluster caloSeed(*iRetunedSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iRetunedSuperCluster.clustersBegin(); iBC != iRetunedSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC)); 
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) retunedSuperCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) retunedSuperCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      } 

      // The global retunedSuperCluster indexing for EE has an offset = nretunedSuperClusterEB
      iSC = (retunedSuperClusterEB.product())->size();
      int iSC_tmp=-1;
      //std::cout << "retunedSuperClustersEE size: " << (retunedSuperClusterEE.product())->size() << std::endl;
      for(const auto& iRetunedSuperCluster : *(retunedSuperClusterEE.product())){    

          iSC_tmp++;
        
          retunedSuperCluster_seedRawId.push_back(iRetunedSuperCluster.seed()->seed().rawId());
          retunedSuperCluster_rawEnergy.push_back(reduceFloat(iRetunedSuperCluster.rawEnergy(),nBits_));
          retunedSuperCluster_rawESEnergy.push_back(reduceFloat(iRetunedSuperCluster.preshowerEnergy(),nBits_));
          retunedSuperCluster_energy.push_back(reduceFloat(iRetunedSuperCluster.energy(),nBits_));
          retunedSuperCluster_eta.push_back(reduceFloat(iRetunedSuperCluster.eta(),nBits_));
          retunedSuperCluster_phi.push_back(reduceFloat(iRetunedSuperCluster.phi(),nBits_));
          retunedSuperCluster_etaWidth.push_back(reduceFloat(iRetunedSuperCluster.etaWidth(),nBits_));
          retunedSuperCluster_phiWidth.push_back(reduceFloat(iRetunedSuperCluster.phiWidth(),nBits_));
          retunedSuperCluster_R.push_back(reduceFloat(iRetunedSuperCluster.position().R(),nBits_)); 
          retunedSuperCluster_nPFClusters.push_back(iRetunedSuperCluster.clusters().size());  
          math::XYZPoint caloPos = iRetunedSuperCluster.seed()->position(); 
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          retunedSuperCluster_ieta.push_back(ee_id.ix());
          retunedSuperCluster_iphi.push_back(ee_id.iy());
          retunedSuperCluster_iz.push_back(ee_id.zside());   

          if(saveShowerShapes_){ 
             reco::CaloCluster caloBC(*iRetunedSuperCluster.seed()); 
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), topology);  
             retunedSuperCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             retunedSuperCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             retunedSuperCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     retunedSuperCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             retunedSuperCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             retunedSuperCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             retunedSuperCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             retunedSuperCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             retunedSuperCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             retunedSuperCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             retunedSuperCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             retunedSuperCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             retunedSuperCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             retunedSuperCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             retunedSuperCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             retunedSuperCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iRetunedSuperCluster.rawEnergy(),nBits_));
             retunedSuperCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             retunedSuperCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             retunedSuperCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             retunedSuperCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             retunedSuperCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             retunedSuperCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             retunedSuperCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             retunedSuperCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             retunedSuperCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             retunedSuperCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             retunedSuperCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             retunedSuperCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             retunedSuperCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             retunedSuperCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             retunedSuperCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             retunedSuperCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             retunedSuperCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             retunedSuperCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             retunedSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iRetunedSuperCluster.rawEnergy(),nBits_));
             retunedSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 
          }

          //compute scores  
          if(saveGenParticles_ && isMC_){
             dR_genScore.clear();
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             retunedSuperCluster_dR_genScore[iSC] = dR_genScore;        
             retunedSuperCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_dR_genScore, 999., false, 0., iSC));
          } 

          if(saveCaloParticlesPU_ && isMC_){ 
             simPU_nSharedXtals.clear();
             simEnergy_sharedXtalsPU.clear(); 
             recoEnergy_sharedXtalsPU.clear();
             simEnergy_noHitsFraction_sharedXtalsPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iRetunedSuperCluster,&hitsAndEnergies_CaloPartPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[8]);   
             } 
             retunedSuperCluster_simPU_nSharedXtals[iSC] = simPU_nSharedXtals[0];   
             retunedSuperCluster_simEnergy_sharedXtalsPU[iSC] = simEnergy_sharedXtalsPU[0];   
             retunedSuperCluster_recoEnergy_sharedXtalsPU[iSC] = recoEnergy_sharedXtalsPU[0];     
             retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU[iSC] = simEnergy_noHitsFraction_sharedXtalsPU[0];   
             retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsPU[0];        
          }

          if(saveCaloParticlesOOTPU_ && isMC_){ 
             simOOTPU_nSharedXtals.clear();
             simEnergy_sharedXtalsOOTPU.clear(); 
             recoEnergy_sharedXtalsOOTPU.clear();
             simEnergy_noHitsFraction_sharedXtalsOOTPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsOOTPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartOOTPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iRetunedSuperCluster,&hitsAndEnergies_CaloPartOOTPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simOOTPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsOOTPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsOOTPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[8]);   
             } 
             retunedSuperCluster_simOOTPU_nSharedXtals[iSC] = simOOTPU_nSharedXtals[0];   
             retunedSuperCluster_simEnergy_sharedXtalsOOTPU[iSC] = simEnergy_sharedXtalsOOTPU[0];   
             retunedSuperCluster_recoEnergy_sharedXtalsOOTPU[iSC] = recoEnergy_sharedXtalsOOTPU[0];     
             retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = simEnergy_noHitsFraction_sharedXtalsOOTPU[0];   
             retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsOOTPU[0];        
          }

          if(saveCaloParticles_ && isMC_){ 
             dR_simScore.clear();
             sim_nSharedXtals.clear();
             sim_fraction_noHitsFraction.clear();
             sim_fraction.clear();
             recoToSim_fraction.clear();
             recoToSim_fraction_sharedXtals.clear();  
             simEnergy_sharedXtals.clear(); 
             recoEnergy_sharedXtals.clear(); 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 GlobalPoint caloParticle_position = caloParts_position.at(iCalo);
                 std::vector<double> scores = getScores(&iRetunedSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             retunedSuperCluster_nXtals.push_back((iRetunedSuperCluster.hitsAndFractions()).size());   
             retunedSuperCluster_dR_simScore[iSC] = dR_simScore;  
             retunedSuperCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             retunedSuperCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             retunedSuperCluster_sim_fraction[iSC] = sim_fraction;   
             retunedSuperCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             retunedSuperCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             retunedSuperCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             retunedSuperCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             retunedSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_dR_simScore, 999., false, 0., iSC));
             retunedSuperCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_nSharedXtals, -999., true, 0., iSC));
             retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             retunedSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction, -999., true, 0., iSC));
             retunedSuperCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoEnergy_sharedXtals, -999., true, 0., iSC)); 
          }

          if(iRetunedSuperCluster.preshowerClusters().isAvailable()){
              for(unsigned int iPC=0; iPC<iRetunedSuperCluster.preshowerClusters().size(); iPC++){
                  if(!iRetunedSuperCluster.preshowerClusters()[iPC].isAvailable()) { continue; } 
                  retunedSuperCluster_psCluster_energy[iSC_tmp].push_back(reduceFloat(iRetunedSuperCluster.preshowerClusters()[iPC]->energy(),nBits_));
                  retunedSuperCluster_psCluster_eta[iSC_tmp].push_back(reduceFloat(iRetunedSuperCluster.preshowerClusters()[iPC]->eta(),nBits_));
                  retunedSuperCluster_psCluster_phi[iSC_tmp].push_back(reduceFloat(iRetunedSuperCluster.preshowerClusters()[iPC]->phi(),nBits_));   
              }
          } 

          if(savePFCluster_ && saveRetunedSC_){   
             //save clusters and retunedSuperClusters mutual info
             reco::CaloCluster caloSeed(*iRetunedSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iRetunedSuperCluster.clustersBegin(); iBC != iRetunedSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC));  
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) retunedSuperCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) retunedSuperCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      }

      //save pfCluster_retunedSuperClustersIndex
      if(savePFCluster_ && saveRetunedSC_){
         for(unsigned int iSC=0; iSC<retunedSuperCluster_pfClustersIndex.size(); iSC++)
             for(unsigned int iPF=0; iPF<retunedSuperCluster_pfClustersIndex.at(iSC).size(); iPF++)
                 if(retunedSuperCluster_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_retunedSuperClustersIndex[retunedSuperCluster_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
      } 
      
      //save inverse of matchings
      if(saveGenParticles_ && isMC_){ 
         fillParticleMatchedIndex(&genParticle_retunedSuperCluster_dR_genScore_MatchedIndex,&retunedSuperCluster_dR_genScore_MatchedIndex);
      } 
      if(saveCaloParticles_ && isMC_){ 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex,&retunedSuperCluster_dR_simScore_MatchedIndex);
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_nSharedXtals_MatchedIndex,&retunedSuperCluster_sim_nSharedXtals_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex,&retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex,&retunedSuperCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_recoToSim_fraction_MatchedIndex,&retunedSuperCluster_recoToSim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex,&retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex,&retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex,&retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex);  
      } 
   }

   //Save deepSuperClusters 
   //std::cout << "-----> deepSuperClusters <-----" << std::endl;  
   if(saveDeepSC_){
      int iSC=0;
      //std::cout << "deepSuperClustersEB size: " << (deepSuperClusterEB.product())->size() << std::endl;
      for(const auto& iDeepSuperCluster : *(deepSuperClusterEB.product())){  

          deepSuperCluster_seedRawId.push_back(iDeepSuperCluster.seed()->seed().rawId());
          deepSuperCluster_rawEnergy.push_back(reduceFloat(iDeepSuperCluster.rawEnergy(),nBits_));
          deepSuperCluster_rawESEnergy.push_back(reduceFloat(iDeepSuperCluster.preshowerEnergy(),nBits_));
          deepSuperCluster_energy.push_back(reduceFloat(iDeepSuperCluster.energy(),nBits_));
          deepSuperCluster_eta.push_back(reduceFloat(iDeepSuperCluster.eta(),nBits_));
          deepSuperCluster_phi.push_back(reduceFloat(iDeepSuperCluster.phi(),nBits_));
          deepSuperCluster_etaWidth.push_back(reduceFloat(iDeepSuperCluster.etaWidth(),nBits_));
          deepSuperCluster_phiWidth.push_back(reduceFloat(iDeepSuperCluster.phiWidth(),nBits_));
          deepSuperCluster_R.push_back(reduceFloat(iDeepSuperCluster.position().R(),nBits_));
          deepSuperCluster_nPFClusters.push_back(iDeepSuperCluster.clusters().size());
          math::XYZPoint caloPos = iDeepSuperCluster.seed()->position();
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          deepSuperCluster_ieta.push_back(eb_id.ieta());
          deepSuperCluster_iphi.push_back(eb_id.iphi());
          deepSuperCluster_iz.push_back(0);   
 
          if(saveShowerShapes_){
             reco::CaloCluster caloBC(*iDeepSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), topology);  
             deepSuperCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             deepSuperCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             deepSuperCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     deepSuperCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             deepSuperCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             deepSuperCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             deepSuperCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             deepSuperCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             deepSuperCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             deepSuperCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             deepSuperCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             deepSuperCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             deepSuperCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             deepSuperCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             deepSuperCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             deepSuperCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iDeepSuperCluster.rawEnergy(),nBits_));
             deepSuperCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             deepSuperCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             deepSuperCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             deepSuperCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             deepSuperCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             deepSuperCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             deepSuperCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             deepSuperCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             deepSuperCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             deepSuperCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             deepSuperCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             deepSuperCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             deepSuperCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             deepSuperCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             deepSuperCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             deepSuperCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             deepSuperCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             deepSuperCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             deepSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iDeepSuperCluster.rawEnergy(),nBits_));
             deepSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             deepSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             deepSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 
          } 
         
          //compute scores  
          if(saveGenParticles_ && isMC_){
             dR_genScore.clear();
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             deepSuperCluster_dR_genScore[iSC] = dR_genScore;        
             deepSuperCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_dR_genScore, 999., false, 0., iSC));
          } 

          if(saveCaloParticlesPU_ && isMC_){ 
             simPU_nSharedXtals.clear();
             simEnergy_sharedXtalsPU.clear(); 
             recoEnergy_sharedXtalsPU.clear();
             simEnergy_noHitsFraction_sharedXtalsPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iDeepSuperCluster,&hitsAndEnergies_CaloPartPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[8]);   
             } 
             deepSuperCluster_simPU_nSharedXtals[iSC] = simPU_nSharedXtals[0];   
             deepSuperCluster_simEnergy_sharedXtalsPU[iSC] = simEnergy_sharedXtalsPU[0];   
             deepSuperCluster_recoEnergy_sharedXtalsPU[iSC] = recoEnergy_sharedXtalsPU[0];     
             deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU[iSC] = simEnergy_noHitsFraction_sharedXtalsPU[0];   
             deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsPU[0];        
          }

          if(saveCaloParticlesOOTPU_ && isMC_){ 
             simOOTPU_nSharedXtals.clear();
             simEnergy_sharedXtalsOOTPU.clear(); 
             recoEnergy_sharedXtalsOOTPU.clear();
             simEnergy_noHitsFraction_sharedXtalsOOTPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsOOTPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartOOTPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iDeepSuperCluster,&hitsAndEnergies_CaloPartOOTPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simOOTPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsOOTPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsOOTPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[8]);   
             } 
             deepSuperCluster_simOOTPU_nSharedXtals[iSC] = simOOTPU_nSharedXtals[0];   
             deepSuperCluster_simEnergy_sharedXtalsOOTPU[iSC] = simEnergy_sharedXtalsOOTPU[0];   
             deepSuperCluster_recoEnergy_sharedXtalsOOTPU[iSC] = recoEnergy_sharedXtalsOOTPU[0];     
             deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = simEnergy_noHitsFraction_sharedXtalsOOTPU[0];   
             deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsOOTPU[0];        
          }

          if(saveCaloParticles_ && isMC_){ 
             dR_simScore.clear();
             sim_nSharedXtals.clear();
             sim_fraction_noHitsFraction.clear();
             sim_fraction.clear();
             recoToSim_fraction.clear();
             recoToSim_fraction_sharedXtals.clear();  
             simEnergy_sharedXtals.clear(); 
             recoEnergy_sharedXtals.clear(); 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 GlobalPoint caloParticle_position = caloParts_position.at(iCalo);
                 std::vector<double> scores = getScores(&iDeepSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             deepSuperCluster_nXtals.push_back((iDeepSuperCluster.hitsAndFractions()).size());   
             deepSuperCluster_dR_simScore[iSC] = dR_simScore;  
             deepSuperCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             deepSuperCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             deepSuperCluster_sim_fraction[iSC] = sim_fraction;   
             deepSuperCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             deepSuperCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             deepSuperCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             deepSuperCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             deepSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_dR_simScore, 999., false, 0., iSC));
             deepSuperCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_nSharedXtals, -999., true, 0., iSC));
             deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             deepSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction, -999., true, 0., iSC));
             deepSuperCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             deepSuperCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoEnergy_sharedXtals, -999., true, 0., iSC));         
          }

          if(savePFCluster_ && saveDeepSC_){   
             //save clusters and deepSuperClusters mutual info
             reco::CaloCluster caloSeed(*iDeepSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iDeepSuperCluster.clustersBegin(); iBC != iDeepSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC)); 
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) deepSuperCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) deepSuperCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      } 

      // The global deepSuperCluster indexing for EE has an offset = ndeepSuperClusterEB
      iSC = (deepSuperClusterEB.product())->size();
      int iSC_tmp=-1;
      //std::cout << "deepSuperClustersEE size: " << (deepSuperClusterEE.product())->size() << std::endl;
      for(const auto& iDeepSuperCluster : *(deepSuperClusterEE.product())){    

          iSC_tmp++;
        
          deepSuperCluster_seedRawId.push_back(iDeepSuperCluster.seed()->seed().rawId());
          deepSuperCluster_rawEnergy.push_back(reduceFloat(iDeepSuperCluster.rawEnergy(),nBits_));     
          deepSuperCluster_rawESEnergy.push_back(reduceFloat(iDeepSuperCluster.preshowerEnergy(),nBits_));
          deepSuperCluster_energy.push_back(reduceFloat(iDeepSuperCluster.energy(),nBits_));
          deepSuperCluster_eta.push_back(reduceFloat(iDeepSuperCluster.eta(),nBits_));
          deepSuperCluster_phi.push_back(reduceFloat(iDeepSuperCluster.phi(),nBits_));
          deepSuperCluster_etaWidth.push_back(reduceFloat(iDeepSuperCluster.etaWidth(),nBits_));
          deepSuperCluster_phiWidth.push_back(reduceFloat(iDeepSuperCluster.phiWidth(),nBits_));
          deepSuperCluster_R.push_back(reduceFloat(iDeepSuperCluster.position().R(),nBits_)); 
          deepSuperCluster_nPFClusters.push_back(iDeepSuperCluster.clusters().size());  
          math::XYZPoint caloPos = iDeepSuperCluster.seed()->position(); 
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          deepSuperCluster_ieta.push_back(ee_id.ix());
          deepSuperCluster_iphi.push_back(ee_id.iy());
          deepSuperCluster_iz.push_back(ee_id.zside());   

          if(saveShowerShapes_){ 
             reco::CaloCluster caloBC(*iDeepSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), topology);  
             deepSuperCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             deepSuperCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             deepSuperCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     deepSuperCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             deepSuperCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             deepSuperCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             deepSuperCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             deepSuperCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             deepSuperCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             deepSuperCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             deepSuperCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             deepSuperCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             deepSuperCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             deepSuperCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             deepSuperCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             deepSuperCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iDeepSuperCluster.rawEnergy(),nBits_));
             deepSuperCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             deepSuperCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             deepSuperCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             deepSuperCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             deepSuperCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             deepSuperCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             deepSuperCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             deepSuperCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             deepSuperCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             deepSuperCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             deepSuperCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             deepSuperCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             deepSuperCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             deepSuperCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             deepSuperCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             deepSuperCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             deepSuperCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             deepSuperCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             deepSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iDeepSuperCluster.rawEnergy(),nBits_));
             deepSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             deepSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             deepSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_));
          }

          //compute scores  
          if(saveGenParticles_ && isMC_){
             dR_genScore.clear();
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             deepSuperCluster_dR_genScore[iSC] = dR_genScore;        
             deepSuperCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_dR_genScore, 999., false, 0., iSC));
          } 

          if(saveCaloParticlesPU_ && isMC_){ 
             simPU_nSharedXtals.clear();
             simEnergy_sharedXtalsPU.clear(); 
             recoEnergy_sharedXtalsPU.clear();
             simEnergy_noHitsFraction_sharedXtalsPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iDeepSuperCluster,&hitsAndEnergies_CaloPartPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsPU.push_back(scores[8]);   
             } 
             deepSuperCluster_simPU_nSharedXtals[iSC] = simPU_nSharedXtals[0];   
             deepSuperCluster_simEnergy_sharedXtalsPU[iSC] = simEnergy_sharedXtalsPU[0];   
             deepSuperCluster_recoEnergy_sharedXtalsPU[iSC] = recoEnergy_sharedXtalsPU[0];     
             deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU[iSC] = simEnergy_noHitsFraction_sharedXtalsPU[0];   
             deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsPU[0];        
          }

          if(saveCaloParticlesOOTPU_ && isMC_){ 
             simOOTPU_nSharedXtals.clear();
             simEnergy_sharedXtalsOOTPU.clear(); 
             recoEnergy_sharedXtalsOOTPU.clear();
             simEnergy_noHitsFraction_sharedXtalsOOTPU.clear(); 
             recoEnergy_noHitsFraction_sharedXtalsOOTPU.clear();
             for(unsigned int iCalo=0; iCalo<hitsAndEnergies_CaloPartOOTPU.size(); iCalo++){
                 std::vector<double> scores = getScores(&iDeepSuperCluster,&hitsAndEnergies_CaloPartOOTPU.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 simOOTPU_nSharedXtals.push_back(scores[0]);  
                 simEnergy_sharedXtalsOOTPU.push_back(scores[5]);  
                 recoEnergy_sharedXtalsOOTPU.push_back(scores[6]); 
                 simEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[7]);  
                 recoEnergy_noHitsFraction_sharedXtalsOOTPU.push_back(scores[8]);   
             } 
             deepSuperCluster_simOOTPU_nSharedXtals[iSC] = simOOTPU_nSharedXtals[0];   
             deepSuperCluster_simEnergy_sharedXtalsOOTPU[iSC] = simEnergy_sharedXtalsOOTPU[0];   
             deepSuperCluster_recoEnergy_sharedXtalsOOTPU[iSC] = recoEnergy_sharedXtalsOOTPU[0];     
             deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = simEnergy_noHitsFraction_sharedXtalsOOTPU[0];   
             deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU[iSC] = recoEnergy_noHitsFraction_sharedXtalsOOTPU[0];        
          }

          if(saveCaloParticles_ && isMC_){
             dR_simScore.clear();
             sim_nSharedXtals.clear();
             sim_fraction_noHitsFraction.clear();
             sim_fraction.clear();
             recoToSim_fraction.clear();
             recoToSim_fraction_sharedXtals.clear();  
             simEnergy_sharedXtals.clear(); 
             recoEnergy_sharedXtals.clear(); 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 GlobalPoint caloParticle_position = caloParts_position.at(iCalo);
                 std::vector<double> scores = getScores(&iDeepSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             deepSuperCluster_nXtals.push_back((iDeepSuperCluster.hitsAndFractions()).size());   
             deepSuperCluster_dR_simScore[iSC] = dR_simScore;  
             deepSuperCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             deepSuperCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             deepSuperCluster_sim_fraction[iSC] = sim_fraction;   
             deepSuperCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             deepSuperCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             deepSuperCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             deepSuperCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             deepSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_dR_simScore, 999., false, 0., iSC));
             deepSuperCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_nSharedXtals, -999., true, 0., iSC));
             deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             deepSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction, -999., true, 0., iSC));
             deepSuperCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             deepSuperCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoEnergy_sharedXtals, -999., true, 0., iSC)); 
          }

          if(iDeepSuperCluster.preshowerClusters().isAvailable()){
              for(unsigned int iPC=0; iPC<iDeepSuperCluster.preshowerClusters().size(); iPC++){
                  if(!iDeepSuperCluster.preshowerClusters()[iPC].isAvailable()) { continue; } 
                  deepSuperCluster_psCluster_energy[iSC_tmp].push_back(reduceFloat(iDeepSuperCluster.preshowerClusters()[iPC]->energy(),nBits_));
                  deepSuperCluster_psCluster_eta[iSC_tmp].push_back(reduceFloat(iDeepSuperCluster.preshowerClusters()[iPC]->eta(),nBits_));
                  deepSuperCluster_psCluster_phi[iSC_tmp].push_back(reduceFloat(iDeepSuperCluster.preshowerClusters()[iPC]->phi(),nBits_));   
              }
          } 

          if(savePFCluster_ && saveDeepSC_){   
             //save clusters and deepSuperClusters mutual info
             reco::CaloCluster caloSeed(*iDeepSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iDeepSuperCluster.clustersBegin(); iBC != iDeepSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC));  
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) deepSuperCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) deepSuperCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      }

      //save pfCluster_deepSuperClustersIndex
      if(savePFCluster_ && saveDeepSC_){
         for(unsigned int iSC=0; iSC<deepSuperCluster_pfClustersIndex.size(); iSC++)
             for(unsigned int iPF=0; iPF<deepSuperCluster_pfClustersIndex.at(iSC).size(); iPF++)
                 if(deepSuperCluster_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_deepSuperClustersIndex[deepSuperCluster_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
      } 

      //save inverse of matchings
      if(saveGenParticles_ && isMC_){ 
         fillParticleMatchedIndex(&genParticle_deepSuperCluster_dR_genScore_MatchedIndex,&deepSuperCluster_dR_genScore_MatchedIndex);
      } 
      if(saveCaloParticles_ && isMC_){ 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_dR_simScore_MatchedIndex,&deepSuperCluster_dR_simScore_MatchedIndex);
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_nSharedXtals_MatchedIndex,&deepSuperCluster_sim_nSharedXtals_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex,&deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_MatchedIndex,&deepSuperCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_recoToSim_fraction_MatchedIndex,&deepSuperCluster_recoToSim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex,&deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_simEnergy_sharedXtals_MatchedIndex,&deepSuperCluster_simEnergy_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex,&deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex);  
      }
   }

   //Save unClustered pfRechits 
   pfRechit_unClustered.clear();
   if(savePFRechits_ && savePFCluster_ && saveSuperCluster_){
      for(const auto& iPFRechit : *(pfRecHits.product())){

          DetId pf_id(iPFRechit.detId());
          bool pfRecHit_isMatched_ = false;
          
          for(unsigned int iPFCl = 0; iPFCl < hitsAndEnergies_PFCluster.size(); iPFCl++){ 
              for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++){ 
                  if(iPFRechit.detId() == hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId()) pfRecHit_isMatched_ = true;
                  break;   
              }
          }

          if(pf_id.subdetId()==EcalBarrel){
             for(unsigned int iSC = 0; iSC < hitsAndEnergies_SuperClusterEB.size(); iSC++){
                 for(unsigned int i = 0; i < hitsAndEnergies_SuperClusterEB.at(iSC).size(); i++){ 
                     if(iPFRechit.detId() == hitsAndEnergies_SuperClusterEB.at(iSC).at(i).first.rawId()) pfRecHit_isMatched_ = true;
                     break;                   
                 }
             }       
          }else if(pf_id.subdetId()==EcalEndcap){
             for(unsigned int iSC = 0; iSC < hitsAndEnergies_SuperClusterEE.size(); iSC++){
                 for(unsigned int i = 0; i < hitsAndEnergies_SuperClusterEE.at(iSC).size(); i++){ 
                     if(iPFRechit.detId() == hitsAndEnergies_SuperClusterEE.at(iSC).at(i).first.rawId()) pfRecHit_isMatched_ = true;
                     break;                   
                 }
             }        
          }

          if(pfRecHit_isMatched_) continue;

          pfRechit_unClustered.push_back(pf_id); 

          cell = geometry->getPosition(pf_id); 
          pfRecHit_unClustered_energy.push_back(reduceFloat(iPFRechit.energy(),nBits_));    
          pfRecHit_unClustered_eta.push_back(reduceFloat(cell.eta(),nBits_));  
          pfRecHit_unClustered_phi.push_back(reduceFloat(cell.phi(),nBits_)); 
          if(pf_id.subdetId()==EcalBarrel){ 
             EBDetId eb_id(pf_id);  
             pfRecHit_unClustered_ieta.push_back(eb_id.ieta());  
             pfRecHit_unClustered_iphi.push_back(eb_id.iphi());  
             pfRecHit_unClustered_iz.push_back(0);     
          }else if(pf_id.subdetId()==EcalEndcap){
             int iz=-99;
             EEDetId ee_id(pf_id);  
             if(ee_id.zside()<0) iz=-1;
             if(ee_id.zside()>0) iz=1; 
             pfRecHit_unClustered_ieta.push_back(ee_id.ix());  
             pfRecHit_unClustered_iphi.push_back(ee_id.iy());  
             pfRecHit_unClustered_iz.push_back(iz);    
          } 
      }   
   }  
   
   //Save noPF rechits 
   if(saveRechits_){
      for(const auto& iRechit : *(recHitsEB.product())){
          
          DetId rechit_id(iRechit.detid());

          bool rechit_isMatched_;
          for(unsigned int iPFCl = 0; iPFCl < hitsAndEnergies_PFCluster.size(); iPFCl++){ 
              for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++){ 
                  if(rechit_id.rawId() == hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId()) rechit_isMatched_ = true;
                  break;   
              }
          }
          if(rechit_isMatched_) continue;  
           
          std::vector<DetId>::iterator it = std::find(pfRechit_unClustered.begin(), pfRechit_unClustered.end(), rechit_id);   
          if (it != pfRechit_unClustered.end()) continue;  
          
          cell = geometry->getPosition(rechit_id); 
          recHit_noPF_energy.push_back(reduceFloat(iRechit.energy(),nBits_));    
          recHit_noPF_eta.push_back(reduceFloat(cell.eta(),nBits_));  
          recHit_noPF_phi.push_back(reduceFloat(cell.phi(),nBits_)); 

          EBDetId eb_id(rechit_id);  
          recHit_noPF_ieta.push_back(eb_id.ieta());  
          recHit_noPF_iphi.push_back(eb_id.iphi());  
          recHit_noPF_iz.push_back(0);     
      }
      for(const auto& iRechit : *(recHitsEE.product())){

          DetId rechit_id(iRechit.detid());

          bool rechit_isMatched_;
          for(unsigned int iPFCl = 0; iPFCl < hitsAndEnergies_PFCluster.size(); iPFCl++){ 
              for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++){ 
                  if(rechit_id.rawId() == hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId()) rechit_isMatched_ = true;
                  break;   
              }
          }
          if(rechit_isMatched_) continue; 
          
          std::vector<DetId>::iterator it = std::find(pfRechit_unClustered.begin(), pfRechit_unClustered.end(), rechit_id);   
          if (it != pfRechit_unClustered.end()) continue;  
          
          cell = geometry->getPosition(rechit_id); 
          recHit_noPF_energy.push_back(reduceFloat(iRechit.energy(),nBits_));    
          recHit_noPF_eta.push_back(reduceFloat(cell.eta(),nBits_));  
          recHit_noPF_phi.push_back(reduceFloat(cell.phi(),nBits_)); 

          int iz=-99;
          EEDetId ee_id(rechit_id);  
          if(ee_id.zside()<0) iz=-1;
          if(ee_id.zside()>0) iz=1; 
          recHit_noPF_ieta.push_back(ee_id.ix());  
          recHit_noPF_iphi.push_back(ee_id.iy());  
          recHit_noPF_iz.push_back(iz);    
      }   
   }  
   //fill tree for each event
   tree->Fill();
}

void RecoSimDumper::beginJob()
{

}

void RecoSimDumper::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void RecoSimDumper::setTree(TTree* tree)
{
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/I");
   tree->Branch("runId", &runId, "runId/I");
   tree->Branch("rho", &rho, "rho/F"); 
   tree->Branch("nVtx", &nVtx, "nVtx/I");
   if(isMC_){ 
    tree->Branch("truePU", &truePU, "truePU/D");
    tree->Branch("obsPU", &obsPU, "obsPU/D");
   }
   if(saveGenParticles_ && isMC_){
      tree->Branch("genParticle_size", &genParticle_size, "genParticle_size/I");  
      tree->Branch("genParticle_genMotherIndex","std::vector<int>",&genParticle_genMotherIndex); 
      tree->Branch("genParticle_pdgId","std::vector<int>",&genParticle_pdgId);
      tree->Branch("genParticle_status","std::vector<int>",&genParticle_status); 
      tree->Branch("genParticle_statusFlag","std::vector<int>",&genParticle_statusFlag); 
      tree->Branch("genParticle_energy","std::vector<float>",&genParticle_energy);
      tree->Branch("genParticle_pt","std::vector<float>",&genParticle_pt);
      tree->Branch("genParticle_eta","std::vector<float>",&genParticle_eta);
      tree->Branch("genParticle_phi","std::vector<float>",&genParticle_phi);
      if(savePFCluster_) tree->Branch("genParticle_pfCluster_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_pfCluster_dR_genScore_MatchedIndex);
      if(saveSuperCluster_) tree->Branch("genParticle_superCluster_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_superCluster_dR_genScore_MatchedIndex);
      if(saveRetunedSC_) tree->Branch("genParticle_retunedSuperCluster_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_retunedSuperCluster_dR_genScore_MatchedIndex); 
      if(saveDeepSC_) tree->Branch("genParticle_deepSuperCluster_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_deepSuperCluster_dR_genScore_MatchedIndex); 
   }
   if(saveCaloParticlesPU_ && isMC_){
      tree->Branch("caloParticlePU_nHitsWithES", &caloParticlePU_nHitsWithES, "caloParticlePU_nHitsWithES/I");  
      tree->Branch("caloParticlePU_nHits", &caloParticlePU_nHits, "caloParticlePU_nHits/I");  
      tree->Branch("caloParticlePU_totEnergyWithES", &caloParticlePU_totEnergyWithES, "caloParticlePU_totEnergyWithES/F");  
      tree->Branch("caloParticlePU_totEnergy", &caloParticlePU_totEnergy, "caloParticlePU_totEnergy/F");  
      if(saveSimhitsPU_){ 
         tree->Branch("caloParticlePU_xtalEnergy", "std::vector<float>", &caloParticlePU_xtalEnergy);  
         tree->Branch("caloParticlePU_xtalEta", "std::vector<float>", &caloParticlePU_xtalEta);   
         tree->Branch("caloParticlePU_xtalPhi", "std::vector<float>", &caloParticlePU_xtalPhi);                
         tree->Branch("caloParticlePU_xtalIeta", "std::vector<int>", &caloParticlePU_xtalIeta);   
         tree->Branch("caloParticlePU_xtalIphi", "std::vector<int>", &caloParticlePU_xtalIphi);  
         tree->Branch("caloParticlePU_xtalIz", "std::vector<int>", &caloParticlePU_xtalIz);              
         tree->Branch("caloParticlePU_xtalIplane", "std::vector<int>", &caloParticlePU_xtalIplane);   
      }  
   }
   if(saveCaloParticlesOOTPU_ && isMC_){
      tree->Branch("caloParticleOOTPU_nHitsWithES", &caloParticlePU_nHitsWithES, "caloParticlePU_nHitsWithES/I");  
      tree->Branch("caloParticleOOTPU_nHits", &caloParticlePU_nHits, "caloParticlePU_nHits/I");  
      tree->Branch("caloParticleOOTPU_totEnergyWithES", &caloParticlePU_totEnergyWithES, "caloParticlePU_totEnergyWithES/F");  
      tree->Branch("caloParticleOOTPU_totEnergy", &caloParticlePU_totEnergy, "caloParticlePU_totEnergy/F");  
      if(saveSimhitsPU_){ 
         tree->Branch("caloParticleOOTPU_xtalEnergy", "std::vector<float>", &caloParticleOOTPU_xtalEnergy);  
         tree->Branch("caloParticleOOTPU_xtalEta", "std::vector<float>", &caloParticleOOTPU_xtalEta);   
         tree->Branch("caloParticleOOTPU_xtalPhi", "std::vector<float>", &caloParticleOOTPU_xtalPhi);                
         tree->Branch("caloParticleOOTPU_xtalIeta", "std::vector<int>", &caloParticleOOTPU_xtalIeta);   
         tree->Branch("caloParticleOOTPU_xtalIphi", "std::vector<int>", &caloParticleOOTPU_xtalIphi);  
         tree->Branch("caloParticleOOTPU_xtalIz", "std::vector<int>", &caloParticleOOTPU_xtalIz);              
         tree->Branch("caloParticleOOTPU_xtalIplane", "std::vector<int>", &caloParticleOOTPU_xtalIplane);   
      }   
   }
   if(saveCaloParticles_ && isMC_){
      tree->Branch("caloParticle_size", &caloParticle_size, "caloParticle_size/I"); 
      tree->Branch("caloParticle_index","std::vector<int>",&caloParticle_index); 
      tree->Branch("caloParticle_nXtals","std::vector<int>",&caloParticle_nXtals);   
      tree->Branch("caloParticle_pdgId","std::vector<int>",&caloParticle_pdgId); 
      tree->Branch("caloParticle_status","std::vector<int>",&caloParticle_status); 
      tree->Branch("caloParticle_charge","std::vector<int>",&caloParticle_charge); 
      tree->Branch("caloParticle_genEnergy","std::vector<float>",&caloParticle_genEnergy);
      tree->Branch("caloParticle_simEnergyWithES","std::vector<float>",&caloParticle_simEnergyWithES); 
      tree->Branch("caloParticle_simEnergy","std::vector<float>",&caloParticle_simEnergy);        
      tree->Branch("caloParticle_simEnergyGoodStatus","std::vector<float>",&caloParticle_simEnergyGoodStatus);  
      tree->Branch("caloParticle_genPt","std::vector<float>",&caloParticle_genPt);
      tree->Branch("caloParticle_simPt","std::vector<float>",&caloParticle_simPt);
      tree->Branch("caloParticle_genEta","std::vector<float>",&caloParticle_genEta);
      tree->Branch("caloParticle_simEta","std::vector<float>",&caloParticle_simEta);
      tree->Branch("caloParticle_genPhi","std::vector<float>",&caloParticle_genPhi);
      tree->Branch("caloParticle_partonIndex","std::vector<int>",&caloParticle_partonIndex);
      tree->Branch("caloParticle_partonPdgId","std::vector<int>",&caloParticle_partonPdgId);
      tree->Branch("caloParticle_partonCharge","std::vector<int>",&caloParticle_partonCharge);
      tree->Branch("caloParticle_partonEnergy","std::vector<float>",&caloParticle_partonEnergy);
      tree->Branch("caloParticle_partonPt","std::vector<float>",&caloParticle_partonPt);
      tree->Branch("caloParticle_partonEta","std::vector<float>",&caloParticle_partonEta);
      tree->Branch("caloParticle_partonPhi","std::vector<float>",&caloParticle_partonPhi);
      tree->Branch("caloParticle_genMotherPdgId","std::vector<int>",&caloParticle_genMotherPdgId);
      tree->Branch("caloParticle_genMotherStatus","std::vector<int>",&caloParticle_genMotherStatus);
      tree->Branch("caloParticle_genMotherCharge","std::vector<int>",&caloParticle_genMotherCharge);
      tree->Branch("caloParticle_genMotherEnergy","std::vector<float>",&caloParticle_genMotherEnergy);
      tree->Branch("caloParticle_genMotherPt","std::vector<float>",&caloParticle_genMotherPt);
      tree->Branch("caloParticle_genMotherEta","std::vector<float>",&caloParticle_genMotherEta);
      tree->Branch("caloParticle_genMotherPhi","std::vector<float>",&caloParticle_genMotherPhi);
      tree->Branch("caloParticle_simPhi","std::vector<float>",&caloParticle_simPhi);
      tree->Branch("caloParticle_simIeta","std::vector<int>",&caloParticle_simIeta);
      tree->Branch("caloParticle_simIphi","std::vector<int>",&caloParticle_simIphi);
      tree->Branch("caloParticle_simIz","std::vector<int>",&caloParticle_simIz);
      tree->Branch("caloParticle_nSharedXtals", "std::vector<int>", &caloParticle_nSharedXtals);
      tree->Branch("caloParticle_sharedIndex1", "std::vector<int>", &caloParticle_sharedIndex1);
      tree->Branch("caloParticle_sharedIndex2", "std::vector<int>", &caloParticle_sharedIndex2);
      tree->Branch("caloParticle_sharedEnergyFrac1", "std::vector<float>", &caloParticle_sharedEnergyFrac1);
      tree->Branch("caloParticle_sharedEnergyFrac2", "std::vector<float>", &caloParticle_sharedEnergyFrac2);
      if(savePFCluster_){   
         tree->Branch("caloParticle_pfCluster_dR_simScore_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_dR_simScore_MatchedIndex);
         tree->Branch("caloParticle_pfCluster_sim_nSharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_sim_nSharedXtals_MatchedIndex);
         tree->Branch("caloParticle_pfCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_sim_fraction_noHitsFraction_MatchedIndex); 
         tree->Branch("caloParticle_pfCluster_sim_fraction_MatchedIndex", "std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_MatchedIndex);      
         tree->Branch("caloParticle_pfCluster_recoToSim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_recoToSim_fraction_MatchedIndex);   
         tree->Branch("caloParticle_pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex);   
         tree->Branch("caloParticle_pfCluster_simEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_simEnergy_sharedXtals_MatchedIndex);  
         tree->Branch("caloParticle_pfCluster_recoEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_recoEnergy_sharedXtals_MatchedIndex);         
      }
      if(saveSuperCluster_){    
         tree->Branch("caloParticle_superCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >", &caloParticle_superCluster_dR_simScore_MatchedIndex);
         tree->Branch("caloParticle_superCluster_sim_nSharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_sim_nSharedXtals_MatchedIndex);
         tree->Branch("caloParticle_superCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_sim_fraction_noHitsFraction_MatchedIndex); 
         tree->Branch("caloParticle_superCluster_sim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_sim_fraction_MatchedIndex);      
         tree->Branch("caloParticle_superCluster_recoToSim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_recoToSim_fraction_MatchedIndex);   
         tree->Branch("caloParticle_superCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_recoToSim_fraction_sharedXtals_MatchedIndex);   
         tree->Branch("caloParticle_superCluster_simEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_simEnergy_sharedXtals_MatchedIndex);  
         tree->Branch("caloParticle_superCluster_recoEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_recoEnergy_sharedXtals_MatchedIndex);
         if(saveRetunedSC_){
            tree->Branch("caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex);
            tree->Branch("caloParticle_retunedSuperCluster_sim_nSharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_sim_nSharedXtals_MatchedIndex);
            tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex); 
            tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex);      
            tree->Branch("caloParticle_retunedSuperCluster_recoToSim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_recoToSim_fraction_MatchedIndex);   
            tree->Branch("caloParticle_retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);   
            tree->Branch("caloParticle_retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex);  
            tree->Branch("caloParticle_retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex);
         } 
         if(saveDeepSC_){
            tree->Branch("caloParticle_deepSuperCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_dR_simScore_MatchedIndex);
            tree->Branch("caloParticle_deepSuperCluster_sim_nSharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_sim_nSharedXtals_MatchedIndex);
            tree->Branch("caloParticle_deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex); 
            tree->Branch("caloParticle_deepSuperCluster_sim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_sim_fraction_MatchedIndex);      
            tree->Branch("caloParticle_deepSuperCluster_recoToSim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_recoToSim_fraction_MatchedIndex);   
            tree->Branch("caloParticle_deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);   
            tree->Branch("caloParticle_deepSuperCluster_simEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_simEnergy_sharedXtals_MatchedIndex);  
            tree->Branch("caloParticle_deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex); 
         } 
      }
      if(saveSimhits_){
         tree->Branch("simHit_energy","std::vector<std::vector<float> >",&simHit_energy);
         tree->Branch("simHit_eta","std::vector<std::vector<float> >",&simHit_eta);
         tree->Branch("simHit_phi","std::vector<std::vector<float> >",&simHit_phi);
         tree->Branch("simHit_ieta","std::vector<std::vector<int> >",&simHit_ieta);
         tree->Branch("simHit_iphi","std::vector<std::vector<int> >",&simHit_iphi);
         tree->Branch("simHit_iz","std::vector<std::vector<int> >",&simHit_iz);
         tree->Branch("simHit_iplane","std::vector<std::vector<int> >",&simHit_iplane);
         tree->Branch("simHit_chStatus","std::vector<std::vector<int> >",&simHit_chStatus);
      }
   }
   if(saveRechits_){
      tree->Branch("recHit_noPF_energy","std::vector<float>",&recHit_noPF_energy);
      tree->Branch("recHit_noPF_eta","std::vector<float>",&recHit_noPF_eta); 
      tree->Branch("recHit_noPF_phi","std::vector<float>",&recHit_noPF_phi);
      tree->Branch("recHit_noPF_ieta","std::vector<int>",&recHit_noPF_ieta); 
      tree->Branch("recHit_noPF_iphi","std::vector<int>",&recHit_noPF_iphi);
      tree->Branch("recHit_noPF_iz","std::vector<int>",&recHit_noPF_iz);     
   }
   if(savePFRechits_ && savePFCluster_ && saveSuperCluster_){ 
      tree->Branch("pfRecHit_unClustered_energy","std::vector<float>",&pfRecHit_unClustered_energy);
      tree->Branch("pfRecHit_unClustered_eta","std::vector<float>",&pfRecHit_unClustered_eta); 
      tree->Branch("pfRecHit_unClustered_phi","std::vector<float>",&pfRecHit_unClustered_phi);
      tree->Branch("pfRecHit_unClustered_ieta","std::vector<int>",&pfRecHit_unClustered_ieta); 
      tree->Branch("pfRecHit_unClustered_iphi","std::vector<int>",&pfRecHit_unClustered_iphi);
      tree->Branch("pfRecHit_unClustered_iz","std::vector<int>",&pfRecHit_unClustered_iz);     
   }
   if(savePFCluster_){
      tree->Branch("pfCluster_rawEnergy","std::vector<float>",&pfCluster_rawEnergy);
      tree->Branch("pfCluster_energy","std::vector<float>",&pfCluster_energy);
      tree->Branch("pfCluster_rawPt","std::vector<float>",&pfCluster_rawPt); 
      tree->Branch("pfCluster_pt","std::vector<float>",&pfCluster_pt);  
      tree->Branch("pfCluster_eta","std::vector<float>",&pfCluster_eta);
      tree->Branch("pfCluster_phi","std::vector<float>",&pfCluster_phi); 
      tree->Branch("pfCluster_ieta","std::vector<int>",&pfCluster_ieta);
      tree->Branch("pfCluster_iphi","std::vector<int>",&pfCluster_iphi);   
      tree->Branch("pfCluster_iz","std::vector<int>",&pfCluster_iz);
      tree->Branch("pfCluster_nXtals","std::vector<double>",&pfCluster_nXtals);  
      if(saveShowerShapes_){  
         tree->Branch("pfCluster_etaWidth","std::vector<float>",&pfCluster_etaWidth); 
         tree->Branch("pfCluster_phiWidth","std::vector<float>",&pfCluster_phiWidth); 
         tree->Branch("pfCluster_e5x5","std::vector<float>",&pfCluster_e5x5);
         tree->Branch("pfCluster_e2x2Ratio","std::vector<float>",&pfCluster_e2x2Ratio);
         tree->Branch("pfCluster_e3x3Ratio","std::vector<float>",&pfCluster_e3x3Ratio);
         tree->Branch("pfCluster_eMaxRatio","std::vector<float>",&pfCluster_eMaxRatio);
         tree->Branch("pfCluster_e2ndRatio","std::vector<float>",&pfCluster_e2ndRatio);
         tree->Branch("pfCluster_eTopRatio","std::vector<float>",&pfCluster_eTopRatio);
         tree->Branch("pfCluster_eRightRatio","std::vector<float>",&pfCluster_eRightRatio);
         tree->Branch("pfCluster_eBottomRatio","std::vector<float>",&pfCluster_eBottomRatio);
         tree->Branch("pfCluster_eLeftRatio","std::vector<float>",&pfCluster_eLeftRatio);
         tree->Branch("pfCluster_e2x5MaxRatio","std::vector<float>",&pfCluster_e2x5MaxRatio);
         tree->Branch("pfCluster_e2x5TopRatio","std::vector<float>",&pfCluster_e2x5TopRatio);
         tree->Branch("pfCluster_e2x5RightRatio","std::vector<float>",&pfCluster_e2x5RightRatio);
         tree->Branch("pfCluster_e2x5BottomRatio","std::vector<float>",&pfCluster_e2x5BottomRatio); 
         tree->Branch("pfCluster_e2x5LeftRatio","std::vector<float>",&pfCluster_e2x5LeftRatio); 
         tree->Branch("pfCluster_swissCross","std::vector<float>",&pfCluster_swissCross); 
         tree->Branch("pfCluster_r9","std::vector<float>",&pfCluster_r9);
         tree->Branch("pfCluster_sigmaIetaIeta","std::vector<float>",&pfCluster_sigmaIetaIeta);
         tree->Branch("pfCluster_sigmaIetaIphi","std::vector<float>",&pfCluster_sigmaIetaIphi);
         tree->Branch("pfCluster_sigmaIphiIphi","std::vector<float>",&pfCluster_sigmaIphiIphi);
         tree->Branch("pfCluster_full5x5_e5x5","std::vector<float>",&pfCluster_full5x5_e5x5);
         tree->Branch("pfCluster_full5x5_e2x2Ratio","std::vector<float>",&pfCluster_full5x5_e2x2Ratio);
         tree->Branch("pfCluster_full5x5_e3x3Ratio","std::vector<float>",&pfCluster_full5x5_e3x3Ratio);
         tree->Branch("pfCluster_full5x5_eMaxRatio","std::vector<float>",&pfCluster_full5x5_eMaxRatio);
         tree->Branch("pfCluster_full5x5_e2ndRatio","std::vector<float>",&pfCluster_full5x5_e2ndRatio);
         tree->Branch("pfCluster_full5x5_eTopRatio","std::vector<float>",&pfCluster_full5x5_eTopRatio);
         tree->Branch("pfCluster_full5x5_eRightRatio","std::vector<float>",&pfCluster_full5x5_eRightRatio);
         tree->Branch("pfCluster_full5x5_eBottomRatio","std::vector<float>",&pfCluster_full5x5_eBottomRatio);
         tree->Branch("pfCluster_full5x5_eLeftRatio","std::vector<float>",&pfCluster_full5x5_eLeftRatio);
         tree->Branch("pfCluster_full5x5_e2x5MaxRatio","std::vector<float>",&pfCluster_full5x5_e2x5MaxRatio);
         tree->Branch("pfCluster_full5x5_e2x5TopRatio","std::vector<float>",&pfCluster_full5x5_e2x5TopRatio);
         tree->Branch("pfCluster_full5x5_e2x5RightRatio","std::vector<float>",&pfCluster_full5x5_e2x5RightRatio);
         tree->Branch("pfCluster_full5x5_e2x5BottomRatio","std::vector<float>",&pfCluster_full5x5_e2x5BottomRatio); 
         tree->Branch("pfCluster_full5x5_e2x5LeftRatio","std::vector<float>",&pfCluster_full5x5_e2x5LeftRatio); 
         tree->Branch("pfCluster_full5x5_swissCross","std::vector<float>",&pfCluster_full5x5_swissCross); 
         tree->Branch("pfCluster_full5x5_r9","std::vector<float>",&pfCluster_full5x5_r9);
         tree->Branch("pfCluster_full5x5_sigmaIetaIeta","std::vector<float>",&pfCluster_full5x5_sigmaIetaIeta);
         tree->Branch("pfCluster_full5x5_sigmaIetaIphi","std::vector<float>",&pfCluster_full5x5_sigmaIetaIphi);
         tree->Branch("pfCluster_full5x5_sigmaIphiIphi","std::vector<float>",&pfCluster_full5x5_sigmaIphiIphi);
      }      
      if(isMC_ && saveCaloParticles_ && saveCaloParticlesPU_){
         tree->Branch("pfCluster_rawEnergyUncalib","std::vector<float>",&pfCluster_rawEnergyUncalib);
         tree->Branch("pfCluster_noise","std::vector<float>",&pfCluster_noise);
         tree->Branch("pfCluster_noiseUncalib","std::vector<float>",&pfCluster_noiseUncalib);
         tree->Branch("pfCluster_noiseNoFractions","std::vector<float>",&pfCluster_noiseNoFractions);
         tree->Branch("pfCluster_noiseUncalibNoFractions","std::vector<float>",&pfCluster_noiseUncalibNoFractions);
         tree->Branch("pfCluster_noiseDB","std::vector<float>",&pfCluster_noiseDB);
         tree->Branch("pfCluster_noiseDBUncalib","std::vector<float>",&pfCluster_noiseDBUncalib);
         tree->Branch("pfCluster_noiseDBNoFractions","std::vector<float>",&pfCluster_noiseDBNoFractions);
         tree->Branch("pfCluster_noiseDBUncalibNoFractions","std::vector<float>",&pfCluster_noiseDBUncalibNoFractions);
      }
      if(saveSuperCluster_) tree->Branch("pfCluster_superClustersIndex","std::vector<std::vector<int> >",&pfCluster_superClustersIndex); 
      if(saveRetunedSC_) tree->Branch("pfCluster_retunedSuperClustersIndex","std::vector<std::vector<int> >",&pfCluster_retunedSuperClustersIndex);  
      if(saveDeepSC_) tree->Branch("pfCluster_deepSuperClustersIndex","std::vector<std::vector<int> >",&pfCluster_deepSuperClustersIndex); 
      if(saveCaloParticles_ && isMC_){ 
         tree->Branch("pfCluster_dR_genScore","std::vector<std::vector<double> >",&pfCluster_dR_genScore);
         tree->Branch("pfCluster_dR_simScore","std::vector<std::vector<double> >",&pfCluster_dR_simScore);
         tree->Branch("pfCluster_sim_nSharedXtals","std::vector<std::vector<double> >",&pfCluster_sim_nSharedXtals);
         tree->Branch("pfCluster_sim_fraction_noHitsFraction","std::vector<std::vector<double> >",&pfCluster_sim_fraction_noHitsFraction); 
         tree->Branch("pfCluster_sim_fraction","std::vector<std::vector<double> >",&pfCluster_sim_fraction);
         tree->Branch("pfCluster_recoToSim_fraction","std::vector<std::vector<double> >",&pfCluster_recoToSim_fraction);
         tree->Branch("pfCluster_recoToSim_fraction_sharedXtals","std::vector<std::vector<double> >",&pfCluster_recoToSim_fraction_sharedXtals);
         tree->Branch("pfCluster_simEnergy_sharedXtals","std::vector<std::vector<double> >",&pfCluster_simEnergy_sharedXtals);
         tree->Branch("pfCluster_recoEnergy_sharedXtals","std::vector<std::vector<double> >",&pfCluster_recoEnergy_sharedXtals);  
         tree->Branch("pfCluster_dR_genScore_MatchedIndex","std::vector<int>",&pfCluster_dR_genScore_MatchedIndex);
         tree->Branch("pfCluster_dR_simScore_MatchedIndex","std::vector<int>",&pfCluster_dR_simScore_MatchedIndex);
         tree->Branch("pfCluster_sim_nSharedXtals_MatchedIndex","std::vector<int>",&pfCluster_sim_nSharedXtals_MatchedIndex);
         tree->Branch("pfCluster_sim_fraction_noHitsFraction_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_noHitsFraction_MatchedIndex); 
         tree->Branch("pfCluster_sim_fraction_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_MatchedIndex);
         tree->Branch("pfCluster_recoToSim_fraction_MatchedIndex","std::vector<int>",&pfCluster_recoToSim_fraction_MatchedIndex);
         tree->Branch("pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex","std::vector<int>",&pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex);
         tree->Branch("pfCluster_simEnergy_sharedXtals_MatchedIndex","std::vector<int>",&pfCluster_simEnergy_sharedXtals_MatchedIndex);
         tree->Branch("pfCluster_recoEnergy_sharedXtals_MatchedIndex","std::vector<int>",&pfCluster_recoEnergy_sharedXtals_MatchedIndex);
      }
      if(saveCaloParticlesPU_ && isMC_){ 
         tree->Branch("pfCluster_simPU_nSharedXtals","std::vector<double>",&pfCluster_simPU_nSharedXtals);
         tree->Branch("pfCluster_simEnergy_sharedXtalsPU","std::vector<double>",&pfCluster_simEnergy_sharedXtalsPU);
         tree->Branch("pfCluster_recoEnergy_sharedXtalsPU","std::vector<double>",&pfCluster_recoEnergy_sharedXtalsPU); 
         tree->Branch("pfCluster_simEnergy_noHitsFraction_sharedXtalsPU","std::vector<double>",&pfCluster_simEnergy_noHitsFraction_sharedXtalsPU);
         tree->Branch("pfCluster_recoEnergy_noHitsFraction_sharedXtalsPU","std::vector<double>",&pfCluster_recoEnergy_noHitsFraction_sharedXtalsPU);  
      }
      if(saveCaloParticlesOOTPU_ && isMC_){ 
         tree->Branch("pfCluster_simOOTPU_nSharedXtals","std::vector<double>",&pfCluster_simOOTPU_nSharedXtals);
         tree->Branch("pfCluster_simEnergy_sharedXtalsOOTPU","std::vector<double>",&pfCluster_simEnergy_sharedXtalsOOTPU);
         tree->Branch("pfCluster_recoEnergy_sharedXtalsOOTPU","std::vector<double>",&pfCluster_recoEnergy_sharedXtalsOOTPU);  
         tree->Branch("pfCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU","std::vector<double>",&pfCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU);
         tree->Branch("pfCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU","std::vector<double>",&pfCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU);  
      }
      if(savePFClusterhits_ && savePFCluster_){ 
         tree->Branch("pfClusterHit_fraction","std::vector<std::vector<float> >",&pfClusterHit_fraction);
         tree->Branch("pfClusterHit_rechitEnergy","std::vector<std::vector<float> >",&pfClusterHit_rechitEnergy);
         tree->Branch("pfClusterHit_eta","std::vector<std::vector<float> >",&pfClusterHit_eta);
         tree->Branch("pfClusterHit_phi","std::vector<std::vector<float> >",&pfClusterHit_phi);
         tree->Branch("pfClusterHit_ieta","std::vector<std::vector<int> >",&pfClusterHit_ieta);
         tree->Branch("pfClusterHit_iphi","std::vector<std::vector<int> >",&pfClusterHit_iphi);
         tree->Branch("pfClusterHit_iz","std::vector<std::vector<int> >",&pfClusterHit_iz);
         tree->Branch("pfClusterHit_adcToGeV","std::vector<std::vector<float> >",&pfClusterHit_adcToGeV); 
         tree->Branch("pfClusterHit_laserCorr","std::vector<std::vector<float> >",&pfClusterHit_laserCorr);    
         tree->Branch("pfClusterHit_ic","std::vector<std::vector<float> >",&pfClusterHit_ic);        
         tree->Branch("pfClusterHit_icMC","std::vector<std::vector<float> >",&pfClusterHit_icMC);
         tree->Branch("pfClusterHit_chStatus","std::vector<std::vector<int> >",&pfClusterHit_chStatus);
      }   
   }
   
   if(saveGsfElectrons_){
      tree->Branch("gsfElectron_index","std::vector<int> ",&gsfElectron_index);    
      tree->Branch("gsfElectron_seedRawId","std::vector<uint32_t> ",&gsfElectron_seedRawId);   
      tree->Branch("gsfElectron_isEB","std::vector<bool> ",&gsfElectron_isEB);   
      tree->Branch("gsfElectron_isEE","std::vector<bool> ",&gsfElectron_isEE);
      tree->Branch("gsfElectron_isEBEEGap","std::vector<bool> ",&gsfElectron_isEBEEGap);
      tree->Branch("gsfElectron_isEBEtaGap","std::vector<bool> ",&gsfElectron_isEBEtaGap);
      tree->Branch("gsfElectron_isEBPhiGap","std::vector<bool> ",&gsfElectron_isEBPhiGap);
      tree->Branch("gsfElectron_isEEDeeGap","std::vector<bool> ",&gsfElectron_isEEDeeGap);
      tree->Branch("gsfElectron_isEERingGap","std::vector<bool> ",&gsfElectron_isEERingGap);
      tree->Branch("gsfElectron_isEcalDriven","std::vector<bool> ",&gsfElectron_isEcalDriven);
      tree->Branch("gsfElectron_isTrackerDriven","std::vector<bool> ",&gsfElectron_isTrackerDriven);
      tree->Branch("gsfElectron_classification","std::vector<int> ",&gsfElectron_classification);  
      tree->Branch("gsfElectron_scNPFClusters","std::vector<int> ",&gsfElectron_scNPFClusters); 
      tree->Branch("gsfElectron_ecalSCNPFClusters","std::vector<int> ",&gsfElectron_ecalSCNPFClusters);     
      tree->Branch("gsfElectron_p","std::vector<float> ",&gsfElectron_p);     
      tree->Branch("gsfElectron_pt","std::vector<float> ",&gsfElectron_pt);     
      tree->Branch("gsfElectron_et","std::vector<float> ",&gsfElectron_et);     
      tree->Branch("gsfElectron_energy","std::vector<float> ",&gsfElectron_energy);   
      tree->Branch("gsfElectron_energyErr","std::vector<float> ",&gsfElectron_energyErr);   
      tree->Branch("gsfElectron_ecalEnergy","std::vector<float> ",&gsfElectron_ecalEnergy); 
      tree->Branch("gsfElectron_ecalEnergyErr","std::vector<float> ",&gsfElectron_ecalEnergyErr);  
      tree->Branch("gsfElectron_eta","std::vector<float> ",&gsfElectron_eta); 
      tree->Branch("gsfElectron_phi","std::vector<float> ",&gsfElectron_phi); 
      tree->Branch("gsfElectron_trkEtaMode","std::vector<float> ",&gsfElectron_trkEtaMode);   
      tree->Branch("gsfElectron_trkPhiMode","std::vector<float> ",&gsfElectron_trkPhiMode); 
      tree->Branch("gsfElectron_trkPMode","std::vector<float> ",&gsfElectron_trkPMode);  
      tree->Branch("gsfElectron_trkPModeErr","std::vector<float> ",&gsfElectron_trkPModeErr); 
      tree->Branch("gsfElectron_trkPInn","std::vector<float> ",&gsfElectron_trkPInn); 
      tree->Branch("gsfElectron_trkPtInn","std::vector<float> ",&gsfElectron_trkPtInn);   
      tree->Branch("gsfElectron_trkPVtx","std::vector<float> ",&gsfElectron_trkPVtx); 
      tree->Branch("gsfElectron_trkPOut","std::vector<float> ",&gsfElectron_trkPOut);  
      tree->Branch("gsfElectron_trkChi2","std::vector<float> ",&gsfElectron_trkChi2); 
      tree->Branch("gsfElectron_trkNDof","std::vector<float> ",&gsfElectron_trkNDof);  
      tree->Branch("gsfElectron_trkFbrem","std::vector<float> ",&gsfElectron_trackFbrem);  
      tree->Branch("gsfElectron_superClusterFbrem","std::vector<float> ",&gsfElectron_superClusterFbrem);   
      tree->Branch("gsfElectron_hademTow","std::vector<float> ",&gsfElectron_hademTow); 
      tree->Branch("gsfElectron_hademCone","std::vector<float> ",&gsfElectron_hademCone);  
      tree->Branch("gsfElectron_ecalDrivenSeed","std::vector<float> ",&gsfElectron_ecalDrivenSeed); 
      tree->Branch("gsfElectron_nrSatCrys","std::vector<float> ",&gsfElectron_nrSatCrys); 
      tree->Branch("gsfElectron_scEta","std::vector<float> ",&gsfElectron_scEta);  
      tree->Branch("gsfElectron_scPhi","std::vector<float> ",&gsfElectron_scPhi);  
      tree->Branch("gsfElectron_scEnergy","std::vector<float> ",&gsfElectron_scEnergy);   
      tree->Branch("gsfElectron_scRawEnergy","std::vector<float> ",&gsfElectron_scRawEnergy);   
      tree->Branch("gsfElectron_scRawESEnergy","std::vector<float> ",&gsfElectron_scRawESEnergy); 
      tree->Branch("gsfElectron_scEt","std::vector<float> ",&gsfElectron_scEt); 
      tree->Branch("gsfElectron_scEtaWidth","std::vector<float> ",&gsfElectron_scEtaWidth);   
      tree->Branch("gsfElectron_scPhiWidth","std::vector<float> ",&gsfElectron_scPhiWidth);   
      tree->Branch("gsfElectron_scEoP","std::vector<float> ",&gsfElectron_scEoP);   
      tree->Branch("gsfElectron_ecalSCEta","std::vector<float> ",&gsfElectron_ecalSCEta);  
      tree->Branch("gsfElectron_ecalSCPhi","std::vector<float> ",&gsfElectron_ecalSCPhi);   
      tree->Branch("gsfElectron_ecalSCEnergy","std::vector<float> ",&gsfElectron_ecalSCEnergy);   
      tree->Branch("gsfElectron_ecalSCRawEnergy","std::vector<float> ",&gsfElectron_ecalSCRawEnergy);   
      tree->Branch("gsfElectron_ecalSCRawESEnergy","std::vector<float> ",&gsfElectron_ecalSCRawESEnergy); 
      tree->Branch("gsfElectron_ecalSCEt","std::vector<float> ",&gsfElectron_ecalSCEt); 
      tree->Branch("gsfElectron_ecalSCEtaWidth","std::vector<float> ",&gsfElectron_ecalSCEtaWidth);   
      tree->Branch("gsfElectron_ecalSCPhiWidth","std::vector<float> ",&gsfElectron_ecalSCPhiWidth);   
      tree->Branch("gsfElectron_ecalSCEoP","std::vector<float> ",&gsfElectron_ecalSCEoP);       
      tree->Branch("gsfElectron_scSwissCross","std::vector<float> ",&gsfElectron_scSwissCross);  
      tree->Branch("gsfElectron_scEMax","std::vector<float> ",&gsfElectron_scEMax);  
      tree->Branch("gsfElectron_scE2x2","std::vector<float> ",&gsfElectron_scE2x2);  
      tree->Branch("gsfElectron_scE3x3","std::vector<float> ",&gsfElectron_scE3x3);  
      tree->Branch("gsfElectron_scE5x5","std::vector<float> ",&gsfElectron_scE5x5); 
      tree->Branch("gsfElectron_scR9","std::vector<float> ",&gsfElectron_scR9);  
      tree->Branch("gsfElectron_scSigmaIEtaIEta","std::vector<float> ",&gsfElectron_scSigmaIEtaIEta); 
      tree->Branch("gsfElectron_scSigmaIEtaIPhi","std::vector<float> ",&gsfElectron_scSigmaIEtaIPhi);  
      tree->Branch("gsfElectron_scSigmaIPhiIPhi","std::vector<float> ",&gsfElectron_scSigmaIPhiIPhi); 
      tree->Branch("gsfElectron_full5x5_scEMax","std::vector<float> ",&gsfElectron_full5x5_scEMax);  
      tree->Branch("gsfElectron_full5x5_scE2x2","std::vector<float> ",&gsfElectron_full5x5_scE2x2);  
      tree->Branch("gsfElectron_full5x5_scE3x3","std::vector<float> ",&gsfElectron_full5x5_scE3x3);  
      tree->Branch("gsfElectron_full5x5_scE5x5","std::vector<float> ",&gsfElectron_full5x5_scE5x5); 
      tree->Branch("gsfElectron_full5x5_scR9","std::vector<float> ",&gsfElectron_full5x5_scR9);  
      tree->Branch("gsfElectron_full5x5_scSigmaIEtaIEta","std::vector<float> ",&gsfElectron_full5x5_scSigmaIEtaIEta); 
      tree->Branch("gsfElectron_full5x5_scSigmaIEtaIPhi","std::vector<float> ",&gsfElectron_full5x5_scSigmaIEtaIPhi);  
      tree->Branch("gsfElectron_full5x5_scSigmaIPhiIPhi","std::vector<float> ",&gsfElectron_full5x5_scSigmaIPhiIPhi);  
      tree->Branch("gsfElectron_HoE","std::vector<float> ",&gsfElectron_HoE); 
      tree->Branch("gsfElectron_trkIso03","std::vector<float> ",&gsfElectron_trkIso03); 
      tree->Branch("gsfElectron_ecalIso03","std::vector<float> ",&gsfElectron_ecalIso03); 
      tree->Branch("gsfElectron_hcalIso03","std::vector<float> ",&gsfElectron_hcalIso03);  
      tree->Branch("gsfElectron_trkIso04","std::vector<float> ",&gsfElectron_trkIso04); 
      tree->Branch("gsfElectron_ecalIso04","std::vector<float> ",&gsfElectron_ecalIso04); 
      tree->Branch("gsfElectron_hcalIso04","std::vector<float> ",&gsfElectron_hcalIso04);  
      tree->Branch("gsfElectron_pfPhotonIso","std::vector<float> ",&gsfElectron_pfPhotonIso); 
      tree->Branch("gsfElectron_pfChargedHadronIso","std::vector<float> ",&gsfElectron_pfChargedHadronIso); 
      tree->Branch("gsfElectron_pfNeutralHadronIso","std::vector<float> ",&gsfElectron_pfNeutralHadronIso);  
      tree->Branch("gsfElectron_mva_Isolated","std::vector<float> ",&gsfElectron_mva_Isolated);  
      tree->Branch("gsfElectron_mva_e_pi","std::vector<float> ",&gsfElectron_mva_e_pi);  
      tree->Branch("gsfElectron_dnn_signal_Isolated","std::vector<float> ",&gsfElectron_dnn_signal_Isolated);  
      tree->Branch("gsfElectron_dnn_signal_nonIsolated","std::vector<float> ",&gsfElectron_dnn_signal_nonIsolated);  
      tree->Branch("gsfElectron_dnn_bkg_nonIsolated","std::vector<float> ",&gsfElectron_dnn_bkg_nonIsolated);  
      tree->Branch("gsfElectron_dnn_bkg_Tau","std::vector<float> ",&gsfElectron_dnn_bkg_Tau);  
      tree->Branch("gsfElectron_dnn_bkg_Photon","std::vector<float> ",&gsfElectron_dnn_bkg_Photon);  
   }

   if(saveGedPhotons_){
      tree->Branch("gedPhoton_index","std::vector<int> ",&gedPhoton_index);    
      tree->Branch("gedPhoton_seedRawId","std::vector<uint32_t> ",&gedPhoton_seedRawId);   
      tree->Branch("gedPhoton_isEB","std::vector<bool> ",&gedPhoton_isEB);  
      tree->Branch("gedPhoton_isEE","std::vector<bool> ",&gedPhoton_isEE);  
      tree->Branch("gedPhoton_isEBEEGap","std::vector<bool> ",&gedPhoton_isEBEEGap);
      tree->Branch("gedPhoton_isEBEtaGap","std::vector<bool> ",&gedPhoton_isEBEtaGap);
      tree->Branch("gedPhoton_isEBPhiGap","std::vector<bool> ",&gedPhoton_isEBPhiGap);
      tree->Branch("gedPhoton_isEEDeeGap","std::vector<bool> ",&gedPhoton_isEEDeeGap);
      tree->Branch("gedPhoton_isEERingGap","std::vector<bool> ",&gedPhoton_isEERingGap);  
      tree->Branch("gedPhoton_scNPFClusters","std::vector<int> ",&gedPhoton_scNPFClusters); 
      tree->Branch("gedPhoton_ecalSCNPFClusters","std::vector<int> ",&gedPhoton_ecalSCNPFClusters);     
      tree->Branch("gedPhoton_hasConversionTracks","std::vector<bool> ",&gedPhoton_hasConversionTracks); 
      tree->Branch("gedPhoton_nConversions","std::vector<int> ",&gedPhoton_nConversions);     
      tree->Branch("gedPhoton_nConversionsOneLeg","std::vector<int> ",&gedPhoton_nConversionsOneLeg);      
      tree->Branch("gedPhoton_et","std::vector<float> ",&gedPhoton_et); 
      tree->Branch("gedPhoton_energy","std::vector<float> ",&gedPhoton_energy);   
      tree->Branch("gedPhoton_energyErr","std::vector<float> ",&gedPhoton_energyErr);   
      tree->Branch("gedPhoton_ecalEnergy","std::vector<float> ",&gedPhoton_ecalEnergy); 
      tree->Branch("gedPhoton_ecalEnergyErr","std::vector<float> ",&gedPhoton_ecalEnergyErr);  
      tree->Branch("gedPhoton_eta","std::vector<float> ",&gedPhoton_eta); 
      tree->Branch("gedPhoton_phi","std::vector<float> ",&gedPhoton_phi);   
      tree->Branch("gedPhoton_hademTow","std::vector<float> ",&gedPhoton_hademTow); 
      tree->Branch("gedPhoton_hademCone","std::vector<float> ",&gedPhoton_hademCone);  
      tree->Branch("gedPhoton_nrSatCrys","std::vector<float> ",&gedPhoton_nrSatCrys); 
      tree->Branch("gedPhoton_scEta","std::vector<float> ",&gedPhoton_scEta); 
      tree->Branch("gedPhoton_scPhi","std::vector<float> ",&gedPhoton_scPhi); 
      tree->Branch("gedPhoton_scEnergy","std::vector<float> ",&gedPhoton_scEnergy); 
      tree->Branch("gedPhoton_scRawEnergy","std::vector<float> ",&gedPhoton_scRawEnergy); 
      tree->Branch("gedPhoton_scRawESEnergy","std::vector<float> ",&gedPhoton_scRawESEnergy);
      tree->Branch("gedPhoton_scEt","std::vector<float> ",&gedPhoton_scEt); 
      tree->Branch("gedPhoton_scEtaWidth","std::vector<float> ",&gedPhoton_scEtaWidth);   
      tree->Branch("gedPhoton_scPhiWidth","std::vector<float> ",&gedPhoton_scPhiWidth);   
      tree->Branch("gedPhoton_ecalSCEta","std::vector<float> ",&gedPhoton_ecalSCEta); 
      tree->Branch("gedPhoton_ecalSCPhi","std::vector<float> ",&gedPhoton_ecalSCPhi); 
      tree->Branch("gedPhoton_ecalSCEnergy","std::vector<float> ",&gedPhoton_ecalSCEnergy); 
      tree->Branch("gedPhoton_ecalSCRawEnergy","std::vector<float> ",&gedPhoton_ecalSCRawEnergy); 
      tree->Branch("gedPhoton_ecalSCRawESEnergy","std::vector<float> ",&gedPhoton_ecalSCRawESEnergy);
      tree->Branch("gedPhoton_ecalSCEt","std::vector<float> ",&gedPhoton_ecalSCEt); 
      tree->Branch("gedPhoton_ecalSCEtaWidth","std::vector<float> ",&gedPhoton_ecalSCEtaWidth);   
      tree->Branch("gedPhoton_ecalSCPhiWidth","std::vector<float> ",&gedPhoton_ecalSCPhiWidth);    
      tree->Branch("gedPhoton_scSwissCross","std::vector<float> ",&gedPhoton_scSwissCross);  
      tree->Branch("gedPhoton_scEMax","std::vector<float> ",&gedPhoton_scEMax);  
      tree->Branch("gedPhoton_scE2x2","std::vector<float> ",&gedPhoton_scE2x2);  
      tree->Branch("gedPhoton_scE3x3","std::vector<float> ",&gedPhoton_scE3x3);  
      tree->Branch("gedPhoton_scE5x5","std::vector<float> ",&gedPhoton_scE5x5); 
      tree->Branch("gedPhoton_scR9","std::vector<float> ",&gedPhoton_scR9);  
      tree->Branch("gedPhoton_scSigmaIEtaIEta","std::vector<float> ",&gedPhoton_scSigmaIEtaIEta); 
      tree->Branch("gedPhoton_scSigmaIEtaIPhi","std::vector<float> ",&gedPhoton_scSigmaIEtaIPhi);  
      tree->Branch("gedPhoton_scSigmaIPhiIPhi","std::vector<float> ",&gedPhoton_scSigmaIPhiIPhi); 
      tree->Branch("gedPhoton_full5x5_scEMax","std::vector<float> ",&gedPhoton_full5x5_scEMax);  
      tree->Branch("gedPhoton_full5x5_scE2x2","std::vector<float> ",&gedPhoton_full5x5_scE2x2);  
      tree->Branch("gedPhoton_full5x5_scE3x3","std::vector<float> ",&gedPhoton_full5x5_scE3x3);  
      tree->Branch("gedPhoton_full5x5_scE5x5","std::vector<float> ",&gedPhoton_full5x5_scE5x5); 
      tree->Branch("gedPhoton_full5x5_scR9","std::vector<float> ",&gedPhoton_full5x5_scR9);  
      tree->Branch("gedPhoton_full5x5_scSigmaIEtaIEta","std::vector<float> ",&gedPhoton_full5x5_scSigmaIEtaIEta); 
      tree->Branch("gedPhoton_full5x5_scSigmaIEtaIPhi","std::vector<float> ",&gedPhoton_full5x5_scSigmaIEtaIPhi);  
      tree->Branch("gedPhoton_full5x5_scSigmaIPhiIPhi","std::vector<float> ",&gedPhoton_full5x5_scSigmaIPhiIPhi); 
      tree->Branch("gedPhoton_HoE","std::vector<float> ",&gedPhoton_HoE);
      tree->Branch("gedPhoton_trkIso03","std::vector<float> ",&gedPhoton_trkIso03); 
      tree->Branch("gedPhoton_ecalIso03","std::vector<float> ",&gedPhoton_ecalIso03); 
      tree->Branch("gedPhoton_hcalIso03","std::vector<float> ",&gedPhoton_hcalIso03);  
      tree->Branch("gedPhoton_trkIso04","std::vector<float> ",&gedPhoton_trkIso04); 
      tree->Branch("gedPhoton_ecalIso04","std::vector<float> ",&gedPhoton_ecalIso04); 
      tree->Branch("gedPhoton_hcalIso04","std::vector<float> ",&gedPhoton_hcalIso04);  
      tree->Branch("gedPhoton_pfPhotonIso","std::vector<float> ",&gedPhoton_pfPhotonIso); 
      tree->Branch("gedPhoton_pfChargedHadronIso","std::vector<float> ",&gedPhoton_pfChargedHadronIso); 
      tree->Branch("gedPhoton_pfNeutralHadronIso","std::vector<float> ",&gedPhoton_pfNeutralHadronIso);   
      tree->Branch("gedPhoton_nClusterOutsideMustache","std::vector<int> ",&gedPhoton_nClusterOutsideMustache);   
      tree->Branch("gedPhoton_etOutsideMustache","std::vector<float> ",&gedPhoton_etOutsideMustache);   
      tree->Branch("gedPhoton_pfMVA","std::vector<float> ",&gedPhoton_pfMVA);   
      tree->Branch("gedPhoton_pfDNN","std::vector<float> ",&gedPhoton_pfDNN);   
   }

   if(savePatPhotons_ || savePatElectrons_ || savePatJets_){
      tree->Branch("patMET_sumEt", &patMET_sumEt, "patMET_sumEt/F"); 
      tree->Branch("patMET_et", &patMET_et, "patMET_et/F"); 
   }
   if(savePatElectrons_){
      tree->Branch("mll", &mll, "mll/F"); 
      tree->Branch("patElectron_index","std::vector<int> ",&patElectron_index); 
      tree->Branch("patElectron_seedRawId","std::vector<uint32_t> ",&patElectron_seedRawId);   
      tree->Branch("patElectron_classification","std::vector<int> ",&patElectron_classification); 
      tree->Branch("patElectron_scNPFClusters","std::vector<int> ",&patElectron_scNPFClusters); 
      tree->Branch("patElectron_ecalSCNPFClusters","std::vector<int> ",&patElectron_ecalSCNPFClusters);     
      tree->Branch("patElectron_charge","std::vector<int> ",&patElectron_charge); 
      tree->Branch("patElectron_isEB","std::vector<bool> ",&patElectron_isEB); 
      tree->Branch("patElectron_isEE","std::vector<bool> ",&patElectron_isEE); 
      tree->Branch("patElectron_isEBEEGap","std::vector<bool> ",&patElectron_isEBEEGap);
      tree->Branch("patElectron_isEBEtaGap","std::vector<bool> ",&patElectron_isEBEtaGap);
      tree->Branch("patElectron_isEBPhiGap","std::vector<bool> ",&patElectron_isEBPhiGap);
      tree->Branch("patElectron_isEEDeeGap","std::vector<bool> ",&patElectron_isEEDeeGap);
      tree->Branch("patElectron_isEERingGap","std::vector<bool> ",&patElectron_isEERingGap);
      tree->Branch("patElectron_isEcalDriven","std::vector<bool> ",&patElectron_isEcalDriven);
      tree->Branch("patElectron_isTrackerDriven","std::vector<bool> ",&patElectron_isTrackerDriven);
      tree->Branch("patElectron_passConversionVeto","std::vector<bool> ",&patElectron_passConversionVeto);
      tree->Branch("patElectron_nOverlapPhotons","std::vector<int> ",&patElectron_nOverlapPhotons); 
      tree->Branch("patElectron_overlapPhotonIndices","std::vector<std::vector<int> >",&patElectron_overlapPhotonIndices); 
      tree->Branch("patElectron_hasOverlapJet","std::vector<bool> ",&patElectron_hasOverlapJet); 
      tree->Branch("patElectron_overlapJetIndex","std::vector<int> ",&patElectron_overlapJetIndex); 
      tree->Branch("patElectron_eta","std::vector<float> ",&patElectron_eta); 
      tree->Branch("patElectron_phi","std::vector<float> ",&patElectron_phi); 
      tree->Branch("patElectron_p","std::vector<float> ",&patElectron_p); 
      tree->Branch("patElectron_pt","std::vector<float> ",&patElectron_pt); 
      tree->Branch("patElectron_pIn","std::vector<float> ",&patElectron_pIn); 
      tree->Branch("patElectron_pOut","std::vector<float> ",&patElectron_pOut); 
      tree->Branch("patElectron_pAtCalo","std::vector<float> ",&patElectron_pAtCalo); 
      tree->Branch("patElectron_deltaEtaIn","std::vector<float> ",&patElectron_deltaEtaIn); 
      tree->Branch("patElectron_deltaPhiIn","std::vector<float> ",&patElectron_deltaPhiIn); 
      tree->Branch("patElectron_deltaEtaSeedClusterAtCalo","std::vector<float> ",&patElectron_deltaEtaSeedClusterAtCalo);
      tree->Branch("patElectron_deltaPhiSeedClusterAtCalo","std::vector<float> ",&patElectron_deltaPhiSeedClusterAtCalo);  
      tree->Branch("patElectron_deltaEtaEleClusterAtCalo","std::vector<float> ",&patElectron_deltaEtaEleClusterAtCalo); 
      tree->Branch("patElectron_deltaPhiEleClusterAtCalo","std::vector<float> ",&patElectron_deltaPhiEleClusterAtCalo); 
      tree->Branch("patElectron_misHits","std::vector<int> ",&patElectron_misHits); 
      tree->Branch("patElectron_nAmbiguousGsfTracks","std::vector<int> ",&patElectron_nAmbiguousGsfTracks); 
      tree->Branch("patElectron_trackFbrem","std::vector<float> ",&patElectron_trackFbrem);
      tree->Branch("patElectron_superClusterFbrem","std::vector<float> ",&patElectron_superClusterFbrem);
      tree->Branch("patElectron_dz","std::vector<float> ",&patElectron_dz);
      tree->Branch("patElectron_dzError","std::vector<float> ",&patElectron_dzError);
      tree->Branch("patElectron_dxy","std::vector<float> ",&patElectron_dxy);
      tree->Branch("patElectron_dxyError","std::vector<float> ",&patElectron_dxyError);
      tree->Branch("patElectron_energy","std::vector<float> ",&patElectron_energy);   
      tree->Branch("patElectron_energyErr","std::vector<float> ",&patElectron_energyErr);   
      tree->Branch("patElectron_ecalEnergy","std::vector<float> ",&patElectron_ecalEnergy); 
      tree->Branch("patElectron_ecalEnergyErr","std::vector<float> ",&patElectron_ecalEnergyErr);  
      tree->Branch("patElectron_et","std::vector<float> ",&patElectron_et); 
      tree->Branch("patElectron_mt","std::vector<float> ",&patElectron_mt); 
      tree->Branch("patElectron_dphiMET","std::vector<float> ",&patElectron_dphiMET); 
      tree->Branch("patElectron_scEta","std::vector<float> ",&patElectron_scEta); 
      tree->Branch("patElectron_scPhi","std::vector<float> ",&patElectron_scPhi);
      tree->Branch("patElectron_scEnergy","std::vector<float> ",&patElectron_scEnergy); 
      tree->Branch("patElectron_scRawEnergy","std::vector<float> ",&patElectron_scRawEnergy); 
      tree->Branch("patElectron_scRawESEnergy","std::vector<float> ",&patElectron_scRawESEnergy); 
      tree->Branch("patElectron_scEt","std::vector<float> ",&patElectron_scEt); 
      tree->Branch("patElectron_scPhiWidth","std::vector<float> ",&patElectron_scPhiWidth);  
      tree->Branch("patElectron_scEtaWidth","std::vector<float> ",&patElectron_scEtaWidth);  
      tree->Branch("patElectron_scEoP","std::vector<float> ",&patElectron_scEoP);  
      tree->Branch("patElectron_ecalSCEta","std::vector<float> ",&patElectron_ecalSCEta); 
      tree->Branch("patElectron_ecalSCPhi","std::vector<float> ",&patElectron_ecalSCPhi); 
      tree->Branch("patElectron_ecalSCEnergy","std::vector<float> ",&patElectron_ecalSCEnergy); 
      tree->Branch("patElectron_ecalSCRawEnergy","std::vector<float> ",&patElectron_ecalSCRawEnergy); 
      tree->Branch("patElectron_ecalSCRawESEnergy","std::vector<float> ",&patElectron_ecalSCRawESEnergy); 
      tree->Branch("patElectron_ecalSCEt","std::vector<float> ",&patElectron_ecalSCEt); 
      tree->Branch("patElectron_ecalSCPhiWidth","std::vector<float> ",&patElectron_ecalSCPhiWidth);  
      tree->Branch("patElectron_ecalSCEtaWidth","std::vector<float> ",&patElectron_ecalSCEtaWidth);  
      tree->Branch("patElectron_ecalSCEoP","std::vector<float> ",&patElectron_ecalSCEoP);   
      tree->Branch("patElectron_scSwissCross","std::vector<float> ",&patElectron_scSwissCross);  
      tree->Branch("patElectron_scE2x2","std::vector<float> ",&patElectron_scE2x2);  
      tree->Branch("patElectron_scE3x3","std::vector<float> ",&patElectron_scE3x3);  
      tree->Branch("patElectron_scE5x5","std::vector<float> ",&patElectron_scE5x5);  
      tree->Branch("patElectron_scEMax","std::vector<float> ",&patElectron_scEMax); 
      tree->Branch("patElectron_scR9","std::vector<float> ",&patElectron_scR9); 
      tree->Branch("patElectron_scSigmaIEtaIEta","std::vector<float> ",&patElectron_scSigmaIEtaIEta);  
      tree->Branch("patElectron_scSigmaIEtaIPhi","std::vector<float> ",&patElectron_scSigmaIEtaIPhi);  
      tree->Branch("patElectron_scSigmaIPhiIPhi","std::vector<float> ",&patElectron_scSigmaIPhiIPhi);  
      tree->Branch("patElectron_full5x5_scE2x2","std::vector<float> ",&patElectron_full5x5_scE2x2);  
      tree->Branch("patElectron_full5x5_scE3x3","std::vector<float> ",&patElectron_full5x5_scE3x3);  
      tree->Branch("patElectron_full5x5_scE5x5","std::vector<float> ",&patElectron_full5x5_scE5x5);  
      tree->Branch("patElectron_full5x5_scEMax","std::vector<float> ",&patElectron_full5x5_scEMax);  
      tree->Branch("patElectron_full5x5_scR9","std::vector<float> ",&patElectron_full5x5_scR9);  
      tree->Branch("patElectron_full5x5_scSigmaIEtaIEta","std::vector<float> ",&patElectron_full5x5_scSigmaIEtaIEta);  
      tree->Branch("patElectron_full5x5_scSigmaIEtaIPhi","std::vector<float> ",&patElectron_full5x5_scSigmaIEtaIPhi);  
      tree->Branch("patElectron_full5x5_scSigmaIPhiIPhi","std::vector<float> ",&patElectron_full5x5_scSigmaIPhiIPhi);  
      tree->Branch("patElectron_HoE","std::vector<float> ",&patElectron_HoE);   
      tree->Branch("patElectron_trkIso03","std::vector<float> ",&patElectron_trkIso03); 
      tree->Branch("patElectron_ecalIso03","std::vector<float> ",&patElectron_ecalIso03); 
      tree->Branch("patElectron_hcalIso03","std::vector<float> ",&patElectron_hcalIso03);  
      tree->Branch("patElectron_trkIso04","std::vector<float> ",&patElectron_trkIso04); 
      tree->Branch("patElectron_ecalIso04","std::vector<float> ",&patElectron_ecalIso04); 
      tree->Branch("patElectron_hcalIso04","std::vector<float> ",&patElectron_hcalIso04);  
      tree->Branch("patElectron_pfPhotonIso","std::vector<float> ",&patElectron_pfPhotonIso); 
      tree->Branch("patElectron_pfChargedHadronIso","std::vector<float> ",&patElectron_pfChargedHadronIso); 
      tree->Branch("patElectron_pfNeutralHadronIso","std::vector<float> ",&patElectron_pfNeutralHadronIso);  
      tree->Branch("patElectron_mva_Isolated","std::vector<float> ",&patElectron_mva_Isolated);  
      tree->Branch("patElectron_mva_e_pi","std::vector<float> ",&patElectron_mva_e_pi);  
      tree->Branch("patElectron_dnn_signal_Isolated","std::vector<float> ",&patElectron_dnn_signal_Isolated);  
      tree->Branch("patElectron_dnn_signal_nonIsolated","std::vector<float> ",&patElectron_dnn_signal_nonIsolated);  
      tree->Branch("patElectron_dnn_bkg_nonIsolated","std::vector<float> ",&patElectron_dnn_bkg_nonIsolated);  
      tree->Branch("patElectron_dnn_bkg_Tau","std::vector<float> ",&patElectron_dnn_bkg_Tau);  
      tree->Branch("patElectron_dnn_bkg_Photon","std::vector<float> ",&patElectron_dnn_bkg_Photon);  
      tree->Branch("patElectron_egmCutBasedElectronIDVeto","std::vector<int> ",&patElectron_egmCutBasedElectronIDVeto);  
      tree->Branch("patElectron_egmCutBasedElectronIDloose","std::vector<int> ",&patElectron_egmCutBasedElectronIDloose);  
      tree->Branch("patElectron_egmCutBasedElectronIDmedium","std::vector<int> ",&patElectron_egmCutBasedElectronIDmedium);
      tree->Branch("patElectron_egmCutBasedElectronIDtight","std::vector<int> ",&patElectron_egmCutBasedElectronIDtight);  
      tree->Branch("patElectron_egmMVAElectronIDloose","std::vector<int> ",&patElectron_egmMVAElectronIDloose);  
      tree->Branch("patElectron_egmMVAElectronIDmedium","std::vector<int> ",&patElectron_egmMVAElectronIDmedium);
      tree->Branch("patElectron_egmMVAElectronIDtight","std::vector<int> ",&patElectron_egmMVAElectronIDtight);  
      tree->Branch("patElectron_egmMVAElectronIDlooseNoIso","std::vector<int> ",&patElectron_egmMVAElectronIDlooseNoIso);  
      tree->Branch("patElectron_egmMVAElectronIDmediumNoIso","std::vector<int> ",&patElectron_egmMVAElectronIDmediumNoIso);
      tree->Branch("patElectron_egmMVAElectronIDtightNoIso","std::vector<int> ",&patElectron_egmMVAElectronIDtightNoIso);  
      tree->Branch("patElectron_heepElectronID","std::vector<int> ",&patElectron_heepElectronID);  
   }
   if(savePatPhotons_){
      tree->Branch("patPhoton_index","std::vector<int> ",&patPhoton_index); 
      tree->Branch("patPhoton_seedRawId","std::vector<uint32_t> ",&patPhoton_seedRawId);   
      tree->Branch("patPhoton_scNPFClusters","std::vector<int> ",&patPhoton_scNPFClusters); 
      tree->Branch("patPhoton_ecalSCNPFClusters","std::vector<int> ",&patPhoton_ecalSCNPFClusters);     
      tree->Branch("patPhoton_passElectronVeto","std::vector<bool> ",&patPhoton_passElectronVeto); 
      tree->Branch("patPhoton_hasPixelSeed","std::vector<bool> ",&patPhoton_hasPixelSeed); 
      tree->Branch("patPhoton_hasConversionTracks","std::vector<bool> ",&patPhoton_hasConversionTracks); 
      tree->Branch("patPhoton_nConversions","std::vector<int> ",&patPhoton_nConversions);     
      tree->Branch("patPhoton_nConversionsOneLeg","std::vector<int> ",&patPhoton_nConversionsOneLeg);     
      tree->Branch("patPhoton_isEB","std::vector<bool> ",&patPhoton_isEB); 
      tree->Branch("patPhoton_isEE","std::vector<bool> ",&patPhoton_isEE); 
      tree->Branch("patPhoton_isEBEEGap","std::vector<bool> ",&patPhoton_isEBEEGap);
      tree->Branch("patPhoton_isEBEtaGap","std::vector<bool> ",&patPhoton_isEBEtaGap);
      tree->Branch("patPhoton_isEBPhiGap","std::vector<bool> ",&patPhoton_isEBPhiGap);
      tree->Branch("patPhoton_isEEDeeGap","std::vector<bool> ",&patPhoton_isEEDeeGap);
      tree->Branch("patPhoton_isEERingGap","std::vector<bool> ",&patPhoton_isEERingGap);  
      tree->Branch("patPhoton_hasOverlapElectron","std::vector<bool> ",&patPhoton_hasOverlapElectron); 
      tree->Branch("patPhoton_overlapElectronIndex","std::vector<int> ",&patPhoton_overlapElectronIndex); 
      tree->Branch("patPhoton_hasOverlapJet","std::vector<bool> ",&patPhoton_hasOverlapJet); 
      tree->Branch("patPhoton_overlapJetIndex","std::vector<int> ",&patPhoton_overlapJetIndex);  
      tree->Branch("patPhoton_eta","std::vector<float> ",&patPhoton_eta); 
      tree->Branch("patPhoton_phi","std::vector<float> ",&patPhoton_phi); 
      tree->Branch("patPhoton_energy","std::vector<float> ",&patPhoton_energy);   
      tree->Branch("patPhoton_energyErr","std::vector<float> ",&patPhoton_energyErr);   
      tree->Branch("patPhoton_ecalEnergy","std::vector<float> ",&patPhoton_ecalEnergy); 
      tree->Branch("patPhoton_ecalEnergyErr","std::vector<float> ",&patPhoton_ecalEnergyErr);  
      tree->Branch("patPhoton_et","std::vector<float> ",&patPhoton_et); 
      tree->Branch("patPhoton_mt","std::vector<float> ",&patPhoton_mt); 
      tree->Branch("patPhoton_dphiMET","std::vector<float> ",&patPhoton_dphiMET); 
      tree->Branch("patPhoton_scEta","std::vector<float> ",&patPhoton_scEta); 
      tree->Branch("patPhoton_scPhi","std::vector<float> ",&patPhoton_scPhi);  
      tree->Branch("patPhoton_scEnergy","std::vector<float> ",&patPhoton_scEnergy); 
      tree->Branch("patPhoton_scRawEnergy","std::vector<float> ",&patPhoton_scRawEnergy); 
      tree->Branch("patPhoton_scRawESEnergy","std::vector<float> ",&patPhoton_scRawESEnergy); 
      tree->Branch("patPhoton_scEt","std::vector<float> ",&patPhoton_scEt); 
      tree->Branch("patPhoton_scPhiWidth","std::vector<float> ",&patPhoton_scPhiWidth);  
      tree->Branch("patPhoton_scEtaWidth","std::vector<float> ",&patPhoton_scEtaWidth);  
      tree->Branch("patPhoton_ecalSCEta","std::vector<float> ",&patPhoton_ecalSCEta); 
      tree->Branch("patPhoton_ecalSCPhi","std::vector<float> ",&patPhoton_ecalSCPhi);
      tree->Branch("patPhoton_ecalSCEnergy","std::vector<float> ",&patPhoton_ecalSCEnergy); 
      tree->Branch("patPhoton_ecalSCRawEnergy","std::vector<float> ",&patPhoton_ecalSCRawEnergy); 
      tree->Branch("patPhoton_ecalSCRawESEnergy","std::vector<float> ",&patPhoton_ecalSCRawESEnergy); 
      tree->Branch("patPhoton_ecalSCEt","std::vector<float> ",&patPhoton_ecalSCEt); 
      tree->Branch("patPhoton_ecalSCPhiWidth","std::vector<float> ",&patPhoton_ecalSCPhiWidth);  
      tree->Branch("patPhoton_ecalSCEtaWidth","std::vector<float> ",&patPhoton_ecalSCEtaWidth); 
      tree->Branch("patPhoton_scSwissCross","std::vector<float> ",&patPhoton_scSwissCross);  
      tree->Branch("patPhoton_scE2x2","std::vector<float> ",&patPhoton_scE2x2);  
      tree->Branch("patPhoton_scE3x3","std::vector<float> ",&patPhoton_scE3x3);  
      tree->Branch("patPhoton_scE5x5","std::vector<float> ",&patPhoton_scE5x5);  
      tree->Branch("patPhoton_scEMax","std::vector<float> ",&patPhoton_scEMax); 
      tree->Branch("patPhoton_scR9","std::vector<float> ",&patPhoton_scR9); 
      tree->Branch("patPhoton_scSigmaIEtaIEta","std::vector<float> ",&patPhoton_scSigmaIEtaIEta);  
      tree->Branch("patPhoton_scSigmaIEtaIPhi","std::vector<float> ",&patPhoton_scSigmaIEtaIPhi);  
      tree->Branch("patPhoton_scSigmaIPhiIPhi","std::vector<float> ",&patPhoton_scSigmaIPhiIPhi);  
      tree->Branch("patPhoton_full5x5_scE2x2","std::vector<float> ",&patPhoton_full5x5_scE2x2);  
      tree->Branch("patPhoton_full5x5_scE3x3","std::vector<float> ",&patPhoton_full5x5_scE3x3);  
      tree->Branch("patPhoton_full5x5_scE5x5","std::vector<float> ",&patPhoton_full5x5_scE5x5);  
      tree->Branch("patPhoton_full5x5_scEMax","std::vector<float> ",&patPhoton_full5x5_scEMax);  
      tree->Branch("patPhoton_full5x5_scR9","std::vector<float> ",&patPhoton_full5x5_scR9);  
      tree->Branch("patPhoton_full5x5_scSigmaIEtaIEta","std::vector<float> ",&patPhoton_full5x5_scSigmaIEtaIEta);  
      tree->Branch("patPhoton_full5x5_scSigmaIEtaIPhi","std::vector<float> ",&patPhoton_full5x5_scSigmaIEtaIPhi);  
      tree->Branch("patPhoton_full5x5_scSigmaIPhiIPhi","std::vector<float> ",&patPhoton_full5x5_scSigmaIPhiIPhi);  
      tree->Branch("patPhoton_HoE","std::vector<float> ",&patPhoton_HoE);   
      tree->Branch("patPhoton_trkIso03","std::vector<float> ",&patPhoton_trkIso03);  
      tree->Branch("patPhoton_ecalIso03","std::vector<float> ",&patPhoton_ecalIso03);  
      tree->Branch("patPhoton_hcalIso03","std::vector<float> ",&patPhoton_hcalIso03);  
      tree->Branch("patPhoton_trkIso04","std::vector<float> ",&patPhoton_trkIso04);  
      tree->Branch("patPhoton_ecalIso04","std::vector<float> ",&patPhoton_ecalIso04);  
      tree->Branch("patPhoton_hcalIso04","std::vector<float> ",&patPhoton_hcalIso04);  
      tree->Branch("patPhoton_patParticleIso","std::vector<float> ",&patPhoton_patParticleIso);  
      tree->Branch("patPhoton_pfChargedHadronIso","std::vector<float> ",&patPhoton_pfChargedHadronIso);  
      tree->Branch("patPhoton_pfNeutralHadronIso","std::vector<float> ",&patPhoton_pfNeutralHadronIso);   
      tree->Branch("patPhoton_pfPhotonIso","std::vector<float> ",&patPhoton_pfPhotonIso);  
      tree->Branch("patPhoton_pfPuChargedHadronIso","std::vector<float> ",&patPhoton_pfPuChargedHadronIso);   
      tree->Branch("patPhoton_nClusterOutsideMustache","std::vector<int> ",&patPhoton_nClusterOutsideMustache);   
      tree->Branch("patPhoton_etOutsideMustache","std::vector<float> ",&patPhoton_etOutsideMustache);   
      tree->Branch("patPhoton_pfMVA","std::vector<float> ",&patPhoton_pfMVA);   
      tree->Branch("patPhoton_pfDNN","std::vector<float> ",&patPhoton_pfDNN);
      tree->Branch("patPhoton_egmCutBasedPhotonIDloose","std::vector<int> ",&patPhoton_egmCutBasedPhotonIDloose);  
      tree->Branch("patPhoton_egmCutBasedPhotonIDmedium","std::vector<int> ",&patPhoton_egmCutBasedPhotonIDmedium);
      tree->Branch("patPhoton_egmCutBasedPhotonIDtight","std::vector<int> ",&patPhoton_egmCutBasedPhotonIDtight);   
      tree->Branch("patPhoton_egmMVAPhotonIDmedium","std::vector<int> ",&patPhoton_egmMVAPhotonIDmedium);
      tree->Branch("patPhoton_egmMVAPhotonIDtight","std::vector<int> ",&patPhoton_egmMVAPhotonIDtight);  
   }

   if(savePatJets_){
      tree->Branch("patJet_index","std::vector<int> ",&patJet_index); 
      tree->Branch("patJet_isCaloJet","std::vector<bool> ",&patJet_isCaloJet); 
      tree->Branch("patJet_isJPTJet","std::vector<bool> ",&patJet_isJPTJet); 
      tree->Branch("patJet_isPFJet","std::vector<bool> ",&patJet_isPFJet); 
      tree->Branch("patJet_isBasicJet","std::vector<bool> ",&patJet_isBasicJet); 
      tree->Branch("patJet_charge","std::vector<float> ",&patJet_charge);  
      tree->Branch("patJet_energy","std::vector<float> ",&patJet_energy);
      tree->Branch("patJet_uncorrectedEnergy","std::vector<float> ",&patJet_uncorrectedEnergy);    
      tree->Branch("patJet_pt","std::vector<float> ",&patJet_pt); 
      tree->Branch("patJet_uncorrectedPt","std::vector<float> ",&patJet_uncorrectedPt); 
      tree->Branch("patJet_eta","std::vector<float> ",&patJet_eta);  
      tree->Branch("patJet_phi","std::vector<float> ",&patJet_phi);  
      tree->Branch("patJet_area","std::vector<float> ",&patJet_area);  
      tree->Branch("patJet_energyFractionHadronic","std::vector<float> ",&patJet_energyFractionHadronic);  
      tree->Branch("patJet_hadEnergyInHB","std::vector<float> ",&patJet_hadEnergyInHB);  
      tree->Branch("patJet_hadEnergyInHO","std::vector<float> ",&patJet_hadEnergyInHO); 
      tree->Branch("patJet_hadEnergyInHE","std::vector<float> ",&patJet_hadEnergyInHE);  
      tree->Branch("patJet_hadEnergyInHF","std::vector<float> ",&patJet_hadEnergyInHF);  
      tree->Branch("patJet_emEnergyInEB","std::vector<float> ",&patJet_emEnergyInEB);  
      tree->Branch("patJet_emEnergyInEE","std::vector<float> ",&patJet_emEnergyInEE); 
      tree->Branch("patJet_emEnergyInHF","std::vector<float> ",&patJet_emEnergyInHF);  
      tree->Branch("patJet_chargedHadronEnergyFraction","std::vector<float> ",&patJet_chargedHadronEnergyFraction);  
      tree->Branch("patJet_neutralHadronEnergyFraction","std::vector<float> ",&patJet_neutralHadronEnergyFraction);  
      tree->Branch("patJet_chargedEmEnergyFraction","std::vector<float> ",&patJet_chargedEmEnergyFraction); 
      tree->Branch("patJet_neutralEmEnergyFraction","std::vector<float> ",&patJet_neutralEmEnergyFraction);  
      tree->Branch("patJet_nOverlapMuons","std::vector<int> ",&patJet_nOverlapMuons); 
      tree->Branch("patJet_nOverlapTaus","std::vector<int> ",&patJet_nOverlapTaus); 
      tree->Branch("patJet_nOverlapElectrons","std::vector<int> ",&patJet_nOverlapElectrons); 
      tree->Branch("patJet_overlapElectronIndices","std::vector<std::vector<int> >",&patJet_overlapElectronIndices); 
      tree->Branch("patJet_nOverlapPhotons","std::vector<int> ",&patJet_nOverlapPhotons); 
      tree->Branch("patJet_overlapPhotonIndices","std::vector<std::vector<int> >",&patJet_overlapPhotonIndices); 
      tree->Branch("patJet_photonEnergy","std::vector<float> ",&patJet_photonEnergy);  
      tree->Branch("patJet_photonEnergyFraction","std::vector<float> ",&patJet_photonEnergyFraction);  
      tree->Branch("patJet_electronEnergy","std::vector<float> ",&patJet_electronEnergy); 
      tree->Branch("patJet_electronEnergyFraction","std::vector<float> ",&patJet_electronEnergyFraction);  
      tree->Branch("patJet_muonEnergy","std::vector<float> ",&patJet_muonEnergy);  
      tree->Branch("patJet_muonEnergyFraction","std::vector<float> ",&patJet_muonEnergyFraction);  
      tree->Branch("patJet_HFHadronEnergy","std::vector<float> ",&patJet_HFHadronEnergy); 
      tree->Branch("patJet_HFHadronEnergyFraction","std::vector<float> ",&patJet_HFHadronEnergyFraction);  
      tree->Branch("patJet_chargedMuEnergy","std::vector<float> ",&patJet_chargedMuEnergy);  
      tree->Branch("patJet_chargedMuEnergyFraction","std::vector<float> ",&patJet_chargedMuEnergyFraction);  
      tree->Branch("patJet_hoEnergy","std::vector<float> ",&patJet_hoEnergy);   
      tree->Branch("patJet_hoEnergyFraction","std::vector<float> ",&patJet_hoEnergyFraction);  
      tree->Branch("patJet_nCandidates","std::vector<int> ",&patJet_nCandidates);
      tree->Branch("patJet_nCandInEcal","std::vector<int> ",&patJet_nCandInEcal);
      tree->Branch("patJet_nCandInEcalWithCharge","std::vector<int> ",&patJet_nCandInEcalWithCharge);
      tree->Branch("patJet_candInEcal_charge","std::vector<std::vector<float> >",&patJet_candInEcal_charge);
      tree->Branch("patJet_candInEcal_ecalEnergy","std::vector<std::vector<float> >",&patJet_candInEcal_ecalEnergy);
      tree->Branch("patJet_candInEcal_ecalEnergyFraction","std::vector<std::vector<float> >",&patJet_candInEcal_ecalEnergyFraction);
      tree->Branch("patJet_candInEcal_hcalEnergy","std::vector<std::vector<float> >",&patJet_candInEcal_hcalEnergy);
      tree->Branch("patJet_candInEcal_hcalEnergyFraction","std::vector<std::vector<float> >",&patJet_candInEcal_hcalEnergyFraction);
      tree->Branch("patJet_candInEcal_eta","std::vector<std::vector<float> >",&patJet_candInEcal_eta);
      tree->Branch("patJet_candInEcal_phi","std::vector<std::vector<float> >",&patJet_candInEcal_phi);
      tree->Branch("patJet_bTagScore_pfDeepCSV","std::vector<float> ",&patJet_bTagScore_pfDeepCSV);   
      tree->Branch("patJet_puIDScore","std::vector<float> ",&patJet_puIDScore); 
      tree->Branch("patJet_puID","std::vector<int> ",&patJet_puID);     
      tree->Branch("patJet_qgL","std::vector<float> ",&patJet_qgL);     
      tree->Branch("patJet_jmeCutBasedPFJetIDloose","std::vector<bool> ",&patJet_jmeCutBasedPFJetIDloose);   
      tree->Branch("patJet_jmeCutBasedPFJetIDtight","std::vector<bool> ",&patJet_jmeCutBasedPFJetIDtight);  
      tree->Branch("patJet_jmeCutBasedPFJetIDtightLepVeto","std::vector<bool> ",&patJet_jmeCutBasedPFJetIDtightLepVeto);       
   } 
 
   if(saveSuperCluster_){
      tree->Branch("superCluster_seedRawId","std::vector<uint32_t> ",&superCluster_seedRawId);     
      tree->Branch("superCluster_rawEnergy","std::vector<float> ",&superCluster_rawEnergy);     
      tree->Branch("superCluster_rawESEnergy","std::vector<float> ",&superCluster_rawESEnergy);     
      tree->Branch("superCluster_energy","std::vector<float> ",&superCluster_energy);
      tree->Branch("superCluster_eta","std::vector<float>",&superCluster_eta);
      tree->Branch("superCluster_phi","std::vector<float>",&superCluster_phi);  
      tree->Branch("superCluster_etaWidth","std::vector<float>",&superCluster_etaWidth);
      tree->Branch("superCluster_phiWidth","std::vector<float>",&superCluster_phiWidth);  
      tree->Branch("superCluster_R","std::vector<float>",&superCluster_R);   
      tree->Branch("superCluster_nPFClusters","std::vector<int>",&superCluster_nPFClusters);   
      tree->Branch("superCluster_ieta","std::vector<int>",&superCluster_ieta);
      tree->Branch("superCluster_iphi","std::vector<int>",&superCluster_iphi);  
      tree->Branch("superCluster_iz","std::vector<int>",&superCluster_iz);   
      tree->Branch("superCluster_nXtals","std::vector<int>",&superCluster_nXtals);   
      if(savePFCluster_) tree->Branch("superCluster_seedIndex","std::vector<int>",&superCluster_seedIndex);     
      if(savePFCluster_) tree->Branch("superCluster_pfClustersIndex","std::vector<std::vector<int> >",&superCluster_pfClustersIndex); 
      tree->Branch("superCluster_psCluster_energy", "std::vector<std::vector<float> >", &superCluster_psCluster_energy);
      tree->Branch("superCluster_psCluster_eta", "std::vector<std::vector<float> >", &superCluster_psCluster_eta);
      tree->Branch("superCluster_psCluster_phi", "std::vector<std::vector<float> >", &superCluster_psCluster_phi); 
      if(saveCaloParticles_ && isMC_){
         tree->Branch("superCluster_dR_genScore","std::vector<std::vector<double> >",&superCluster_dR_genScore);
         tree->Branch("superCluster_dR_simScore","std::vector<std::vector<double> >",&superCluster_dR_simScore);
         tree->Branch("superCluster_sim_nSharedXtals","std::vector<std::vector<double> >",&superCluster_sim_nSharedXtals);
         tree->Branch("superCluster_sim_fraction_noHitsFraction","std::vector<std::vector<double> >",&superCluster_sim_fraction_noHitsFraction); 
         tree->Branch("superCluster_sim_fraction","std::vector<std::vector<double> >",&superCluster_sim_fraction);
         tree->Branch("superCluster_recoToSim_fraction","std::vector<std::vector<double> >",&superCluster_recoToSim_fraction);
         tree->Branch("superCluster_recoToSim_fraction_sharedXtals","std::vector<std::vector<double> >",&superCluster_recoToSim_fraction_sharedXtals);
         tree->Branch("superCluster_simEnergy_sharedXtals","std::vector<std::vector<double> >",&superCluster_simEnergy_sharedXtals);
         tree->Branch("superCluster_recoEnergy_sharedXtals","std::vector<std::vector<double> >",&superCluster_recoEnergy_sharedXtals);  
         tree->Branch("superCluster_dR_genScore_MatchedIndex","std::vector<int>",&superCluster_dR_genScore_MatchedIndex);
         tree->Branch("superCluster_dR_simScore_MatchedIndex","std::vector<int>",&superCluster_dR_simScore_MatchedIndex);
         tree->Branch("superCluster_sim_nSharedXtals_MatchedIndex","std::vector<int>",&superCluster_sim_nSharedXtals_MatchedIndex);
         tree->Branch("superCluster_sim_fraction_noHitsFraction_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_noHitsFraction_MatchedIndex); 
         tree->Branch("superCluster_sim_fraction_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_MatchedIndex);
         tree->Branch("superCluster_recoToSim_fraction_MatchedIndex","std::vector<int>",&superCluster_recoToSim_fraction_MatchedIndex);
         tree->Branch("superCluster_recoToSim_fraction_sharedXtals_MatchedIndex","std::vector<int>",&superCluster_recoToSim_fraction_sharedXtals_MatchedIndex);
         tree->Branch("superCluster_simEnergy_sharedXtals_MatchedIndex","std::vector<int>",&superCluster_simEnergy_sharedXtals_MatchedIndex);
         tree->Branch("superCluster_recoEnergy_sharedXtals_MatchedIndex","std::vector<int>",&superCluster_recoEnergy_sharedXtals_MatchedIndex);
      }    
      if(saveCaloParticlesPU_ && isMC_){ 
         tree->Branch("superCluster_simPU_nSharedXtals","std::vector<double>",&superCluster_simPU_nSharedXtals);
         tree->Branch("superCluster_simEnergy_sharedXtalsPU","std::vector<double>",&superCluster_simEnergy_sharedXtalsPU);
         tree->Branch("superCluster_recoEnergy_sharedXtalsPU","std::vector<double>",&superCluster_recoEnergy_sharedXtalsPU);  
         tree->Branch("superCluster_simEnergy_noHitsFraction_sharedXtalsPU","std::vector<double>",&superCluster_simEnergy_noHitsFraction_sharedXtalsPU);
         tree->Branch("superCluster_recoEnergy_noHitsFraction_sharedXtalsPU","std::vector<double>",&superCluster_recoEnergy_noHitsFraction_sharedXtalsPU);  
      }
      if(saveCaloParticlesOOTPU_ && isMC_){ 
         tree->Branch("superCluster_simOOTPU_nSharedXtals","std::vector<double>",&superCluster_simOOTPU_nSharedXtals);
         tree->Branch("superCluster_simEnergy_sharedXtalsOOTPU","std::vector<double>",&superCluster_simEnergy_sharedXtalsOOTPU);
         tree->Branch("superCluster_recoEnergy_sharedXtalsOOTPU","std::vector<double>",&superCluster_recoEnergy_sharedXtalsOOTPU);  
         tree->Branch("superCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU","std::vector<double>",&superCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU);
         tree->Branch("superCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU","std::vector<double>",&superCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU);  
      }
   }
   if(saveRetunedSC_){   
      tree->Branch("retunedSuperCluster_seedRawId","std::vector<uint32_t> ",&retunedSuperCluster_seedRawId);  
      tree->Branch("retunedSuperCluster_rawEnergy","std::vector<float> ",&retunedSuperCluster_rawEnergy);
      tree->Branch("retunedSuperCluster_rawESEnergy","std::vector<float> ",&retunedSuperCluster_rawESEnergy);
      tree->Branch("retunedSuperCluster_energy","std::vector<float> ",&retunedSuperCluster_energy);
      tree->Branch("retunedSuperCluster_eta","std::vector<float>",&retunedSuperCluster_eta);
      tree->Branch("retunedSuperCluster_phi","std::vector<float>",&retunedSuperCluster_phi);  
      tree->Branch("retunedSuperCluster_etaWidth","std::vector<float>",&retunedSuperCluster_etaWidth);
      tree->Branch("retunedSuperCluster_phiWidth","std::vector<float>",&retunedSuperCluster_phiWidth);  
      tree->Branch("retunedSuperCluster_R","std::vector<float>",&retunedSuperCluster_R);   
      tree->Branch("retunedSuperCluster_nPFClusters","std::vector<int>",&retunedSuperCluster_nPFClusters);   
      tree->Branch("retunedSuperCluster_ieta","std::vector<int>",&retunedSuperCluster_ieta);
      tree->Branch("retunedSuperCluster_iphi","std::vector<int>",&retunedSuperCluster_iphi);  
      tree->Branch("retunedSuperCluster_iz","std::vector<int>",&retunedSuperCluster_iz);    
      tree->Branch("retunedSuperCluster_nXtals","std::vector<int>",&retunedSuperCluster_nXtals);  
      if(savePFCluster_) tree->Branch("retunedSuperCluster_seedIndex","std::vector<int>",&retunedSuperCluster_seedIndex);     
      if(savePFCluster_) tree->Branch("retunedSuperCluster_pfClustersIndex","std::vector<std::vector<int> >",&retunedSuperCluster_pfClustersIndex); 
      tree->Branch("retunedSuperCluster_psCluster_energy", "std::vector<std::vector<float> >", &retunedSuperCluster_psCluster_energy);
      tree->Branch("retunedSuperCluster_psCluster_eta", "std::vector<std::vector<float> >", &retunedSuperCluster_psCluster_eta);
      tree->Branch("retunedSuperCluster_psCluster_phi", "std::vector<std::vector<float> >", &retunedSuperCluster_psCluster_phi); 
      if(saveCaloParticles_ && isMC_){
         tree->Branch("retunedSuperCluster_dR_genScore","std::vector<std::vector<double> >",&retunedSuperCluster_dR_genScore);
         tree->Branch("retunedSuperCluster_dR_simScore","std::vector<std::vector<double> >",&retunedSuperCluster_dR_simScore);
         tree->Branch("retunedSuperCluster_sim_nSharedXtals","std::vector<std::vector<double> >",&retunedSuperCluster_sim_nSharedXtals);
         tree->Branch("retunedSuperCluster_sim_fraction_noHitsFraction","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction_noHitsFraction); 
         tree->Branch("retunedSuperCluster_sim_fraction","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction);
         tree->Branch("retunedSuperCluster_recoToSim_fraction","std::vector<std::vector<double> >",&retunedSuperCluster_recoToSim_fraction);
         tree->Branch("retunedSuperCluster_recoToSim_fraction_sharedXtals","std::vector<std::vector<double> >",&retunedSuperCluster_recoToSim_fraction_sharedXtals);
         tree->Branch("retunedSuperCluster_simEnergy_sharedXtals","std::vector<std::vector<double> >",&retunedSuperCluster_simEnergy_sharedXtals);
         tree->Branch("retunedSuperCluster_recoEnergy_sharedXtals","std::vector<std::vector<double> >",&retunedSuperCluster_recoEnergy_sharedXtals);  
         tree->Branch("retunedSuperCluster_dR_genScore_MatchedIndex","std::vector<int>",&retunedSuperCluster_dR_genScore_MatchedIndex);
         tree->Branch("retunedSuperCluster_dR_simScore_MatchedIndex","std::vector<int>",&retunedSuperCluster_dR_simScore_MatchedIndex);
         tree->Branch("retunedSuperCluster_sim_nSharedXtals_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_nSharedXtals_MatchedIndex);
         tree->Branch("retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<int>", &retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex); 
         tree->Branch("retunedSuperCluster_sim_fraction_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_fraction_MatchedIndex);
         tree->Branch("retunedSuperCluster_recoToSim_fraction_MatchedIndex","std::vector<int>",&retunedSuperCluster_recoToSim_fraction_MatchedIndex);
         tree->Branch("retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<int>", &retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);
         tree->Branch("retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex","std::vector<int>",&retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex);
         tree->Branch("retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex","std::vector<int>",&retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex);
      } 
      if(saveCaloParticlesPU_ && isMC_){ 
         tree->Branch("retunedSuperCluster_simPU_nSharedXtals","std::vector<double>",&retunedSuperCluster_simPU_nSharedXtals);
         tree->Branch("retunedSuperCluster_simEnergy_sharedXtalsPU","std::vector<double>",&retunedSuperCluster_simEnergy_sharedXtalsPU);
         tree->Branch("retunedSuperCluster_recoEnergy_sharedXtalsPU","std::vector<double>",&retunedSuperCluster_recoEnergy_sharedXtalsPU);  
         tree->Branch("retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU","std::vector<double>",&retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU);
         tree->Branch("retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU","std::vector<double>",&retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU);  
      }
      if(saveCaloParticlesOOTPU_ && isMC_){ 
         tree->Branch("retunedSuperCluster_simOOTPU_nSharedXtals","std::vector<double>",&retunedSuperCluster_simOOTPU_nSharedXtals);
         tree->Branch("retunedSuperCluster_simEnergy_sharedXtalsOOTPU","std::vector<double>",&retunedSuperCluster_simEnergy_sharedXtalsOOTPU);
         tree->Branch("retunedSuperCluster_recoEnergy_sharedXtalsOOTPU","std::vector<double>",&retunedSuperCluster_recoEnergy_sharedXtalsOOTPU);  
         tree->Branch("retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU","std::vector<double>",&retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU);
         tree->Branch("retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU","std::vector<double>",&retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU);  
      }
   }    
   if(saveDeepSC_){
      tree->Branch("deepSuperCluster_seedRawId","std::vector<uint32_t> ",&deepSuperCluster_seedRawId); 
      tree->Branch("deepSuperCluster_rawEnergy","std::vector<float> ",&deepSuperCluster_rawEnergy);
      tree->Branch("deepSuperCluster_rawESEnergy","std::vector<float> ",&deepSuperCluster_rawESEnergy);
      tree->Branch("deepSuperCluster_energy","std::vector<float> ",&deepSuperCluster_energy);
      tree->Branch("deepSuperCluster_eta","std::vector<float>",&deepSuperCluster_eta);
      tree->Branch("deepSuperCluster_phi","std::vector<float>",&deepSuperCluster_phi);  
      tree->Branch("deepSuperCluster_etaWidth","std::vector<float>",&deepSuperCluster_etaWidth);
      tree->Branch("deepSuperCluster_phiWidth","std::vector<float>",&deepSuperCluster_phiWidth);  
      tree->Branch("deepSuperCluster_R","std::vector<float>",&deepSuperCluster_R);   
      tree->Branch("deepSuperCluster_nPFClusters","std::vector<int>",&deepSuperCluster_nPFClusters);   
      tree->Branch("deepSuperCluster_ieta","std::vector<int>",&deepSuperCluster_ieta);
      tree->Branch("deepSuperCluster_iphi","std::vector<int>",&deepSuperCluster_iphi);  
      tree->Branch("deepSuperCluster_iz","std::vector<int>",&deepSuperCluster_iz);   
      tree->Branch("deepSuperCluster_nXtals","std::vector<int>",&deepSuperCluster_nXtals);    
      if(savePFCluster_) tree->Branch("deepSuperCluster_seedIndex","std::vector<int>",&deepSuperCluster_seedIndex);     
      if(savePFCluster_) tree->Branch("deepSuperCluster_pfClustersIndex","std::vector<std::vector<int> >",&deepSuperCluster_pfClustersIndex); 
      tree->Branch("deepSuperCluster_psCluster_energy", "std::vector<std::vector<float> >", &deepSuperCluster_psCluster_energy);
      tree->Branch("deepSuperCluster_psCluster_eta", "std::vector<std::vector<float> >", &deepSuperCluster_psCluster_eta);
      tree->Branch("deepSuperCluster_psCluster_phi", "std::vector<std::vector<float> >", &deepSuperCluster_psCluster_phi);
      if(saveCaloParticles_ && isMC_){
         tree->Branch("deepSuperCluster_dR_genScore","std::vector<std::vector<double> >",&deepSuperCluster_dR_genScore);
         tree->Branch("deepSuperCluster_dR_simScore","std::vector<std::vector<double> >",&deepSuperCluster_dR_simScore);
         tree->Branch("deepSuperCluster_sim_nSharedXtals","std::vector<std::vector<double> >",&deepSuperCluster_sim_nSharedXtals);
         tree->Branch("deepSuperCluster_sim_fraction_noHitsFraction","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction_noHitsFraction); 
         tree->Branch("deepSuperCluster_sim_fraction","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction);
         tree->Branch("deepSuperCluster_recoToSim_fraction","std::vector<std::vector<double> >",&deepSuperCluster_recoToSim_fraction);
         tree->Branch("deepSuperCluster_recoToSim_fraction_sharedXtals","std::vector<std::vector<double> >",&deepSuperCluster_recoToSim_fraction_sharedXtals);
         tree->Branch("deepSuperCluster_simEnergy_sharedXtals","std::vector<std::vector<double> >",&deepSuperCluster_simEnergy_sharedXtals);
         tree->Branch("deepSuperCluster_recoEnergy_sharedXtals","std::vector<std::vector<double> >",&deepSuperCluster_recoEnergy_sharedXtals);  
         tree->Branch("deepSuperCluster_dR_genScore_MatchedIndex","std::vector<int>",&deepSuperCluster_dR_genScore_MatchedIndex);
         tree->Branch("deepSuperCluster_dR_simScore_MatchedIndex","std::vector<int>",&deepSuperCluster_dR_simScore_MatchedIndex);
         tree->Branch("deepSuperCluster_sim_nSharedXtals_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_nSharedXtals_MatchedIndex);
         tree->Branch("deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<int>", &deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex); 
         tree->Branch("deepSuperCluster_sim_fraction_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_fraction_MatchedIndex);
         tree->Branch("deepSuperCluster_recoToSim_fraction_MatchedIndex","std::vector<int>",&deepSuperCluster_recoToSim_fraction_MatchedIndex);
         tree->Branch("deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<int>", &deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);
         tree->Branch("deepSuperCluster_simEnergy_sharedXtals_MatchedIndex","std::vector<int>",&deepSuperCluster_simEnergy_sharedXtals_MatchedIndex);
         tree->Branch("deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex","std::vector<int>",&deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex);
      }   
      if(saveCaloParticlesPU_ && isMC_){ 
         tree->Branch("deepSuperCluster_simPU_nSharedXtals","std::vector<double>",&deepSuperCluster_simPU_nSharedXtals);
         tree->Branch("deepSuperCluster_simEnergy_sharedXtalsPU","std::vector<double>",&deepSuperCluster_simEnergy_sharedXtalsPU);
         tree->Branch("deepSuperCluster_recoEnergy_sharedXtalsPU","std::vector<double>",&deepSuperCluster_recoEnergy_sharedXtalsPU);  
         tree->Branch("deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU","std::vector<double>",&deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU);
         tree->Branch("deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU","std::vector<double>",&deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU);  
      }
      if(saveCaloParticlesOOTPU_ && isMC_){ 
         tree->Branch("deepSuperCluster_simOOTPU_nSharedXtals","std::vector<double>",&deepSuperCluster_simOOTPU_nSharedXtals);
         tree->Branch("deepSuperCluster_simEnergy_sharedXtalsOOTPU","std::vector<double>",&deepSuperCluster_simEnergy_sharedXtalsOOTPU);
         tree->Branch("deepSuperCluster_recoEnergy_sharedXtalsOOTPU","std::vector<double>",&deepSuperCluster_recoEnergy_sharedXtalsOOTPU); 
         tree->Branch("deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU","std::vector<double>",&deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU);
         tree->Branch("deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU","std::vector<double>",&deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU);  
      }   
   }
   if(saveSuperCluster_ && saveShowerShapes_){  
      tree->Branch("superCluster_e5x5","std::vector<float>",&superCluster_e5x5);
      tree->Branch("superCluster_e2x2Ratio","std::vector<float>",&superCluster_e2x2Ratio);
      tree->Branch("superCluster_e3x3Ratio","std::vector<float>",&superCluster_e3x3Ratio);
      tree->Branch("superCluster_eMaxRatio","std::vector<float>",&superCluster_eMaxRatio);
      tree->Branch("superCluster_e2ndRatio","std::vector<float>",&superCluster_e2ndRatio);
      tree->Branch("superCluster_eTopRatio","std::vector<float>",&superCluster_eTopRatio);
      tree->Branch("superCluster_eRightRatio","std::vector<float>",&superCluster_eRightRatio);
      tree->Branch("superCluster_eBottomRatio","std::vector<float>",&superCluster_eBottomRatio);
      tree->Branch("superCluster_eLeftRatio","std::vector<float>",&superCluster_eLeftRatio);
      tree->Branch("superCluster_e2x5MaxRatio","std::vector<float>",&superCluster_e2x5MaxRatio);
      tree->Branch("superCluster_e2x5TopRatio","std::vector<float>",&superCluster_e2x5TopRatio);
      tree->Branch("superCluster_e2x5RightRatio","std::vector<float>",&superCluster_e2x5RightRatio);
      tree->Branch("superCluster_e2x5BottomRatio","std::vector<float>",&superCluster_e2x5BottomRatio); 
      tree->Branch("superCluster_e2x5LeftRatio","std::vector<float>",&superCluster_e2x5LeftRatio); 
      tree->Branch("superCluster_swissCross","std::vector<float>",&superCluster_swissCross); 
      tree->Branch("superCluster_r9","std::vector<float>",&superCluster_r9);
      tree->Branch("superCluster_sigmaIetaIeta","std::vector<float>",&superCluster_sigmaIetaIeta);
      tree->Branch("superCluster_sigmaIetaIphi","std::vector<float>",&superCluster_sigmaIetaIphi);
      tree->Branch("superCluster_sigmaIphiIphi","std::vector<float>",&superCluster_sigmaIphiIphi);
      tree->Branch("superCluster_full5x5_e5x5","std::vector<float>",&superCluster_full5x5_e5x5);
      tree->Branch("superCluster_full5x5_e2x2Ratio","std::vector<float>",&superCluster_full5x5_e2x2Ratio);
      tree->Branch("superCluster_full5x5_e3x3Ratio","std::vector<float>",&superCluster_full5x5_e3x3Ratio);
      tree->Branch("superCluster_full5x5_eMaxRatio","std::vector<float>",&superCluster_full5x5_eMaxRatio);
      tree->Branch("superCluster_full5x5_e2ndRatio","std::vector<float>",&superCluster_full5x5_e2ndRatio);
      tree->Branch("superCluster_full5x5_eTopRatio","std::vector<float>",&superCluster_full5x5_eTopRatio);
      tree->Branch("superCluster_full5x5_eRightRatio","std::vector<float>",&superCluster_full5x5_eRightRatio);
      tree->Branch("superCluster_full5x5_eBottomRatio","std::vector<float>",&superCluster_full5x5_eBottomRatio);
      tree->Branch("superCluster_full5x5_eLeftRatio","std::vector<float>",&superCluster_full5x5_eLeftRatio);
      tree->Branch("superCluster_full5x5_e2x5MaxRatio","std::vector<float>",&superCluster_full5x5_e2x5MaxRatio);
      tree->Branch("superCluster_full5x5_e2x5TopRatio","std::vector<float>",&superCluster_full5x5_e2x5TopRatio);
      tree->Branch("superCluster_full5x5_e2x5RightRatio","std::vector<float>",&superCluster_full5x5_e2x5RightRatio);
      tree->Branch("superCluster_full5x5_e2x5BottomRatio","std::vector<float>",&superCluster_full5x5_e2x5BottomRatio); 
      tree->Branch("superCluster_full5x5_e2x5LeftRatio","std::vector<float>",&superCluster_full5x5_e2x5LeftRatio); 
      tree->Branch("superCluster_full5x5_swissCross","std::vector<float>",&superCluster_full5x5_swissCross); 
      tree->Branch("superCluster_full5x5_r9","std::vector<float>",&superCluster_full5x5_r9);
      tree->Branch("superCluster_full5x5_sigmaIetaIeta","std::vector<float>",&superCluster_full5x5_sigmaIetaIeta);
      tree->Branch("superCluster_full5x5_sigmaIetaIphi","std::vector<float>",&superCluster_full5x5_sigmaIetaIphi);
      tree->Branch("superCluster_full5x5_sigmaIphiIphi","std::vector<float>",&superCluster_full5x5_sigmaIphiIphi);
      if(saveRetunedSC_){
         tree->Branch("retunedSuperCluster_e5x5","std::vector<float>",&retunedSuperCluster_e5x5);
         tree->Branch("retunedSuperCluster_e2x2Ratio","std::vector<float>",&retunedSuperCluster_e2x2Ratio);
         tree->Branch("retunedSuperCluster_e3x3Ratio","std::vector<float>",&retunedSuperCluster_e3x3Ratio);
         tree->Branch("retunedSuperCluster_eMaxRatio","std::vector<float>",&retunedSuperCluster_eMaxRatio);
         tree->Branch("retunedSuperCluster_e2ndRatio","std::vector<float>",&retunedSuperCluster_e2ndRatio);
         tree->Branch("retunedSuperCluster_eTopRatio","std::vector<float>",&retunedSuperCluster_eTopRatio);
         tree->Branch("retunedSuperCluster_eRightRatio","std::vector<float>",&retunedSuperCluster_eRightRatio);
         tree->Branch("retunedSuperCluster_eBottomRatio","std::vector<float>",&retunedSuperCluster_eBottomRatio);
         tree->Branch("retunedSuperCluster_eLeftRatio","std::vector<float>",&retunedSuperCluster_eLeftRatio);
         tree->Branch("retunedSuperCluster_e2x5MaxRatio","std::vector<float>",&retunedSuperCluster_e2x5MaxRatio);
         tree->Branch("retunedSuperCluster_e2x5TopRatio","std::vector<float>",&retunedSuperCluster_e2x5TopRatio);
         tree->Branch("retunedSuperCluster_e2x5RightRatio","std::vector<float>",&retunedSuperCluster_e2x5RightRatio);
         tree->Branch("retunedSuperCluster_e2x5BottomRatio","std::vector<float>",&retunedSuperCluster_e2x5BottomRatio); 
         tree->Branch("retunedSuperCluster_e2x5LeftRatio","std::vector<float>",&retunedSuperCluster_e2x5LeftRatio); 
         tree->Branch("retunedSuperCluster_swissCross","std::vector<float>",&retunedSuperCluster_swissCross); 
         tree->Branch("retunedSuperCluster_r9","std::vector<float>",&retunedSuperCluster_r9);
         tree->Branch("retunedSuperCluster_sigmaIetaIeta","std::vector<float>",&retunedSuperCluster_sigmaIetaIeta);
         tree->Branch("retunedSuperCluster_sigmaIetaIphi","std::vector<float>",&retunedSuperCluster_sigmaIetaIphi);
         tree->Branch("retunedSuperCluster_sigmaIphiIphi","std::vector<float>",&retunedSuperCluster_sigmaIphiIphi);
         tree->Branch("retunedSuperCluster_full5x5_e5x5","std::vector<float>",&retunedSuperCluster_full5x5_e5x5);
         tree->Branch("retunedSuperCluster_full5x5_e2x2Ratio","std::vector<float>",&retunedSuperCluster_full5x5_e2x2Ratio);
         tree->Branch("retunedSuperCluster_full5x5_e3x3Ratio","std::vector<float>",&retunedSuperCluster_full5x5_e3x3Ratio);
         tree->Branch("retunedSuperCluster_full5x5_eMaxRatio","std::vector<float>",&retunedSuperCluster_full5x5_eMaxRatio);
         tree->Branch("retunedSuperCluster_full5x5_e2ndRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2ndRatio);
         tree->Branch("retunedSuperCluster_full5x5_eTopRatio","std::vector<float>",&retunedSuperCluster_full5x5_eTopRatio);
         tree->Branch("retunedSuperCluster_full5x5_eRightRatio","std::vector<float>",&retunedSuperCluster_full5x5_eRightRatio);
         tree->Branch("retunedSuperCluster_full5x5_eBottomRatio","std::vector<float>",&retunedSuperCluster_full5x5_eBottomRatio);
         tree->Branch("retunedSuperCluster_full5x5_eLeftRatio","std::vector<float>",&retunedSuperCluster_full5x5_eLeftRatio);
         tree->Branch("retunedSuperCluster_full5x5_e2x5MaxRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2x5MaxRatio);
         tree->Branch("retunedSuperCluster_full5x5_e2x5TopRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2x5TopRatio);
         tree->Branch("retunedSuperCluster_full5x5_e2x5RightRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2x5RightRatio);
         tree->Branch("retunedSuperCluster_full5x5_e2x5BottomRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2x5BottomRatio); 
         tree->Branch("retunedSuperCluster_full5x5_e2x5LeftRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2x5LeftRatio); 
         tree->Branch("retunedSuperCluster_full5x5_swissCross","std::vector<float>",&retunedSuperCluster_full5x5_swissCross); 
         tree->Branch("retunedSuperCluster_full5x5_r9","std::vector<float>",&retunedSuperCluster_full5x5_r9);
         tree->Branch("retunedSuperCluster_full5x5_sigmaIetaIeta","std::vector<float>",&retunedSuperCluster_full5x5_sigmaIetaIeta);
         tree->Branch("retunedSuperCluster_full5x5_sigmaIetaIphi","std::vector<float>",&retunedSuperCluster_full5x5_sigmaIetaIphi);
         tree->Branch("retunedSuperCluster_full5x5_sigmaIphiIphi","std::vector<float>",&retunedSuperCluster_full5x5_sigmaIphiIphi);
      }
      if(saveDeepSC_){
         tree->Branch("deepSuperCluster_e5x5","std::vector<float>",&deepSuperCluster_e5x5);
         tree->Branch("deepSuperCluster_e2x2Ratio","std::vector<float>",&deepSuperCluster_e2x2Ratio);
         tree->Branch("deepSuperCluster_e3x3Ratio","std::vector<float>",&deepSuperCluster_e3x3Ratio);
         tree->Branch("deepSuperCluster_eMaxRatio","std::vector<float>",&deepSuperCluster_eMaxRatio);
         tree->Branch("deepSuperCluster_e2ndRatio","std::vector<float>",&deepSuperCluster_e2ndRatio);
         tree->Branch("deepSuperCluster_eTopRatio","std::vector<float>",&deepSuperCluster_eTopRatio);
         tree->Branch("deepSuperCluster_eRightRatio","std::vector<float>",&deepSuperCluster_eRightRatio);
         tree->Branch("deepSuperCluster_eBottomRatio","std::vector<float>",&deepSuperCluster_eBottomRatio);
         tree->Branch("deepSuperCluster_eLeftRatio","std::vector<float>",&deepSuperCluster_eLeftRatio);
         tree->Branch("deepSuperCluster_e2x5MaxRatio","std::vector<float>",&deepSuperCluster_e2x5MaxRatio);
         tree->Branch("deepSuperCluster_e2x5TopRatio","std::vector<float>",&deepSuperCluster_e2x5TopRatio);
         tree->Branch("deepSuperCluster_e2x5RightRatio","std::vector<float>",&deepSuperCluster_e2x5RightRatio);
         tree->Branch("deepSuperCluster_e2x5BottomRatio","std::vector<float>",&deepSuperCluster_e2x5BottomRatio); 
         tree->Branch("deepSuperCluster_e2x5LeftRatio","std::vector<float>",&deepSuperCluster_e2x5LeftRatio); 
         tree->Branch("deepSuperCluster_swissCross","std::vector<float>",&deepSuperCluster_swissCross); 
         tree->Branch("deepSuperCluster_r9","std::vector<float>",&deepSuperCluster_r9);
         tree->Branch("deepSuperCluster_sigmaIetaIeta","std::vector<float>",&deepSuperCluster_sigmaIetaIeta);
         tree->Branch("deepSuperCluster_sigmaIetaIphi","std::vector<float>",&deepSuperCluster_sigmaIetaIphi);
         tree->Branch("deepSuperCluster_sigmaIphiIphi","std::vector<float>",&deepSuperCluster_sigmaIphiIphi);
         tree->Branch("deepSuperCluster_full5x5_e5x5","std::vector<float>",&deepSuperCluster_full5x5_e5x5);
         tree->Branch("deepSuperCluster_full5x5_e2x2Ratio","std::vector<float>",&deepSuperCluster_full5x5_e2x2Ratio);
         tree->Branch("deepSuperCluster_full5x5_e3x3Ratio","std::vector<float>",&deepSuperCluster_full5x5_e3x3Ratio);
         tree->Branch("deepSuperCluster_full5x5_eMaxRatio","std::vector<float>",&deepSuperCluster_full5x5_eMaxRatio);
         tree->Branch("deepSuperCluster_full5x5_e2ndRatio","std::vector<float>",&deepSuperCluster_full5x5_e2ndRatio);
         tree->Branch("deepSuperCluster_full5x5_eTopRatio","std::vector<float>",&deepSuperCluster_full5x5_eTopRatio);
         tree->Branch("deepSuperCluster_full5x5_eRightRatio","std::vector<float>",&deepSuperCluster_full5x5_eRightRatio);
         tree->Branch("deepSuperCluster_full5x5_eBottomRatio","std::vector<float>",&deepSuperCluster_full5x5_eBottomRatio);
         tree->Branch("deepSuperCluster_full5x5_eLeftRatio","std::vector<float>",&deepSuperCluster_full5x5_eLeftRatio);
         tree->Branch("deepSuperCluster_full5x5_e2x5MaxRatio","std::vector<float>",&deepSuperCluster_full5x5_e2x5MaxRatio);
         tree->Branch("deepSuperCluster_full5x5_e2x5TopRatio","std::vector<float>",&deepSuperCluster_full5x5_e2x5TopRatio);
         tree->Branch("deepSuperCluster_full5x5_e2x5RightRatio","std::vector<float>",&deepSuperCluster_full5x5_e2x5RightRatio);
         tree->Branch("deepSuperCluster_full5x5_e2x5BottomRatio","std::vector<float>",&deepSuperCluster_full5x5_e2x5BottomRatio); 
         tree->Branch("deepSuperCluster_full5x5_e2x5LeftRatio","std::vector<float>",&deepSuperCluster_full5x5_e2x5LeftRatio); 
         tree->Branch("deepSuperCluster_full5x5_swissCross","std::vector<float>",&deepSuperCluster_full5x5_swissCross); 
         tree->Branch("deepSuperCluster_full5x5_r9","std::vector<float>",&deepSuperCluster_full5x5_r9);
         tree->Branch("deepSuperCluster_full5x5_sigmaIetaIeta","std::vector<float>",&deepSuperCluster_full5x5_sigmaIetaIeta);
         tree->Branch("deepSuperCluster_full5x5_sigmaIetaIphi","std::vector<float>",&deepSuperCluster_full5x5_sigmaIetaIphi);
         tree->Branch("deepSuperCluster_full5x5_sigmaIphiIphi","std::vector<float>",&deepSuperCluster_full5x5_sigmaIphiIphi); 
      } 
   }
}

void RecoSimDumper::setVectors(int nGenParticles, int nCaloParticles, int nPFClusters, int nSuperClustersEB, int nSuperClustersEE, int nRetunedSuperClustersEB, int nRetunedSuperClustersEE, int nDeepSuperClustersEB, int nDeepSuperClustersEE)
{
   genParticle_pfCluster_dR_genScore_MatchedIndex.clear();
   genParticle_superCluster_dR_genScore_MatchedIndex.clear();
   genParticle_retunedSuperCluster_dR_genScore_MatchedIndex.clear();
   genParticle_deepSuperCluster_dR_genScore_MatchedIndex.clear();
   if(saveGenParticles_){
      genParticle_pfCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
      genParticle_superCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
      genParticle_retunedSuperCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
      genParticle_deepSuperCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
   }

   caloParticle_index.clear();
   caloParticle_nXtals.clear();
   caloParticle_pdgId.clear(); 
   caloParticle_status.clear(); 
   caloParticle_charge.clear(); 
   caloParticle_genEnergy.clear(); 
   caloParticle_simEnergy.clear(); 
   caloParticle_simEnergyGoodStatus.clear();  
   caloParticle_simEnergyWithES.clear(); 
   caloParticle_genPt.clear(); 
   caloParticle_simPt.clear(); 
   caloParticle_genEta.clear(); 
   caloParticle_simEta.clear(); 
   caloParticle_genPhi.clear(); 
   caloParticle_simPhi.clear(); 
   caloParticle_partonIndex.clear(); 
   caloParticle_partonPdgId.clear();  
   caloParticle_partonCharge.clear();  
   caloParticle_partonEnergy.clear(); 
   caloParticle_partonPt.clear(); 
   caloParticle_partonEta.clear(); 
   caloParticle_partonPhi.clear();  
   caloParticle_genMotherPdgId.clear();
   caloParticle_genMotherStatus.clear();  
   caloParticle_genMotherCharge.clear();  
   caloParticle_genMotherEnergy.clear(); 
   caloParticle_genMotherPt.clear(); 
   caloParticle_genMotherEta.clear(); 
   caloParticle_genMotherPhi.clear();  
   caloParticle_simIeta.clear(); 
   caloParticle_simIphi.clear(); 
   caloParticle_simIz.clear(); 
   caloParticle_nSharedXtals.clear();
   caloParticle_sharedIndex1.clear();
   caloParticle_sharedIndex2.clear();
   caloParticle_sharedEnergyFrac1.clear();
   caloParticle_sharedEnergyFrac2.clear();
   caloParticle_pfCluster_dR_simScore_MatchedIndex.clear(); 
   caloParticle_pfCluster_sim_nSharedXtals_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_MatchedIndex.clear(); 
   caloParticle_pfCluster_recoToSim_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   caloParticle_pfCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   caloParticle_pfCluster_recoEnergy_sharedXtals_MatchedIndex.clear();     
   caloParticle_superCluster_dR_simScore_MatchedIndex.clear(); 
   caloParticle_superCluster_sim_nSharedXtals_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_MatchedIndex.clear(); 
   caloParticle_superCluster_recoToSim_fraction_MatchedIndex.clear();
   caloParticle_superCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   caloParticle_superCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   caloParticle_superCluster_recoEnergy_sharedXtals_MatchedIndex.clear();  
   caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex.clear(); 
   caloParticle_retunedSuperCluster_sim_nSharedXtals_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex.clear(); 
   caloParticle_retunedSuperCluster_recoToSim_fraction_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   caloParticle_retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   caloParticle_retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex.clear();       
   caloParticle_deepSuperCluster_dR_simScore_MatchedIndex.clear(); 
   caloParticle_deepSuperCluster_sim_nSharedXtals_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.clear(); 
   caloParticle_deepSuperCluster_recoToSim_fraction_MatchedIndex.clear();
   caloParticle_deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   caloParticle_deepSuperCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   caloParticle_deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex.clear(); 
   if(saveCaloParticles_){ 
      caloParticle_pfCluster_dR_simScore_MatchedIndex.resize(nCaloParticles); 
      caloParticle_pfCluster_sim_nSharedXtals_MatchedIndex.resize(nCaloParticles);
      caloParticle_pfCluster_sim_fraction_noHitsFraction_MatchedIndex.resize(nCaloParticles);
      caloParticle_pfCluster_sim_fraction_MatchedIndex.resize(nCaloParticles); 
      caloParticle_pfCluster_recoToSim_fraction_MatchedIndex.resize(nCaloParticles);
      caloParticle_pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex.resize(nCaloParticles);  
      caloParticle_pfCluster_simEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles); 
      caloParticle_pfCluster_recoEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles);     
      caloParticle_superCluster_dR_simScore_MatchedIndex.resize(nCaloParticles); 
      caloParticle_superCluster_sim_nSharedXtals_MatchedIndex.resize(nCaloParticles);
      caloParticle_superCluster_sim_fraction_noHitsFraction_MatchedIndex.resize(nCaloParticles);
      caloParticle_superCluster_sim_fraction_MatchedIndex.resize(nCaloParticles); 
      caloParticle_superCluster_recoToSim_fraction_MatchedIndex.resize(nCaloParticles);
      caloParticle_superCluster_recoToSim_fraction_sharedXtals_MatchedIndex.resize(nCaloParticles);  
      caloParticle_superCluster_simEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles); 
      caloParticle_superCluster_recoEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles);  
      caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex.resize(nCaloParticles); 
      caloParticle_retunedSuperCluster_sim_nSharedXtals_MatchedIndex.resize(nCaloParticles);
      caloParticle_retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.resize(nCaloParticles);
      caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex.resize(nCaloParticles); 
      caloParticle_retunedSuperCluster_recoToSim_fraction_MatchedIndex.resize(nCaloParticles);
      caloParticle_retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.resize(nCaloParticles);  
      caloParticle_retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles); 
      caloParticle_retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles);       
      caloParticle_deepSuperCluster_dR_simScore_MatchedIndex.resize(nCaloParticles); 
      caloParticle_deepSuperCluster_sim_nSharedXtals_MatchedIndex.resize(nCaloParticles);
      caloParticle_deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.resize(nCaloParticles);
      caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.resize(nCaloParticles); 
      caloParticle_deepSuperCluster_recoToSim_fraction_MatchedIndex.resize(nCaloParticles);
      caloParticle_deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.resize(nCaloParticles);  
      caloParticle_deepSuperCluster_simEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles); 
      caloParticle_deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles);  
   }
 
   simHit_energy.clear();
   simHit_eta.clear();
   simHit_phi.clear();
   simHit_ieta.clear();
   simHit_iphi.clear();
   simHit_iz.clear();
   simHit_iplane.clear();
   simHit_chStatus.clear();
   if(saveSimhits_ && saveCaloParticles_){
      simHit_energy.resize(nCaloParticles);
      simHit_eta.resize(nCaloParticles);
      simHit_phi.resize(nCaloParticles);
      simHit_ieta.resize(nCaloParticles);
      simHit_iphi.resize(nCaloParticles);
      simHit_iz.resize(nCaloParticles);
      simHit_iplane.resize(nCaloParticles);
      simHit_chStatus.resize(nCaloParticles);
   }

   recHit_noPF_energy.clear();
   recHit_noPF_eta.clear();
   recHit_noPF_phi.clear();
   recHit_noPF_ieta.clear();
   recHit_noPF_iphi.clear();
   recHit_noPF_iz.clear();  

   pfRecHit_unClustered_energy.clear();
   pfRecHit_unClustered_eta.clear();
   pfRecHit_unClustered_phi.clear();
   pfRecHit_unClustered_ieta.clear();
   pfRecHit_unClustered_iphi.clear();
   pfRecHit_unClustered_iz.clear();

   pfCluster_rawEnergy.clear();
   pfCluster_rawEnergyUncalib.clear();
   pfCluster_energy.clear();
   pfCluster_rawPt.clear();
   pfCluster_pt.clear();
   pfCluster_eta.clear();
   pfCluster_phi.clear();
   pfCluster_ieta.clear();
   pfCluster_iphi.clear();
   pfCluster_iz.clear();
   pfCluster_noise.clear();
   pfCluster_noiseUncalib.clear();
   pfCluster_noiseNoFractions.clear();
   pfCluster_noiseUncalibNoFractions.clear();
   pfCluster_noiseDB.clear();
   pfCluster_noiseDBUncalib.clear();
   pfCluster_noiseDBNoFractions.clear();
   pfCluster_noiseDBUncalibNoFractions.clear();
   pfCluster_superClustersIndex.clear();
   pfCluster_retunedSuperClustersIndex.clear();
   pfCluster_deepSuperClustersIndex.clear(); 
   pfCluster_etaWidth.clear();
   pfCluster_phiWidth.clear();
   pfCluster_e5x5.clear();
   pfCluster_e2x2Ratio.clear();
   pfCluster_e3x3Ratio.clear();
   pfCluster_eMaxRatio.clear();
   pfCluster_e2ndRatio.clear();
   pfCluster_eTopRatio.clear();
   pfCluster_eRightRatio.clear();
   pfCluster_eBottomRatio.clear();
   pfCluster_eLeftRatio.clear();
   pfCluster_e2x5MaxRatio.clear();
   pfCluster_e2x5TopRatio.clear();
   pfCluster_e2x5RightRatio.clear();
   pfCluster_e2x5BottomRatio.clear();
   pfCluster_e2x5LeftRatio.clear();
   pfCluster_swissCross.clear();
   pfCluster_r9.clear();
   pfCluster_sigmaIetaIeta.clear(); 
   pfCluster_sigmaIetaIphi.clear(); 
   pfCluster_sigmaIphiIphi.clear(); 
   pfCluster_full5x5_e5x5.clear();
   pfCluster_full5x5_e2x2Ratio.clear();
   pfCluster_full5x5_e3x3Ratio.clear();
   pfCluster_full5x5_eMaxRatio.clear();
   pfCluster_full5x5_e2ndRatio.clear();
   pfCluster_full5x5_eTopRatio.clear();
   pfCluster_full5x5_eRightRatio.clear();
   pfCluster_full5x5_eBottomRatio.clear();
   pfCluster_full5x5_eLeftRatio.clear();
   pfCluster_full5x5_e2x5MaxRatio.clear();
   pfCluster_full5x5_e2x5TopRatio.clear();
   pfCluster_full5x5_e2x5RightRatio.clear();
   pfCluster_full5x5_e2x5BottomRatio.clear();
   pfCluster_full5x5_e2x5LeftRatio.clear();
   pfCluster_full5x5_swissCross.clear();
   pfCluster_full5x5_r9.clear();
   pfCluster_full5x5_sigmaIetaIeta.clear(); 
   pfCluster_full5x5_sigmaIetaIphi.clear(); 
   pfCluster_full5x5_sigmaIphiIphi.clear();    
   pfCluster_nXtals.clear();
   pfCluster_dR_genScore_MatchedIndex.clear();
   pfCluster_dR_simScore_MatchedIndex.clear();
   pfCluster_sim_nSharedXtals_MatchedIndex.clear();
   pfCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   pfCluster_sim_fraction_MatchedIndex.clear(); 
   pfCluster_recoToSim_fraction_MatchedIndex.clear();
   pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   pfCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   pfCluster_recoEnergy_sharedXtals_MatchedIndex.clear();     
   pfCluster_dR_genScore.resize(nPFClusters);
   pfCluster_dR_simScore.resize(nPFClusters);
   pfCluster_sim_nSharedXtals.resize(nPFClusters);
   pfCluster_sim_fraction_noHitsFraction.resize(nPFClusters);
   pfCluster_sim_fraction.resize(nPFClusters);
   pfCluster_recoToSim_fraction.resize(nPFClusters);
   pfCluster_recoToSim_fraction_sharedXtals.resize(nPFClusters);
   pfCluster_simEnergy_sharedXtals.resize(nPFClusters);
   pfCluster_recoEnergy_sharedXtals.resize(nPFClusters);   
   pfCluster_simPU_nSharedXtals.resize(nPFClusters);
   pfCluster_simEnergy_sharedXtalsPU.resize(nPFClusters);
   pfCluster_recoEnergy_sharedXtalsPU.resize(nPFClusters);   
   pfCluster_simEnergy_noHitsFraction_sharedXtalsPU.resize(nPFClusters);
   pfCluster_recoEnergy_noHitsFraction_sharedXtalsPU.resize(nPFClusters); 
   pfCluster_simOOTPU_nSharedXtals.resize(nPFClusters);
   pfCluster_simEnergy_sharedXtalsOOTPU.resize(nPFClusters);
   pfCluster_recoEnergy_sharedXtalsOOTPU.resize(nPFClusters);   
   pfCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU.resize(nPFClusters);
   pfCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU.resize(nPFClusters); 
   pfCluster_superClustersIndex.resize(nPFClusters); 
   pfCluster_retunedSuperClustersIndex.resize(nPFClusters);  
   pfCluster_deepSuperClustersIndex.resize(nPFClusters); 

   pfClusterHit_fraction.clear();
   pfClusterHit_rechitEnergy.clear();
   pfClusterHit_eta.clear();
   pfClusterHit_phi.clear();   
   pfClusterHit_ieta.clear();
   pfClusterHit_iphi.clear(); 
   pfClusterHit_iz.clear();     
   pfClusterHit_adcToGeV.clear();
   pfClusterHit_ic.clear();   
   pfClusterHit_icMC.clear();
   pfClusterHit_laserCorr.clear(); 
   pfClusterHit_chStatus.clear();           
   pfClusterHit_fraction.resize(nPFClusters);  
   pfClusterHit_rechitEnergy.resize(nPFClusters);  
   pfClusterHit_eta.resize(nPFClusters);  
   pfClusterHit_phi.resize(nPFClusters);     
   pfClusterHit_ieta.resize(nPFClusters);  
   pfClusterHit_iphi.resize(nPFClusters);   
   pfClusterHit_iz.resize(nPFClusters);    
   pfClusterHit_adcToGeV.resize(nPFClusters);
   pfClusterHit_ic.resize(nPFClusters);
   pfClusterHit_icMC.resize(nPFClusters);    
   pfClusterHit_laserCorr.resize(nPFClusters);   
   pfClusterHit_chStatus.resize(nPFClusters);   

   gsfElectron_index.clear();
   gsfElectron_seedRawId.clear();
   gsfElectron_isEB.clear();
   gsfElectron_isEE.clear();
   gsfElectron_isEBEEGap.clear();
   gsfElectron_isEBEtaGap.clear();
   gsfElectron_isEBPhiGap.clear();
   gsfElectron_isEEDeeGap.clear();
   gsfElectron_isEERingGap.clear();  
   gsfElectron_isEcalDriven.clear();
   gsfElectron_isTrackerDriven.clear(); 
   gsfElectron_scNPFClusters.clear();  
   gsfElectron_ecalSCNPFClusters.clear();  
   gsfElectron_classification.clear(); 
   gsfElectron_p.clear();
   gsfElectron_pt.clear();
   gsfElectron_et.clear();
   gsfElectron_energy.clear();
   gsfElectron_energyErr.clear();
   gsfElectron_ecalEnergy.clear();
   gsfElectron_ecalEnergyErr.clear();
   gsfElectron_eta.clear();
   gsfElectron_phi.clear();
   gsfElectron_trkEtaMode.clear();
   gsfElectron_trkPhiMode.clear();
   gsfElectron_trkPMode.clear();
   gsfElectron_trkPModeErr.clear();
   gsfElectron_trkPInn.clear();
   gsfElectron_trkPtInn.clear();
   gsfElectron_trkPVtx.clear();
   gsfElectron_trkPOut.clear();
   gsfElectron_trkChi2.clear();
   gsfElectron_trkNDof.clear();
   gsfElectron_trackFbrem.clear();
   gsfElectron_superClusterFbrem.clear();
   gsfElectron_hademTow.clear();
   gsfElectron_hademCone.clear();
   gsfElectron_ecalDrivenSeed.clear();
   gsfElectron_nrSatCrys.clear();
   gsfElectron_scEnergy.clear();
   gsfElectron_scRawEnergy.clear();
   gsfElectron_scRawESEnergy.clear();
   gsfElectron_scEt.clear();
   gsfElectron_scEtaWidth.clear();
   gsfElectron_scPhiWidth.clear(); 
   gsfElectron_scEoP.clear(); 
   gsfElectron_ecalSCEnergy.clear();
   gsfElectron_ecalSCRawEnergy.clear();
   gsfElectron_ecalSCRawESEnergy.clear();
   gsfElectron_ecalSCEt.clear();
   gsfElectron_ecalSCEtaWidth.clear();
   gsfElectron_ecalSCPhiWidth.clear(); 
   gsfElectron_ecalSCEoP.clear(); 
   gsfElectron_ecalSCEta.clear();
   gsfElectron_ecalSCPhi.clear(); 
   gsfElectron_scEta.clear();
   gsfElectron_scPhi.clear();
   gsfElectron_scSwissCross.clear();
   gsfElectron_scEMax.clear();
   gsfElectron_scE2x2.clear(); 
   gsfElectron_scE3x3.clear();
   gsfElectron_scE5x5.clear(); 
   gsfElectron_scR9.clear(); 
   gsfElectron_scSigmaIEtaIEta.clear();
   gsfElectron_scSigmaIEtaIPhi.clear();
   gsfElectron_scSigmaIPhiIPhi.clear(); 
   gsfElectron_full5x5_scEMax.clear();
   gsfElectron_full5x5_scE2x2.clear(); 
   gsfElectron_full5x5_scE3x3.clear();
   gsfElectron_full5x5_scE5x5.clear(); 
   gsfElectron_full5x5_scR9.clear(); 
   gsfElectron_full5x5_scSigmaIEtaIEta.clear();
   gsfElectron_full5x5_scSigmaIEtaIPhi.clear();
   gsfElectron_full5x5_scSigmaIPhiIPhi.clear();
   gsfElectron_HoE.clear(); 
   gsfElectron_trkIso03.clear(); 
   gsfElectron_ecalIso03.clear();  
   gsfElectron_hcalIso03.clear();
   gsfElectron_trkIso04.clear(); 
   gsfElectron_ecalIso04.clear();  
   gsfElectron_hcalIso04.clear();
   gsfElectron_pfPhotonIso.clear();     
   gsfElectron_pfChargedHadronIso.clear();     
   gsfElectron_pfNeutralHadronIso.clear(); 
   gsfElectron_mva_Isolated.clear(); 
   gsfElectron_mva_e_pi.clear();
   gsfElectron_dnn_signal_Isolated.clear();
   gsfElectron_dnn_signal_nonIsolated.clear();  
   gsfElectron_dnn_bkg_nonIsolated.clear(); 
   gsfElectron_dnn_bkg_Tau.clear();
   gsfElectron_dnn_bkg_Photon.clear();

   gedPhoton_index.clear();  
   gedPhoton_seedRawId.clear();  
   gedPhoton_isEB.clear();
   gedPhoton_isEE.clear();
   gedPhoton_isEBEEGap.clear();
   gedPhoton_isEBEtaGap.clear();
   gedPhoton_isEBPhiGap.clear();
   gedPhoton_isEEDeeGap.clear();
   gedPhoton_isEERingGap.clear();  
   gedPhoton_scNPFClusters.clear();  
   gedPhoton_ecalSCNPFClusters.clear();   
   gedPhoton_hasConversionTracks.clear();
   gedPhoton_nConversions.clear();
   gedPhoton_nConversionsOneLeg.clear();
   gedPhoton_et.clear();  
   gedPhoton_energy.clear();
   gedPhoton_energyErr.clear();
   gedPhoton_ecalEnergy.clear();
   gedPhoton_ecalEnergyErr.clear();
   gedPhoton_eta.clear();  
   gedPhoton_phi.clear();  
   gedPhoton_hademTow.clear();  
   gedPhoton_hademCone.clear();  
   gedPhoton_nrSatCrys.clear();  
   gedPhoton_scEnergy.clear();  
   gedPhoton_scRawEnergy.clear();  
   gedPhoton_scRawESEnergy.clear();
   gedPhoton_scEt.clear();
   gedPhoton_scEtaWidth.clear();
   gedPhoton_scPhiWidth.clear();  
   gedPhoton_ecalSCEnergy.clear();  
   gedPhoton_ecalSCRawEnergy.clear();  
   gedPhoton_ecalSCRawESEnergy.clear();
   gedPhoton_ecalSCEt.clear();
   gedPhoton_ecalSCEtaWidth.clear();
   gedPhoton_ecalSCPhiWidth.clear(); 
   gedPhoton_ecalSCEta.clear();
   gedPhoton_ecalSCPhi.clear();  
   gedPhoton_scEta.clear();
   gedPhoton_scPhi.clear();     
   gedPhoton_scSwissCross.clear();
   gedPhoton_scEMax.clear();
   gedPhoton_scE2x2.clear(); 
   gedPhoton_scE3x3.clear();
   gedPhoton_scE5x5.clear(); 
   gedPhoton_scR9.clear(); 
   gedPhoton_scSigmaIEtaIEta.clear();
   gedPhoton_scSigmaIEtaIPhi.clear();
   gedPhoton_scSigmaIPhiIPhi.clear(); 
   gedPhoton_full5x5_scEMax.clear();
   gedPhoton_full5x5_scE2x2.clear(); 
   gedPhoton_full5x5_scE3x3.clear();
   gedPhoton_full5x5_scE5x5.clear(); 
   gedPhoton_full5x5_scR9.clear(); 
   gedPhoton_full5x5_scSigmaIEtaIEta.clear();
   gedPhoton_full5x5_scSigmaIEtaIPhi.clear();
   gedPhoton_full5x5_scSigmaIPhiIPhi.clear();
   gedPhoton_HoE.clear();  
   gedPhoton_trkIso03.clear(); 
   gedPhoton_ecalIso03.clear();  
   gedPhoton_hcalIso03.clear();
   gedPhoton_trkIso04.clear(); 
   gedPhoton_ecalIso04.clear();  
   gedPhoton_hcalIso04.clear();
   gedPhoton_pfPhotonIso.clear();     
   gedPhoton_pfChargedHadronIso.clear();     
   gedPhoton_pfNeutralHadronIso.clear(); 
   gedPhoton_nClusterOutsideMustache.clear(); 
   gedPhoton_etOutsideMustache.clear(); 
   gedPhoton_pfMVA.clear(); 
   gedPhoton_pfDNN.clear(); 

   patElectron_index.clear();
   patElectron_seedRawId.clear();
   patElectron_classification.clear();
   patElectron_scNPFClusters.clear();  
   patElectron_ecalSCNPFClusters.clear();   
   patElectron_charge.clear();
   patElectron_isEB.clear();
   patElectron_isEE.clear();
   patElectron_isEBEEGap.clear();
   patElectron_isEBEtaGap.clear();
   patElectron_isEBPhiGap.clear();
   patElectron_isEEDeeGap.clear();
   patElectron_isEERingGap.clear();  
   patElectron_isEcalDriven.clear();
   patElectron_isTrackerDriven.clear(); 
   patElectron_passConversionVeto.clear(); 
   patElectron_nOverlapPhotons.clear();
   patElectron_overlapPhotonIndices.clear(); 
   patElectron_hasOverlapJet.clear();
   patElectron_overlapJetIndex.clear();
   patElectron_eta.clear();
   patElectron_phi.clear();
   patElectron_p.clear();
   patElectron_pt.clear();
   patElectron_pIn.clear();
   patElectron_pOut.clear();
   patElectron_pAtCalo.clear();
   patElectron_deltaEtaIn.clear();
   patElectron_deltaPhiIn.clear();
   patElectron_deltaEtaSeedClusterAtCalo.clear();
   patElectron_deltaEtaEleClusterAtCalo.clear();
   patElectron_deltaPhiEleClusterAtCalo.clear();
   patElectron_deltaPhiSeedClusterAtCalo.clear();
   patElectron_misHits.clear();
   patElectron_nAmbiguousGsfTracks.clear(); 
   patElectron_trackFbrem.clear();
   patElectron_superClusterFbrem.clear();
   patElectron_dz.clear();
   patElectron_dzError.clear();
   patElectron_dxy.clear();
   patElectron_dxyError.clear();
   patElectron_energy.clear();
   patElectron_energyErr.clear();
   patElectron_ecalEnergy.clear();
   patElectron_ecalEnergyErr.clear();
   patElectron_et.clear();
   patElectron_mt.clear();
   patElectron_dphiMET.clear();
   patElectron_scEnergy.clear();
   patElectron_scRawEnergy.clear();
   patElectron_scRawESEnergy.clear();
   patElectron_scEt.clear();
   patElectron_scEtaWidth.clear();
   patElectron_scPhiWidth.clear(); 
   patElectron_scEoP.clear(); 
   patElectron_ecalSCEnergy.clear();
   patElectron_ecalSCRawEnergy.clear();
   patElectron_ecalSCRawESEnergy.clear();
   patElectron_ecalSCEt.clear();
   patElectron_ecalSCEtaWidth.clear();
   patElectron_ecalSCPhiWidth.clear(); 
   patElectron_ecalSCEoP.clear(); 
   patElectron_ecalSCEta.clear();
   patElectron_ecalSCPhi.clear();
   patElectron_scEta.clear();
   patElectron_scPhi.clear();
   patElectron_scSwissCross.clear();
   patElectron_scEMax.clear();
   patElectron_scE2x2.clear(); 
   patElectron_scE3x3.clear();
   patElectron_scE5x5.clear(); 
   patElectron_scR9.clear(); 
   patElectron_scSigmaIEtaIEta.clear();
   patElectron_scSigmaIEtaIPhi.clear();
   patElectron_scSigmaIPhiIPhi.clear(); 
   patElectron_full5x5_scEMax.clear();
   patElectron_full5x5_scE2x2.clear(); 
   patElectron_full5x5_scE3x3.clear();
   patElectron_full5x5_scE5x5.clear(); 
   patElectron_full5x5_scR9.clear(); 
   patElectron_full5x5_scSigmaIEtaIEta.clear();
   patElectron_full5x5_scSigmaIEtaIPhi.clear();
   patElectron_full5x5_scSigmaIPhiIPhi.clear();
   patElectron_HoE.clear();
   patElectron_trkIso03.clear(); 
   patElectron_ecalIso03.clear();  
   patElectron_hcalIso03.clear();
   patElectron_trkIso04.clear(); 
   patElectron_ecalIso04.clear();  
   patElectron_hcalIso04.clear(); 
   patElectron_pfPhotonIso.clear();     
   patElectron_pfChargedHadronIso.clear();     
   patElectron_pfNeutralHadronIso.clear(); 
   patElectron_mva_Isolated.clear(); 
   patElectron_mva_e_pi.clear();
   patElectron_dnn_signal_Isolated.clear();
   patElectron_dnn_signal_nonIsolated.clear();  
   patElectron_dnn_bkg_nonIsolated.clear(); 
   patElectron_dnn_bkg_Tau.clear();
   patElectron_dnn_bkg_Photon.clear();
   patElectron_egmCutBasedElectronIDVeto.clear();
   patElectron_egmCutBasedElectronIDloose.clear();
   patElectron_egmCutBasedElectronIDmedium.clear();
   patElectron_egmCutBasedElectronIDtight.clear();
   patElectron_egmMVAElectronIDloose.clear();
   patElectron_egmMVAElectronIDmedium.clear();
   patElectron_egmMVAElectronIDtight.clear();
   patElectron_egmMVAElectronIDlooseNoIso.clear();
   patElectron_egmMVAElectronIDmediumNoIso.clear();
   patElectron_egmMVAElectronIDtightNoIso.clear();
   patElectron_heepElectronID.clear();

   patPhoton_index.clear();
   patPhoton_seedRawId.clear();
   patPhoton_scNPFClusters.clear();  
   patPhoton_ecalSCNPFClusters.clear(); 
   patPhoton_isEB.clear();
   patPhoton_isEE.clear();
   patPhoton_isEBEEGap.clear();
   patPhoton_isEBEtaGap.clear();
   patPhoton_isEBPhiGap.clear();
   patPhoton_isEEDeeGap.clear();
   patPhoton_isEERingGap.clear();  
   patPhoton_passElectronVeto.clear();
   patPhoton_hasPixelSeed.clear();
   patPhoton_hasConversionTracks.clear();
   patPhoton_nConversions.clear();
   patPhoton_nConversionsOneLeg.clear(); 
   patPhoton_hasOverlapElectron.clear();
   patPhoton_overlapElectronIndex.clear();
   patPhoton_hasOverlapJet.clear();
   patPhoton_overlapJetIndex.clear();
   patPhoton_eta.clear();
   patPhoton_phi.clear();
   patPhoton_energy.clear();
   patPhoton_energyErr.clear();
   patPhoton_ecalEnergy.clear();
   patPhoton_ecalEnergyErr.clear();
   patPhoton_et.clear();
   patPhoton_mt.clear();
   patPhoton_dphiMET.clear();  
   patPhoton_scEnergy.clear();  
   patPhoton_scRawEnergy.clear();  
   patPhoton_scRawESEnergy.clear();
   patPhoton_scEt.clear();
   patPhoton_scEtaWidth.clear();
   patPhoton_scPhiWidth.clear(); 
   patPhoton_ecalSCEnergy.clear();  
   patPhoton_ecalSCRawEnergy.clear();  
   patPhoton_ecalSCRawESEnergy.clear();
   patPhoton_ecalSCEt.clear();
   patPhoton_ecalSCEtaWidth.clear();
   patPhoton_ecalSCPhiWidth.clear();
   patPhoton_ecalSCEta.clear();
   patPhoton_ecalSCPhi.clear();  
   patPhoton_scEta.clear();
   patPhoton_scPhi.clear();   
   patPhoton_scSwissCross.clear();
   patPhoton_scEMax.clear();
   patPhoton_scE2x2.clear(); 
   patPhoton_scE3x3.clear();
   patPhoton_scE5x5.clear(); 
   patPhoton_scR9.clear(); 
   patPhoton_scSigmaIEtaIEta.clear();
   patPhoton_scSigmaIEtaIPhi.clear();
   patPhoton_scSigmaIPhiIPhi.clear(); 
   patPhoton_full5x5_scEMax.clear();
   patPhoton_full5x5_scE2x2.clear(); 
   patPhoton_full5x5_scE3x3.clear();
   patPhoton_full5x5_scE5x5.clear(); 
   patPhoton_full5x5_scR9.clear(); 
   patPhoton_full5x5_scSigmaIEtaIEta.clear();
   patPhoton_full5x5_scSigmaIEtaIPhi.clear();
   patPhoton_full5x5_scSigmaIPhiIPhi.clear();
   patPhoton_HoE.clear();
   patPhoton_trkIso03.clear();
   patPhoton_ecalIso03.clear();
   patPhoton_hcalIso03.clear();
   patPhoton_trkIso04.clear();
   patPhoton_ecalIso04.clear();
   patPhoton_hcalIso04.clear();
   patPhoton_patParticleIso.clear();
   patPhoton_pfChargedHadronIso.clear();
   patPhoton_pfNeutralHadronIso.clear();
   patPhoton_pfPhotonIso.clear();
   patPhoton_pfPuChargedHadronIso.clear();
   patPhoton_nClusterOutsideMustache.clear(); 
   patPhoton_etOutsideMustache.clear(); 
   patPhoton_pfMVA.clear(); 
   patPhoton_pfDNN.clear(); 
   patPhoton_egmCutBasedPhotonIDloose.clear();
   patPhoton_egmCutBasedPhotonIDmedium.clear();
   patPhoton_egmCutBasedPhotonIDtight.clear();
   patPhoton_egmMVAPhotonIDmedium.clear();
   patPhoton_egmMVAPhotonIDtight.clear();

   patJet_index.clear();
   patJet_isCaloJet.clear();
   patJet_isJPTJet.clear();
   patJet_isPFJet.clear();
   patJet_isBasicJet.clear();
   patJet_charge.clear();
   patJet_energy.clear();  
   patJet_eta.clear();    
   patJet_phi.clear(); 
   patJet_area.clear(); 
   patJet_pt.clear();  
   patJet_uncorrectedEnergy.clear();  
   patJet_uncorrectedPt.clear();      
   patJet_energyFractionHadronic.clear();
   patJet_hadEnergyInHB.clear();
   patJet_hadEnergyInHO.clear();
   patJet_hadEnergyInHE.clear();
   patJet_hadEnergyInHF.clear();
   patJet_emEnergyInEB.clear();
   patJet_emEnergyInEE.clear();
   patJet_emEnergyInHF.clear();
   patJet_chargedHadronEnergyFraction.clear();
   patJet_neutralHadronEnergyFraction.clear();
   patJet_chargedEmEnergyFraction.clear();
   patJet_neutralEmEnergyFraction.clear();
   patJet_nOverlapMuons.clear();
   patJet_nOverlapTaus.clear(); 
   patJet_nOverlapElectrons.clear();
   patJet_overlapElectronIndices.clear();
   patJet_nOverlapPhotons.clear();
   patJet_overlapPhotonIndices.clear();
   patJet_photonEnergy.clear();
   patJet_photonEnergyFraction.clear();
   patJet_electronEnergy.clear();
   patJet_electronEnergyFraction.clear();
   patJet_muonEnergy.clear();
   patJet_muonEnergyFraction.clear();
   patJet_HFHadronEnergy.clear();       
   patJet_HFHadronEnergyFraction.clear();
   patJet_HFEMEnergy.clear();   
   patJet_HFEMEnergyFraction.clear();
   patJet_chargedMuEnergy.clear(); 
   patJet_chargedMuEnergyFraction.clear();
   patJet_hoEnergy.clear();  
   patJet_hoEnergyFraction.clear();
   patJet_nCandidates.clear();
   patJet_nCandInEcal.clear();
   patJet_nCandInEcalWithCharge.clear();
   patJet_candInEcal_charge.clear(); 
   patJet_candInEcal_ecalEnergy.clear();
   patJet_candInEcal_ecalEnergyFraction.clear();
   patJet_candInEcal_hcalEnergy.clear();
   patJet_candInEcal_hcalEnergyFraction.clear();
   patJet_candInEcal_eta.clear();
   patJet_candInEcal_phi.clear();  
   patJet_bTagScore_pfDeepCSV.clear();
   patJet_puIDScore.clear();
   patJet_puID.clear(); 
   patJet_qgL.clear(); 
   patJet_jmeCutBasedPFJetIDloose.clear(); 
   patJet_jmeCutBasedPFJetIDtight.clear(); 
   patJet_jmeCutBasedPFJetIDtightLepVeto.clear(); 
   
   superCluster_seedRawId.clear();
   superCluster_rawEnergy.clear(); 
   superCluster_rawESEnergy.clear(); 
   superCluster_energy.clear(); 
   superCluster_eta.clear(); 
   superCluster_phi.clear();  
   superCluster_etaWidth.clear();     
   superCluster_phiWidth.clear(); 
   superCluster_R.clear(); 
   superCluster_nPFClusters.clear(); 
   superCluster_ieta.clear(); 
   superCluster_iphi.clear();
   superCluster_iz.clear(); 
   superCluster_seedIndex.clear(); 
   superCluster_pfClustersIndex.clear(); 
   superCluster_e5x5.clear();
   superCluster_e2x2Ratio.clear();
   superCluster_e3x3Ratio.clear();
   superCluster_eMaxRatio.clear();
   superCluster_e2ndRatio.clear();
   superCluster_eTopRatio.clear();
   superCluster_eRightRatio.clear();
   superCluster_eBottomRatio.clear();
   superCluster_eLeftRatio.clear();
   superCluster_e2x5MaxRatio.clear();
   superCluster_e2x5TopRatio.clear();
   superCluster_e2x5RightRatio.clear();
   superCluster_e2x5BottomRatio.clear();
   superCluster_e2x5LeftRatio.clear();
   superCluster_swissCross.clear();
   superCluster_r9.clear();
   superCluster_sigmaIetaIeta.clear(); 
   superCluster_sigmaIetaIphi.clear(); 
   superCluster_sigmaIphiIphi.clear(); 
   superCluster_full5x5_e5x5.clear();
   superCluster_full5x5_e2x2Ratio.clear();
   superCluster_full5x5_e3x3Ratio.clear();
   superCluster_full5x5_eMaxRatio.clear();
   superCluster_full5x5_e2ndRatio.clear();
   superCluster_full5x5_eTopRatio.clear();
   superCluster_full5x5_eRightRatio.clear();
   superCluster_full5x5_eBottomRatio.clear();
   superCluster_full5x5_eLeftRatio.clear();
   superCluster_full5x5_e2x5MaxRatio.clear();
   superCluster_full5x5_e2x5TopRatio.clear();
   superCluster_full5x5_e2x5RightRatio.clear();
   superCluster_full5x5_e2x5BottomRatio.clear();
   superCluster_full5x5_e2x5LeftRatio.clear();
   superCluster_full5x5_swissCross.clear();
   superCluster_full5x5_r9.clear();
   superCluster_full5x5_sigmaIetaIeta.clear(); 
   superCluster_full5x5_sigmaIetaIphi.clear(); 
   superCluster_full5x5_sigmaIphiIphi.clear();   
   superCluster_psCluster_energy.clear();
   superCluster_psCluster_eta.clear();
   superCluster_psCluster_phi.clear();
   superCluster_nXtals.clear();
   superCluster_dR_genScore_MatchedIndex.clear();
   superCluster_dR_simScore_MatchedIndex.clear();
   superCluster_sim_nSharedXtals_MatchedIndex.clear();
   superCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   superCluster_sim_fraction_MatchedIndex.clear(); 
   superCluster_recoToSim_fraction_MatchedIndex.clear();
   superCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   superCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   superCluster_recoEnergy_sharedXtals_MatchedIndex.clear();    
   if(saveSuperCluster_){
      int nSuperClusters = nSuperClustersEB + nSuperClustersEE;
      superCluster_seedIndex.resize(nSuperClusters);     
      superCluster_pfClustersIndex.resize(nSuperClusters);
      superCluster_psCluster_energy.resize(nSuperClustersEE);
      superCluster_psCluster_eta.resize(nSuperClustersEE);
      superCluster_psCluster_phi.resize(nSuperClustersEE);
      superCluster_dR_genScore.resize(nSuperClusters);
      superCluster_dR_simScore.resize(nSuperClusters);
      superCluster_sim_nSharedXtals.resize(nSuperClusters);
      superCluster_sim_fraction_noHitsFraction.resize(nSuperClusters);
      superCluster_sim_fraction.resize(nSuperClusters);
      superCluster_recoToSim_fraction.resize(nSuperClusters);
      superCluster_recoToSim_fraction_sharedXtals.resize(nSuperClusters);
      superCluster_simEnergy_sharedXtals.resize(nSuperClusters);
      superCluster_recoEnergy_sharedXtals.resize(nSuperClusters);  
      superCluster_simPU_nSharedXtals.resize(nSuperClusters);
      superCluster_simEnergy_sharedXtalsPU.resize(nSuperClusters);
      superCluster_recoEnergy_sharedXtalsPU.resize(nSuperClusters);   
      superCluster_simEnergy_noHitsFraction_sharedXtalsPU.resize(nSuperClusters);
      superCluster_recoEnergy_noHitsFraction_sharedXtalsPU.resize(nSuperClusters); 
      superCluster_simOOTPU_nSharedXtals.resize(nSuperClusters);
      superCluster_simEnergy_sharedXtalsOOTPU.resize(nSuperClusters);
      superCluster_recoEnergy_sharedXtalsOOTPU.resize(nSuperClusters);   
      superCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU.resize(nSuperClusters);
      superCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU.resize(nSuperClusters); 
   }
   
   retunedSuperCluster_seedRawId.clear();
   retunedSuperCluster_rawEnergy.clear(); 
   retunedSuperCluster_rawESEnergy.clear(); 
   retunedSuperCluster_energy.clear(); 
   retunedSuperCluster_eta.clear(); 
   retunedSuperCluster_phi.clear();  
   retunedSuperCluster_etaWidth.clear();     
   retunedSuperCluster_phiWidth.clear(); 
   retunedSuperCluster_R.clear(); 
   retunedSuperCluster_nPFClusters.clear(); 
   retunedSuperCluster_ieta.clear(); 
   retunedSuperCluster_iphi.clear();
   retunedSuperCluster_iz.clear(); 
   retunedSuperCluster_seedIndex.clear(); 
   retunedSuperCluster_pfClustersIndex.clear(); 
   retunedSuperCluster_e5x5.clear();
   retunedSuperCluster_e2x2Ratio.clear();
   retunedSuperCluster_e3x3Ratio.clear();
   retunedSuperCluster_eMaxRatio.clear();
   retunedSuperCluster_e2ndRatio.clear();
   retunedSuperCluster_eTopRatio.clear();
   retunedSuperCluster_eRightRatio.clear();
   retunedSuperCluster_eBottomRatio.clear();
   retunedSuperCluster_eLeftRatio.clear();
   retunedSuperCluster_e2x5MaxRatio.clear();
   retunedSuperCluster_e2x5TopRatio.clear();
   retunedSuperCluster_e2x5RightRatio.clear();
   retunedSuperCluster_e2x5BottomRatio.clear();
   retunedSuperCluster_e2x5LeftRatio.clear();
   retunedSuperCluster_swissCross.clear();
   retunedSuperCluster_r9.clear();
   retunedSuperCluster_sigmaIetaIeta.clear(); 
   retunedSuperCluster_sigmaIetaIphi.clear(); 
   retunedSuperCluster_sigmaIphiIphi.clear(); 
   retunedSuperCluster_full5x5_e5x5.clear();
   retunedSuperCluster_full5x5_e2x2Ratio.clear();
   retunedSuperCluster_full5x5_e3x3Ratio.clear();
   retunedSuperCluster_full5x5_eMaxRatio.clear();
   retunedSuperCluster_full5x5_e2ndRatio.clear();
   retunedSuperCluster_full5x5_eTopRatio.clear();
   retunedSuperCluster_full5x5_eRightRatio.clear();
   retunedSuperCluster_full5x5_eBottomRatio.clear();
   retunedSuperCluster_full5x5_eLeftRatio.clear();
   retunedSuperCluster_full5x5_e2x5MaxRatio.clear();
   retunedSuperCluster_full5x5_e2x5TopRatio.clear();
   retunedSuperCluster_full5x5_e2x5RightRatio.clear();
   retunedSuperCluster_full5x5_e2x5BottomRatio.clear();
   retunedSuperCluster_full5x5_e2x5LeftRatio.clear();
   retunedSuperCluster_full5x5_swissCross.clear();
   retunedSuperCluster_full5x5_r9.clear();
   retunedSuperCluster_full5x5_sigmaIetaIeta.clear(); 
   retunedSuperCluster_full5x5_sigmaIetaIphi.clear(); 
   retunedSuperCluster_full5x5_sigmaIphiIphi.clear();    
   retunedSuperCluster_psCluster_energy.clear();
   retunedSuperCluster_psCluster_eta.clear();
   retunedSuperCluster_psCluster_phi.clear();
   retunedSuperCluster_nXtals.clear(); 
   retunedSuperCluster_dR_genScore_MatchedIndex.clear();
   retunedSuperCluster_dR_simScore_MatchedIndex.clear();
   retunedSuperCluster_sim_nSharedXtals_MatchedIndex.clear();
   retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   retunedSuperCluster_sim_fraction_MatchedIndex.clear(); 
   retunedSuperCluster_recoToSim_fraction_MatchedIndex.clear();
   retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex.clear();  
   if(saveRetunedSC_){
      int nRetunedSuperClusters = nRetunedSuperClustersEB + nRetunedSuperClustersEE;
      retunedSuperCluster_seedIndex.resize(nRetunedSuperClusters);     
      retunedSuperCluster_pfClustersIndex.resize(nRetunedSuperClusters);
      retunedSuperCluster_psCluster_energy.resize(nRetunedSuperClustersEE);
      retunedSuperCluster_psCluster_eta.resize(nRetunedSuperClustersEE);
      retunedSuperCluster_psCluster_phi.resize(nRetunedSuperClustersEE);
      retunedSuperCluster_dR_genScore.resize(nRetunedSuperClusters);
      retunedSuperCluster_dR_simScore.resize(nRetunedSuperClusters);
      retunedSuperCluster_sim_nSharedXtals.resize(nRetunedSuperClusters);
      retunedSuperCluster_sim_fraction_noHitsFraction.resize(nRetunedSuperClusters);
      retunedSuperCluster_sim_fraction.resize(nRetunedSuperClusters);
      retunedSuperCluster_recoToSim_fraction.resize(nRetunedSuperClusters);
      retunedSuperCluster_recoToSim_fraction_sharedXtals.resize(nRetunedSuperClusters);
      retunedSuperCluster_simEnergy_sharedXtals.resize(nRetunedSuperClusters);
      retunedSuperCluster_recoEnergy_sharedXtals.resize(nRetunedSuperClusters);  
      retunedSuperCluster_simPU_nSharedXtals.resize(nRetunedSuperClusters);
      retunedSuperCluster_simEnergy_sharedXtalsPU.resize(nRetunedSuperClusters);
      retunedSuperCluster_recoEnergy_sharedXtalsPU.resize(nRetunedSuperClusters);   
      retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU.resize(nRetunedSuperClusters);
      retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU.resize(nRetunedSuperClusters); 
      retunedSuperCluster_simOOTPU_nSharedXtals.resize(nRetunedSuperClusters);
      retunedSuperCluster_simEnergy_sharedXtalsOOTPU.resize(nRetunedSuperClusters);
      retunedSuperCluster_recoEnergy_sharedXtalsOOTPU.resize(nRetunedSuperClusters);   
      retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU.resize(nRetunedSuperClusters);
      retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU.resize(nRetunedSuperClusters);  
   }

   deepSuperCluster_seedRawId.clear();
   deepSuperCluster_rawEnergy.clear(); 
   deepSuperCluster_rawESEnergy.clear(); 
   deepSuperCluster_energy.clear(); 
   deepSuperCluster_eta.clear(); 
   deepSuperCluster_phi.clear();  
   deepSuperCluster_etaWidth.clear();     
   deepSuperCluster_phiWidth.clear(); 
   deepSuperCluster_R.clear(); 
   deepSuperCluster_nPFClusters.clear(); 
   deepSuperCluster_ieta.clear(); 
   deepSuperCluster_iphi.clear();
   deepSuperCluster_iz.clear(); 
   deepSuperCluster_seedIndex.clear(); 
   deepSuperCluster_pfClustersIndex.clear(); 
   deepSuperCluster_e5x5.clear();
   deepSuperCluster_e2x2Ratio.clear();
   deepSuperCluster_e3x3Ratio.clear();
   deepSuperCluster_eMaxRatio.clear();
   deepSuperCluster_e2ndRatio.clear();
   deepSuperCluster_eTopRatio.clear();
   deepSuperCluster_eRightRatio.clear();
   deepSuperCluster_eBottomRatio.clear();
   deepSuperCluster_eLeftRatio.clear();
   deepSuperCluster_e2x5MaxRatio.clear();
   deepSuperCluster_e2x5TopRatio.clear();
   deepSuperCluster_e2x5RightRatio.clear();
   deepSuperCluster_e2x5BottomRatio.clear();
   deepSuperCluster_e2x5LeftRatio.clear();
   deepSuperCluster_swissCross.clear();
   deepSuperCluster_r9.clear();
   deepSuperCluster_sigmaIetaIeta.clear(); 
   deepSuperCluster_sigmaIetaIphi.clear(); 
   deepSuperCluster_sigmaIphiIphi.clear(); 
   deepSuperCluster_full5x5_e5x5.clear();
   deepSuperCluster_full5x5_e2x2Ratio.clear();
   deepSuperCluster_full5x5_e3x3Ratio.clear();
   deepSuperCluster_full5x5_eMaxRatio.clear();
   deepSuperCluster_full5x5_e2ndRatio.clear();
   deepSuperCluster_full5x5_eTopRatio.clear();
   deepSuperCluster_full5x5_eRightRatio.clear();
   deepSuperCluster_full5x5_eBottomRatio.clear();
   deepSuperCluster_full5x5_eLeftRatio.clear();
   deepSuperCluster_full5x5_e2x5MaxRatio.clear();
   deepSuperCluster_full5x5_e2x5TopRatio.clear();
   deepSuperCluster_full5x5_e2x5RightRatio.clear();
   deepSuperCluster_full5x5_e2x5BottomRatio.clear();
   deepSuperCluster_full5x5_e2x5LeftRatio.clear();
   deepSuperCluster_full5x5_swissCross.clear();
   deepSuperCluster_full5x5_r9.clear();
   deepSuperCluster_full5x5_sigmaIetaIeta.clear(); 
   deepSuperCluster_full5x5_sigmaIetaIphi.clear(); 
   deepSuperCluster_full5x5_sigmaIphiIphi.clear();      
   deepSuperCluster_psCluster_energy.clear();
   deepSuperCluster_psCluster_eta.clear();
   deepSuperCluster_psCluster_phi.clear();
   deepSuperCluster_nXtals.clear();
   deepSuperCluster_dR_genScore_MatchedIndex.clear();
   deepSuperCluster_dR_simScore_MatchedIndex.clear();
   deepSuperCluster_sim_nSharedXtals_MatchedIndex.clear();
   deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   deepSuperCluster_sim_fraction_MatchedIndex.clear(); 
   deepSuperCluster_recoToSim_fraction_MatchedIndex.clear();
   deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   deepSuperCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex.clear();  
   if(saveDeepSC_){ 
      int nDeepSuperClusters = nDeepSuperClustersEB + nDeepSuperClustersEE; 
      deepSuperCluster_seedIndex.resize(nDeepSuperClusters); 
      deepSuperCluster_pfClustersIndex.resize(nDeepSuperClusters);
      deepSuperCluster_psCluster_energy.resize(nDeepSuperClustersEE);
      deepSuperCluster_psCluster_eta.resize(nDeepSuperClustersEE);
      deepSuperCluster_psCluster_phi.resize(nDeepSuperClustersEE);
      deepSuperCluster_dR_genScore.resize(nDeepSuperClusters);
      deepSuperCluster_dR_simScore.resize(nDeepSuperClusters);
      deepSuperCluster_sim_nSharedXtals.resize(nDeepSuperClusters);
      deepSuperCluster_sim_fraction_noHitsFraction.resize(nDeepSuperClusters);
      deepSuperCluster_sim_fraction.resize(nDeepSuperClusters);
      deepSuperCluster_recoToSim_fraction.resize(nDeepSuperClusters);
      deepSuperCluster_recoToSim_fraction_sharedXtals.resize(nDeepSuperClusters);
      deepSuperCluster_simEnergy_sharedXtals.resize(nDeepSuperClusters);
      deepSuperCluster_recoEnergy_sharedXtals.resize(nDeepSuperClusters);   
      deepSuperCluster_simPU_nSharedXtals.resize(nDeepSuperClusters);
      deepSuperCluster_simEnergy_sharedXtalsPU.resize(nDeepSuperClusters);
      deepSuperCluster_recoEnergy_sharedXtalsPU.resize(nDeepSuperClusters);   
      deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU.resize(nDeepSuperClusters);
      deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU.resize(nDeepSuperClusters); 
      deepSuperCluster_simOOTPU_nSharedXtals.resize(nDeepSuperClusters);
      deepSuperCluster_simEnergy_sharedXtalsOOTPU.resize(nDeepSuperClusters);
      deepSuperCluster_recoEnergy_sharedXtalsOOTPU.resize(nDeepSuperClusters);   
      deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU.resize(nDeepSuperClusters);
      deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU.resize(nDeepSuperClusters); 
   }   
}

double RecoSimDumper::ptFast(const double energy, const math::XYZPoint& position, const math::XYZPoint& origin)
{
   const auto v = position - origin;
   return energy * std::sqrt(v.perp2() / v.mag2()); 
}

int RecoSimDumper::getGenMother(const reco::GenParticle* genParticle)
{
   int genMotherIndex=-1;
   if(genParticle->numberOfMothers()>0) genMotherIndex = genParticle->motherRef(0).key();
   return genMotherIndex; 
}

int RecoSimDumper::getGenStatusFlag(const reco::GenParticle* genParticle)
{
   int statusFlag = 
   genParticle->statusFlags().isLastCopyBeforeFSR()                  * 16384 +
   genParticle->statusFlags().isLastCopy()                           * 8192  +
   genParticle->statusFlags().isFirstCopy()                          * 4096  +
   genParticle->statusFlags().fromHardProcessBeforeFSR()             * 2048  +
   genParticle->statusFlags().isDirectHardProcessTauDecayProduct()   * 1024  +
   genParticle->statusFlags().isHardProcessTauDecayProduct()         * 512   +
   genParticle->statusFlags().fromHardProcess()                      * 256   +  
   genParticle->statusFlags().isHardProcess()                        * 128   +
   genParticle->statusFlags().isDirectHadronDecayProduct()           * 64    +
   genParticle->statusFlags().isDirectPromptTauDecayProduct()        * 32    +
   genParticle->statusFlags().isDirectTauDecayProduct()              * 16    +
   genParticle->statusFlags().isPromptTauDecayProduct()              * 8     +
   genParticle->statusFlags().isTauDecayProduct()                    * 4     +
   genParticle->statusFlags().isDecayedLeptonHadron()                * 2     +
   genParticle->statusFlags().isPrompt()                             * 1      ;
   
   return statusFlag; 
}

void RecoSimDumper::printGenStatusFlag(const reco::GenParticle* genParticle)
{
   
   std::cout << " -- Flags: "; 
   if(genParticle->statusFlags().isLastCopyBeforeFSR()) std::cout << "isLastCopyBeforeFSR() && ";
   if(genParticle->statusFlags().isLastCopy()) std::cout << "isLastCopy() && ";     
   if(genParticle->statusFlags().isFirstCopy()) std::cout << "isFirstCopy() && ";
   if(genParticle->statusFlags().fromHardProcessBeforeFSR()) std::cout << "fromHardProcessBeforeFSR() && ";
   if(genParticle->statusFlags().isDirectHardProcessTauDecayProduct()) std::cout << "isDirectHardProcessTauDecayProduct() && ";
   if(genParticle->statusFlags().isHardProcessTauDecayProduct()) std::cout << "isHardProcessTauDecayProduct() && ";
   if(genParticle->statusFlags().fromHardProcess()) std::cout << "fromHardProcess() && "; 
   if(genParticle->statusFlags().isHardProcess()) std::cout << "isHardProcess() && ";
   if(genParticle->statusFlags().isDirectHadronDecayProduct()) std::cout << "isDirectHadronDecayProduct() && ";
   if(genParticle->statusFlags().isDirectPromptTauDecayProduct()) std::cout << "isDirectPromptTauDecayProduct() && ";
   if(genParticle->statusFlags().isDirectTauDecayProduct()) std::cout << "isDirectTauDecayProduct() && ";
   if(genParticle->statusFlags().isPromptTauDecayProduct()) std::cout << "isPromptTauDecayProduct() && ";
   if(genParticle->statusFlags().isTauDecayProduct()) std::cout << "isTauDecayProduct() && ";
   if(genParticle->statusFlags().isDecayedLeptonHadron()) std::cout << "isDecayedLeptonHadron() && ";
   if(genParticle->statusFlags().isPrompt()) std::cout << "isPrompt() ";
   std::cout << "-- " << getGenStatusFlag(genParticle) << " -- ";
}

int RecoSimDumper::getGenParton(const std::vector<reco::GenParticle>* genParticles, const int genIndex)
{
   int genMotherIndex = -1;
   genMotherIndices.clear();
   if(genIndex>=0) genMotherIndex = genIndex;
   while(genMotherIndex>=0)
   {  
      genMotherIndex = getGenMother(&genParticles->at(genMotherIndex));
      if(genMotherIndex<0) break; 
      //printGenStatusFlag(&genParticles->at(genMotherIndex));
      genMotherIndices.push_back(genMotherIndex);
   }
    
   int genPartonIndex = -1;
   if(genMotherIndices.size()>=1) genPartonIndex = genMotherIndices.at(genMotherIndices.size()-1);
   return genPartonIndex; 
}

std::vector<std::pair<DetId,std::pair<float,float>>>* RecoSimDumper::getSharedHitsAndEnergies(const std::vector<std::pair<DetId, float> >* hitsAndEnergies1, const std::vector<std::pair<DetId, float> >* hitsAndEnergies2)
{
    
    std::vector<std::pair<DetId,pair<float,float>>>* sharedHitsAndEnergies = new std::vector<std::pair<DetId,pair<float,float>>>;
    for(unsigned int i = 0; i < hitsAndEnergies1->size(); i++)
    {  
        for(unsigned int j = 0; j < hitsAndEnergies2->size(); j++)
        {
            if(hitsAndEnergies1->at(i).first.rawId()==hitsAndEnergies2->at(j).first.rawId()) sharedHitsAndEnergies->push_back(make_pair(hitsAndEnergies1->at(i).first, make_pair(hitsAndEnergies1->at(i).second,hitsAndEnergies2->at(j).second)));
        }  
    }
    return sharedHitsAndEnergies;
}

std::vector<std::pair<DetId, float> >* RecoSimDumper::getHitsAndEnergiesCaloPart(const CaloParticle* iCaloParticle, float simHitEnergy_cut)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_CaloPart_tmp = new std::vector<std::pair<DetId, float> >;
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_tmp = new std::vector<std::pair<DetId, float> >;
    std::map<DetId, float> HitsAndEnergies_map;
    
    const auto& simClusters = iCaloParticle->simClusters();
    for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
        auto simCluster = simClusters[iSC];  
        auto hits_and_energies = simCluster->hits_and_energies();
        for(unsigned int i = 0; i < hits_and_energies.size(); i++){ 
            if(hits_and_energies[i].second < simHitEnergy_cut) continue; 
            HitsAndEnergies_tmp->push_back(make_pair(DetId(hits_and_energies[i].first),hits_and_energies[i].second));  
        }  
    }

    for(unsigned int i = 0; i < HitsAndEnergies_tmp->size(); i++){  
        if (HitsAndEnergies_map.find(HitsAndEnergies_tmp->at(i).first) == HitsAndEnergies_map.end()) {
            HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]=HitsAndEnergies_tmp->at(i).second;      
        }else{
            HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]=HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]+HitsAndEnergies_tmp->at(i).second; 
        }
    }

    for(auto const& hit : HitsAndEnergies_map) 
         HitsAndEnergies_CaloPart_tmp->push_back(make_pair(hit.first,hit.second));

    return HitsAndEnergies_CaloPart_tmp;
}

std::vector<std::pair<DetId, float> >* RecoSimDumper::getHitsAndEnergiesBC(reco::CaloCluster* iPFCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_tmp = new std::vector<std::pair<DetId, float> >;
    const std::vector<std::pair<DetId,float> > &hitsAndFractions = iPFCluster->hitsAndFractions();
    for(unsigned int i = 0; i < hitsAndFractions.size(); i++){
        if(hitsAndFractions.at(i).first.subdetId()==EcalBarrel){
           HitsAndEnergies_tmp->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*recHitsEB->find(hitsAndFractions.at(i).first)).energy()));
        }else if(hitsAndFractions.at(i).first.subdetId()==EcalEndcap){
           HitsAndEnergies_tmp->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*recHitsEE->find(hitsAndFractions.at(i).first)).energy()));
        }
    }

    return HitsAndEnergies_tmp;
}


std::vector<std::pair<DetId, float> >* RecoSimDumper::getHitsAndEnergiesSC(const reco::SuperCluster* iSuperCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_SuperCluster_tmp = new std::vector<std::pair<DetId, float> >;
    std::map<DetId, float> HitsAndEnergies_map;

    for(reco::CaloCluster_iterator iBC = iSuperCluster->clustersBegin(); iBC != iSuperCluster->clustersEnd(); ++iBC){
        const std::vector<std::pair<DetId,float> > &seedrechits = ( *iBC )->hitsAndFractions();
        for(unsigned int i = 0; i < seedrechits.size(); i++){  
            if(seedrechits.at(i).first.subdetId()==EcalBarrel){   
                if (HitsAndEnergies_map.find(seedrechits.at(i).first) == HitsAndEnergies_map.end()) {
                    HitsAndEnergies_map[seedrechits.at(i).first]=seedrechits.at(i).second * (*recHitsEB->find(seedrechits.at(i).first)).energy();    
                }else{
                    HitsAndEnergies_map[seedrechits.at(i).first]=HitsAndEnergies_map[seedrechits.at(i).first]+seedrechits.at(i).second * (*recHitsEB->find(seedrechits.at(i).first)).energy();
                } 
            }else if(seedrechits.at(i).first.subdetId()==EcalEndcap){   
                if (HitsAndEnergies_map.find(seedrechits.at(i).first) == HitsAndEnergies_map.end()) {
                    HitsAndEnergies_map[seedrechits.at(i).first]=seedrechits.at(i).second * (*recHitsEE->find(seedrechits.at(i).first)).energy();   
                }else{
                    HitsAndEnergies_map[seedrechits.at(i).first]=HitsAndEnergies_map[seedrechits.at(i).first]+seedrechits.at(i).second * (*recHitsEE->find(seedrechits.at(i).first)).energy();
                } 
            }
        }                      
    } 

    for(auto const& hit : HitsAndEnergies_map) 
        HitsAndEnergies_SuperCluster_tmp->push_back(make_pair(hit.first,hit.second));

    return HitsAndEnergies_SuperCluster_tmp;
}

std::pair<double,double> RecoSimDumper::calculateCovariances(const reco::PFCluster* pfCluster, const EcalRecHitCollection* recHits, const CaloSubdetectorGeometry* geometry) {

  double etaWidth = 0.;
  double phiWidth = 0.;
  double numeratorEtaWidth = 0;
  double numeratorPhiWidth = 0;

  double clEnergy = pfCluster->energy();
  double denominator = clEnergy;

  double clEta = pfCluster->position().eta();
  double clPhi = pfCluster->position().phi();

  const std::vector<std::pair<DetId, float> >& detId = pfCluster->hitsAndFractions();
  // Loop over recHits associated with the given SuperCluster
  for (std::vector<std::pair<DetId, float> >::const_iterator hit = detId.begin(); hit != detId.end(); ++hit) {
    EcalRecHitCollection::const_iterator rHit = recHits->find((*hit).first);
    //FIXME: THIS IS JUST A WORKAROUND A FIX SHOULD BE APPLIED
    if (rHit == recHits->end()) {
      continue;
    }
    auto this_cell = geometry->getGeometry(rHit->id());
    if (this_cell == nullptr) {
      //edm::LogInfo("SuperClusterShapeAlgo") << "pointer to the cell in Calculate_Covariances is NULL!";
      continue;
    }
    GlobalPoint position = this_cell->getPosition();
    //take into account energy fractions
    double energyHit = rHit->energy() * hit->second; 

    //form differences
    double dPhi = position.phi() - clPhi;
    if (dPhi > +Geom::pi()) {
      dPhi = Geom::twoPi() - dPhi;
    }
    if (dPhi < -Geom::pi()) {
      dPhi = Geom::twoPi() + dPhi;
    }

    double dEta = position.eta() - clEta;

    if (energyHit > 0) {
      numeratorEtaWidth += energyHit * dEta * dEta;
      numeratorPhiWidth += energyHit * dPhi * dPhi;
    }

    etaWidth = sqrt(numeratorEtaWidth / denominator);
    phiWidth = sqrt(numeratorPhiWidth / denominator);
  }

  return std::make_pair(etaWidth,phiWidth);
}

std::vector<float> RecoSimDumper::getShowerShapes(reco::CaloCluster* caloBC, const EcalRecHitCollection* recHits, const CaloTopology *topology)
{
    std::vector<float> shapes;
    shapes.resize(38); 
    locCov_ = EcalClusterTools::localCovariances(*caloBC, recHits, topology);
    full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(*caloBC, recHits, topology);
    
    float e5x5 = EcalClusterTools::e5x5(*caloBC, recHits, topology); // e5x5
    float e3x3 = EcalClusterTools::e3x3(*caloBC, recHits, topology); // e3x3
    float eMax = EcalClusterTools::eMax(*caloBC, recHits); // eMax
    float eTop = EcalClusterTools::eTop(*caloBC, recHits, topology); // eTop 
    float eRight = EcalClusterTools::eRight(*caloBC, recHits, topology); // eRight
    float eBottom = EcalClusterTools::eBottom(*caloBC, recHits, topology); // eBottom
    float eLeft = EcalClusterTools::eLeft(*caloBC, recHits, topology); // eLeft
    float e4 = eTop + eRight + eBottom + eLeft;

    shapes[0] = e5x5;
    shapes[1] = EcalClusterTools::e2x2(*caloBC, recHits, topology)/e5x5; // e2x2/e5x5
    shapes[2] = EcalClusterTools::e3x3(*caloBC, recHits, topology)/e5x5; // e3x3/e5x5
    shapes[3] = EcalClusterTools::eMax(*caloBC,  recHits)/e5x5; // eMax/e5x5
    shapes[4] = EcalClusterTools::e2nd(*caloBC, recHits)/e5x5; // e2nd/e5x5
    shapes[5] = EcalClusterTools::eTop(*caloBC, recHits, topology)/e5x5; // eTop/e5x5 
    shapes[6] = EcalClusterTools::eRight(*caloBC, recHits, topology)/e5x5; // eRight/e5x5
    shapes[7] = EcalClusterTools::eBottom(*caloBC, recHits, topology)/e5x5; // eBottom/e5x5
    shapes[8] = EcalClusterTools::eLeft(*caloBC, recHits, topology)/e5x5; // eLeft/e5x5
    shapes[9] = EcalClusterTools::e2x5Max(*caloBC, recHits, topology)/e5x5; // e2x5Max/e5x5
    shapes[10] = EcalClusterTools::e2x5Top(*caloBC, recHits, topology)/e5x5; // e2x5Top/e5x5  
    shapes[11] = EcalClusterTools::e2x5Right(*caloBC, recHits, topology)/e5x5; // e2x5Bottom/e5x5  
    shapes[12] = EcalClusterTools::e2x5Bottom(*caloBC, recHits, topology)/e5x5; // e2x5Left/e5x5  
    shapes[13] = EcalClusterTools::e2x5Left(*caloBC, recHits, topology)/e5x5; // e2x5Right/e5x5   
    shapes[14] = 1.-e4/eMax; // swissCross 
    shapes[15] = e3x3/caloBC->energy(); // r9     
    shapes[16] = sqrt(locCov_[0]); // sigmaIetaIeta 
    shapes[17] = locCov_[1]; // sigmaIetaIphi
    shapes[18] = !edm::isFinite(locCov_[2]) ? 0. : sqrt(locCov_[2]); // sigmaIphiIphi 

    //if(std::isnan(shapes.at(17))) std::cout << shapes[0] << " - " << shapes[1] << " - " << shapes[2] << " - " << shapes[3] << " - " << shapes[4] << " - " << shapes[5]  << " - " << shapes[6] << " - " << shapes[7] << " - " << shapes[8] << " - " << shapes[9] << " - " << shapes[10] << " - " << shapes[11] << " - " << shapes[12] << " - " << shapes[13] << " - " << shapes[14] << " - " << shapes[15]  << " - " << shapes[16] << " - " << shapes[17] << " - " << shapes[18] << std::endl;

    // full_5x5 variables
    float full5x5_e5x5 = noZS::EcalClusterTools::e5x5(*caloBC, recHits, topology); // e5x5
    float full5x5_e3x3 = noZS::EcalClusterTools::e3x3(*caloBC, recHits, topology); // e3x3
    float full5x5_eMax = noZS::EcalClusterTools::eMax(*caloBC, recHits); // eMax
    float full5x5_eTop = noZS::EcalClusterTools::eTop(*caloBC, recHits, topology); // eTop 
    float full5x5_eRight = noZS::EcalClusterTools::eRight(*caloBC, recHits, topology); // eRight
    float full5x5_eBottom = noZS::EcalClusterTools::eBottom(*caloBC, recHits, topology); // eBottom
    float full5x5_eLeft = noZS::EcalClusterTools::eLeft(*caloBC, recHits, topology); // eLeft
    float full5x5_e4 = full5x5_eTop + full5x5_eRight + full5x5_eBottom + full5x5_eLeft;

    shapes[19] = full5x5_e5x5;
    shapes[20] = noZS::EcalClusterTools::e2x2(*caloBC, recHits, topology)/full5x5_e5x5; // e2x2/e5x5
    shapes[21] = noZS::EcalClusterTools::e3x3(*caloBC, recHits, topology)/full5x5_e5x5; // e3x3/e5x5
    shapes[22] = noZS::EcalClusterTools::eMax(*caloBC, recHits)/full5x5_e5x5; // eMax/e5x5
    shapes[23] = noZS::EcalClusterTools::e2nd(*caloBC, recHits)/full5x5_e5x5; // e2nd/e5x5
    shapes[24] = noZS::EcalClusterTools::eTop(*caloBC, recHits, topology)/full5x5_e5x5; // eTop/e5x5 
    shapes[25] = noZS::EcalClusterTools::eRight(*caloBC, recHits, topology)/full5x5_e5x5; // eRight/e5x5
    shapes[26] = noZS::EcalClusterTools::eBottom(*caloBC, recHits, topology)/full5x5_e5x5; // eBottom/e5x5
    shapes[27] = noZS::EcalClusterTools::eLeft(*caloBC, recHits, topology)/full5x5_e5x5; // eLeft/e5x5
    shapes[28] = noZS::EcalClusterTools::e2x5Max(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Max/e5x5
    shapes[29] = noZS::EcalClusterTools::e2x5Top(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Top/e5x5  
    shapes[30] = noZS::EcalClusterTools::e2x5Right(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Bottom/e5x5  
    shapes[31] = noZS::EcalClusterTools::e2x5Bottom(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Left/e5x5  
    shapes[32] = noZS::EcalClusterTools::e2x5Left(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Right/e5x5   
    shapes[33] = 1.-full5x5_e4/full5x5_eMax; // swissCross 
    shapes[34] = full5x5_e3x3/caloBC->energy(); // r9
    shapes[35] = sqrt(full5x5_locCov_[0]); // sigmaIetaIeta        
    shapes[36] = full5x5_locCov_[1]; // sigmaIetaIphi          
    shapes[37] = !edm::isFinite(full5x5_locCov_[2]) ? 0. : sqrt(full5x5_locCov_[2]); // sigmaIphiIphi

    for(unsigned iVar=0; iVar<shapes.size(); iVar++)
        if(std::isnan(shapes.at(iVar))) std::cout << "showerShape = " << iVar << " ---> NAN " << std::endl;  

    return shapes; 
}

std::vector<double> RecoSimDumper::getScores(const reco::PFCluster* pfCluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE)
{
    std::vector<double> scores;
    scores.resize(9);

    double nSharedXtals=0;
    double simFraction_noHitsFraction=0.;
    double simFraction=0.; 
    double recoToSim=0.; 
    double recoToSim_shared=0.; 

    double simEnergy=0.;
    double simEnergy_shared=0.;
    double simEnergy_shared_noHitsFraction=0.;
    double recoEnergy=pfCluster->energy();
    double recoEnergy_shared=0.;
    double recoEnergy_shared_noHitsFraction=0.;
   
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart)
        simEnergy+=hit_CaloPart.second;
   
    const std::vector<std::pair<DetId,float> >* hitsAndFractions = &pfCluster->hitsAndFractions();
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){
        for(const std::pair<DetId, float>& hit_Cluster : *hitsAndFractions){     
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
               nSharedXtals+=1.;
               simEnergy_shared+=hit_CaloPart.second*hit_Cluster.second;
               simEnergy_shared_noHitsFraction+=hit_CaloPart.second;
               if(hit_Cluster.first.subdetId()==EcalBarrel){ 
                  recoEnergy_shared+=hit_Cluster.second*(*recHitsEB->find(hit_Cluster.first)).energy();
                  recoEnergy_shared_noHitsFraction+=(*recHitsEB->find(hit_Cluster.first)).energy();
               }
               if(hit_Cluster.first.subdetId()==EcalEndcap){ 
                  recoEnergy_shared+=hit_Cluster.second*(*recHitsEE->find(hit_Cluster.first)).energy(); 
                  recoEnergy_shared_noHitsFraction+=(*recHitsEE->find(hit_Cluster.first)).energy(); 
               } 
            } 
        }
    }

    if(nSharedXtals<0.) nSharedXtals = -1.; 

    if(simEnergy>0.) simFraction_noHitsFraction = simEnergy_shared_noHitsFraction/simEnergy;
    else simFraction_noHitsFraction = -1.; 

    if(simEnergy>0.) simFraction = simEnergy_shared/simEnergy;
    else simFraction = -1.;
 
    if(simEnergy_shared>0.) recoToSim = recoEnergy/simEnergy_shared;
    else recoToSim = -1.; 

    if(simEnergy_shared>0.) recoToSim_shared = recoEnergy_shared/simEnergy_shared;
    else recoToSim_shared = -1.;  

    if(simEnergy_shared<0) simEnergy_shared = -1.;
    if(recoEnergy_shared<0) recoEnergy_shared = -1.;
    if(simEnergy_shared_noHitsFraction<0) simEnergy_shared_noHitsFraction = -1.;
    if(recoEnergy_shared_noHitsFraction<0) recoEnergy_shared_noHitsFraction = -1.;
    
    scores[0] = (double)nSharedXtals;
    scores[1] = simFraction_noHitsFraction;
    scores[2] = simFraction;
    scores[3] = recoToSim;
    scores[4] = recoToSim_shared;
    scores[5] = simEnergy_shared;
    scores[6] = recoEnergy_shared;
    scores[7] = simEnergy_shared_noHitsFraction;
    scores[8] = recoEnergy_shared_noHitsFraction;

    for(unsigned iVar=0; iVar<scores.size(); iVar++)
        if(std::isnan(scores.at(iVar))) std::cout << "score = " << iVar << " ---> NAN " << std::endl; 

    return scores;
}

std::vector<double> RecoSimDumper::getScores(const reco::SuperCluster* superCluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE)
{
    std::vector<double> scores;
    scores.resize(9);

    double nSharedXtals=0;
    double simFraction_noHitsFraction=0.;
    double simFraction=0.; 
    double recoToSim=0.; 
    double recoToSim_shared=0.; 

    double simEnergy=0.;
    double simEnergy_shared=0.;
    double simEnergy_shared_noHitsFraction=0.;
    double recoEnergy=superCluster->energy();
    double recoEnergy_shared=0.;
    double recoEnergy_shared_noHitsFraction=0.; 
   
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart)
        simEnergy+=hit_CaloPart.second;
   
    std::vector<DetId> superCluster_IDs;
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){      
        for(reco::CaloCluster_iterator iBC = superCluster->clustersBegin(); iBC != superCluster->clustersEnd(); ++iBC){
            const std::vector<std::pair<DetId,float> >* hitsAndFractions = &( *iBC )->hitsAndFractions();
            for(const std::pair<DetId, float>& hit_Cluster : *hitsAndFractions){      
                if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
                   simEnergy_shared+=hit_CaloPart.second*hit_Cluster.second;
                   simEnergy_shared_noHitsFraction+=hit_CaloPart.second;
                   if(hit_Cluster.first.subdetId()==EcalBarrel){ 
                      recoEnergy_shared+=hit_Cluster.second*(*recHitsEB->find(hit_Cluster.first)).energy();
                      recoEnergy_shared_noHitsFraction+=(*recHitsEB->find(hit_Cluster.first)).energy();
                   }
                   if(hit_Cluster.first.subdetId()==EcalEndcap){ 
                      recoEnergy_shared+=hit_Cluster.second*(*recHitsEE->find(hit_Cluster.first)).energy(); 
                      recoEnergy_shared_noHitsFraction+=(*recHitsEE->find(hit_Cluster.first)).energy(); 
                   } 
                   if(std::find(superCluster_IDs.begin(),superCluster_IDs.end(),hit_Cluster.first) == superCluster_IDs.end()) superCluster_IDs.push_back(hit_Cluster.first);
                } 
            } 
        }
    }

    nSharedXtals = (int)superCluster_IDs.size();
    if(nSharedXtals<0.) nSharedXtals = -1.; 

    if(simEnergy>0.) simFraction_noHitsFraction = simEnergy_shared_noHitsFraction/simEnergy;
    else simFraction_noHitsFraction = -1.; 

    if(simEnergy>0.) simFraction = simEnergy_shared/simEnergy;
    else simFraction = -1.;
 
    if(simEnergy_shared>0.) recoToSim = recoEnergy/simEnergy_shared;
    else recoToSim = -1.; 

    if(simEnergy_shared>0.) recoToSim_shared = recoEnergy_shared/simEnergy_shared;
    else recoToSim_shared = -1.;  

    if(simEnergy_shared<0) simEnergy_shared = -1.;
    if(recoEnergy_shared<0) recoEnergy_shared = -1.;
    if(simEnergy_shared_noHitsFraction<0) simEnergy_shared_noHitsFraction = -1.;
    if(recoEnergy_shared_noHitsFraction<0) recoEnergy_shared_noHitsFraction = -1.;
    
    scores[0] = (double)nSharedXtals;
    scores[1] = simFraction_noHitsFraction;
    scores[2] = simFraction;
    scores[3] = recoToSim;
    scores[4] = recoToSim_shared;
    scores[5] = simEnergy_shared;
    scores[6] = recoEnergy_shared;
    scores[7] = simEnergy_shared_noHitsFraction;
    scores[8] = recoEnergy_shared_noHitsFraction;

    for(unsigned iVar=0; iVar<scores.size(); iVar++)
        if(std::isnan(scores.at(iVar))) std::cout << "score = " << iVar << " ---> NAN " << std::endl; 

    return scores;
}

std::vector<double> RecoSimDumper::getNoise(const reco::PFCluster* pfCluster, const std::vector<std::vector<std::pair<DetId, float> >> *hits_and_energies_CaloPart, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPartPU, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE, const EcalLaserAlphas* laserAlpha, const EcalLaserAPDPNRatios* laserRatio, const EcalIntercalibConstants* ical, const EcalIntercalibConstants* icalMC, const EcalPedestals* ped, const EcalADCToGeVConstant* adcToGeV, const EcalGainRatios* gr, bool useFractions=true)
{
    std::vector<double> noises;
    noises.resize(5);
    
    double recoE = 0.;
    double simE = 0.;
    double recoE_uncalib = 0.;
    double simE_uncalib = 0.;
    double noiseDB = 0.;
    double noiseDB_uncalib = 0.;
    
    const std::vector<std::pair<DetId,float> >* hitsAndFractions = &pfCluster->hitsAndFractions();  
    for(const std::pair<DetId, float>& hit_Cluster : *hitsAndFractions){ 

        double ic = *ical->getMap().find(hit_Cluster.first.rawId()); 
        double icMC = *icalMC->getMap().find(hit_Cluster.first.rawId());     
        double alpha = *laserAlpha->getMap().find(hit_Cluster.first.rawId());
        double apdpn = (*laserRatio->getLaserMap().find(hit_Cluster.first.rawId())).p2;
        double laserCorr = pow(apdpn,-alpha);
        double pedestal = 1.;
        double agv = 1.;
        double gain12Over6 = (*gr->getMap().find(hit_Cluster.first.rawId())).gain12Over6();
        double gain6Over1 = (*gr->getMap().find(hit_Cluster.first.rawId())).gain6Over1();
        if(hit_Cluster.first.subdetId()==EcalBarrel){ 
           agv = adcToGeV->getEBValue();
           pedestal = (*ped->find(hit_Cluster.first.rawId())).rms_x12;
           if((*recHitsEB->find(hit_Cluster.first)).checkFlag(16)) pedestal = (*ped->find(hit_Cluster.first.rawId())).rms_x6 * gain12Over6;
           if((*recHitsEB->find(hit_Cluster.first)).checkFlag(17)) pedestal = (*ped->find(hit_Cluster.first.rawId())).rms_x1 * gain6Over1;
        }else if(hit_Cluster.first.subdetId()==EcalEndcap){
           agv = adcToGeV->getEEValue();
           pedestal = (*ped->find(hit_Cluster.first.rawId())).rms_x12;
           if((*recHitsEE->find(hit_Cluster.first)).checkFlag(16)) pedestal = (*ped->find(hit_Cluster.first.rawId())).rms_x6 * gain12Over6;
           if((*recHitsEE->find(hit_Cluster.first)).checkFlag(17)) pedestal = (*ped->find(hit_Cluster.first.rawId())).rms_x1 * gain6Over1;
        }

        double fraction = 1.; 
        if(useFractions) fraction = hit_Cluster.second; 
        double icRatio = ic/icMC;

        noiseDB += pedestal*fraction*agv*laserCorr*ic;
        noiseDB_uncalib += pedestal*fraction*agv;

        if(hit_Cluster.first.subdetId()==EcalBarrel){ 
           recoE += fraction*(*recHitsEB->find(hit_Cluster.first)).energy();
           recoE_uncalib += fraction*(*recHitsEB->find(hit_Cluster.first)).energy()/(ic*laserCorr);
        }
        if(hit_Cluster.first.subdetId()==EcalEndcap){ 
           recoE += fraction*(*recHitsEE->find(hit_Cluster.first)).energy();
           recoE_uncalib += fraction*(*recHitsEE->find(hit_Cluster.first)).energy()/(ic*laserCorr);   
        }

        for(unsigned int iCalo=0; iCalo<hits_and_energies_CaloPart->size(); iCalo++){
            for(const std::pair<DetId, float>& hit_CaloPart : hits_and_energies_CaloPart->at(iCalo)){
                if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
                   simE+=icRatio*fraction*hit_CaloPart.second;
                   simE_uncalib+=icRatio*fraction*hit_CaloPart.second/(ic*laserCorr);
                } 
            } 
        }
        for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPartPU){
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
               simE+=icRatio*fraction*hit_CaloPart.second;
               simE_uncalib+=icRatio*fraction*hit_CaloPart.second/(ic*laserCorr);
            } 
        } 
    }  
    
    noises[0] = recoE_uncalib;
    noises[1] = recoE-simE;
    noises[2] = recoE_uncalib-simE_uncalib;
    noises[3] = noiseDB;
    noises[4] = noiseDB_uncalib;

    return noises;
}

int RecoSimDumper::getMatchedIndex(std::vector<std::vector<double>>* scores, double selection, bool useMax, double scale, int iCl)
{
   int matchedIndex = -1; 

   std::vector<double> score;
   for(unsigned int iCalo=0; iCalo<scores->at(iCl).size(); iCalo++){
       if(scores->at(iCl).at(iCalo)>0.) score.push_back(fabs(scale-scores->at(iCl).at(iCalo))); 
       else score.push_back(-1.); 
   }

   if(!useMax){ 
      std::replace(score.begin(),score.end(), -1., 999.);
      if(std::all_of(score.begin(),score.end(),[](double i){return i==-1.;}) || std::all_of(score.begin(),score.end(),[](double i){return i==999.;})) matchedIndex=-1;
      else matchedIndex = std::min_element(score.begin(),score.end()) - score.begin();
   }else{
      std::replace(score.begin(),score.end(), -1., -999.); 
      if(std::all_of(score.begin(),score.end(),[](double i){return i==-1.;}) || std::all_of(score.begin(),score.end(),[](double i){return i==-999.;})) matchedIndex=-1;
      else matchedIndex = std::max_element(score.begin(),score.end()) - score.begin();
   }

   if(matchedIndex==-1) return -1;
   
   if(useMax && score.at(matchedIndex) > selection) return matchedIndex;
   if(!useMax && score.at(matchedIndex) < selection) return matchedIndex; 
   
   return -1; 
}

void RecoSimDumper::fillParticleMatchedIndex(std::vector<std::vector<int>>* particleMatchedIndex, std::vector<int>* clusterMatchedIndex)
{
   for(unsigned int iCl=0; iCl<clusterMatchedIndex->size(); iCl++)
          if(clusterMatchedIndex->at(iCl)>=0) particleMatchedIndex->at(clusterMatchedIndex->at(iCl)).push_back(iCl);
}

GlobalPoint RecoSimDumper::calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES)
{
  double preshowerStartEta = 1.653;
  double preshowerEndEta = 2.6;
  double cl_energy_float = 0;
  double max_e = 0.0;
  double clusterT0 = 0.0;
  DetId id_max; 
 
  // find the seed and max layer
  for(const std::pair<DetId, float>& hit_CP : *hits_and_energies_CP){
    const double rh_energyf = hit_CP.second;
    cl_energy_float += rh_energyf;
    if (rh_energyf > max_e) {
      max_e = rh_energyf;
      id_max = hit_CP.first;
    }
  }

  const CaloSubdetectorGeometry* ecal_geom = nullptr;
  // get seed geometry information
  if(id_max.subdetId()==EcalBarrel){
      ecal_geom = _ebGeom;
      clusterT0 = _param_T0_EB;
  }else if(id_max.subdetId()==EcalEndcap){  
      ecal_geom = _eeGeom;
      clusterT0 = _param_T0_EE;
  }
  if (ecal_geom == nullptr)
     return GlobalPoint(-999999., -999999., -999999.);

  auto center_cell = ecal_geom->getGeometry(id_max);
  const double ctreta = center_cell->etaPos();
  const double actreta = std::abs(ctreta);
  // need to change T0 if in ES
  if (actreta > preshowerStartEta && actreta < preshowerEndEta && useES) {
    if (ctreta > 0 && _esPlus)
      clusterT0 = _param_T0_ES;
    if (ctreta < 0 && _esMinus)
      clusterT0 = _param_T0_ES;
  }

  // floats to reproduce exactly the EGM code
  const float maxDepth = _param_X0 * (clusterT0 + log(cl_energy_float));
  const float maxToFront = center_cell->getPosition().mag();
  // calculate the position
  const double logETot_inv = -log(cl_energy_float);
  double position_norm = 0.0;
  double x(0.0), y(0.0), z(0.0);
 
  for(const std::pair<DetId, float>& hit_CP : *hits_and_energies_CP) {
    if(hit_CP.first.subdetId()!=EcalBarrel && hit_CP.first.subdetId()!=EcalEndcap) continue;
    auto cell = ecal_geom->getGeometry(hit_CP.first);
    if(!cell.get()) continue;
    double weight = 0.0;
    const double rh_energy = hit_CP.second;
    if (rh_energy > 0.0)
      weight = std::max(0.0, (_param_W0 + log(rh_energy) + logETot_inv));
    const float depth = maxDepth + maxToFront - cell->getPosition().mag();
    const GlobalPoint pos = static_cast<const TruncatedPyramid*>(cell.get())->getPosition(depth);
    
    x += weight * pos.x();
    y += weight * pos.y();
    z += weight * pos.z();

    position_norm += weight;
  }

  // FALL BACK to LINEAR WEIGHTS
  if (position_norm == 0.) {
    for(const std::pair<DetId, float>& hit_CP : *hits_and_energies_CP) {
      if(hit_CP.first.subdetId()!=EcalBarrel && hit_CP.first.subdetId()!=EcalEndcap) continue; 
      auto cell = ecal_geom->getGeometry(hit_CP.first);
      if(!cell.get()) continue;
      double weight = 0.0;
      const double rh_energy = hit_CP.second;
      if (rh_energy > 0.0)
        weight = rh_energy / cl_energy_float;

      const float depth = maxDepth + maxToFront - cell->getPosition().mag();
      const GlobalPoint pos = cell->getPosition(depth);

      x += weight * pos.x();
      y += weight * pos.y();
      z += weight * pos.z();

      position_norm += weight;
    }
  }

  if (position_norm < _minAllowedNorm) {
    edm::LogError("WeirdClusterNormalization") << "Cluster too far from seeding cell: set position to (0,0,0).";
    return GlobalPoint(0., 0., 0.);
  } else {
    const double norm_inverse = 1.0 / position_norm;
    x *= norm_inverse;
    y *= norm_inverse;
    z *= norm_inverse;
  }

  return GlobalPoint(x, y, z); 
}

bool RecoSimDumper::passesJedID(const pat::Jet* jet, std::string version, std::string algo, std::string quality) 
{
    float eta      = jet->eta();
    float NHF      = jet->neutralHadronEnergyFraction();
    float NEMF     = jet->neutralEmEnergyFraction();
    float CHF      = jet->chargedHadronEnergyFraction();
    float MUF      = jet->muonEnergyFraction();
    float CEMF     = jet->chargedEmEnergyFraction();
    int   NumConst = jet->chargedMultiplicity()+jet->neutralMultiplicity();
    int   CHM      = jet->chargedMultiplicity();
    int   NumNeutralParticles = jet->neutralMultiplicity();

    jetIdMap.clear();

    //WPs definition here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID

    //LEGACY16 CHS/PUPPI LOOSE/TIGHT/TIGHTLEPVETO
    jetIdMap["2.4"]["LEGACY16"]["CHS"]["LOOSE"]              = ( abs(eta)<=2.4 && NHF<0.99 && NEMF<0.99 && NumConst>1 && CHF>0 && CHM>0 && CEMF<0.99 );
    jetIdMap["2.7"]["LEGACY16"]["CHS"]["LOOSE"]              = ( abs(eta)>2.4 && abs(eta)<=2.7 && NHF<0.99 && NEMF<0.99 && NumConst>1 ); 
    jetIdMap["3.0"]["LEGACY16"]["CHS"]["LOOSE"]              = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 ); 
    jetIdMap["forward"]["LEGACY16"]["CHS"]["LOOSE"]          = ( abs(eta)>3.0 && NEMF<0.90 && NumNeutralParticles>10 );  

    jetIdMap["2.4"]["LEGACY16"]["CHS"]["TIGHT"]              = ( abs(eta)<=2.4 && NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 && CEMF<0.99 ); 
    jetIdMap["2.7"]["LEGACY16"]["CHS"]["TIGHT"]              = ( abs(eta)>2.4 && abs(eta)<=2.7 && NHF<0.90 && NEMF<0.90 && NumConst>1 ); 
    jetIdMap["3.0"]["LEGACY16"]["CHS"]["TIGHT"]              = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 ); 
    jetIdMap["forward"]["LEGACY16"]["CHS"]["TIGHT"]          = ( abs(eta)>3.0 && NEMF<0.90 && NumNeutralParticles>10 );  

    jetIdMap["2.4"]["LEGACY16"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)<=2.4 && NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 && CHF>0 && CHM>0 && CEMF<0.90 ); 
    jetIdMap["2.7"]["LEGACY16"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)>2.4 && abs(eta)<=2.7 && NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 ); 
    jetIdMap["3.0"]["LEGACY16"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 ); 
    jetIdMap["forward"]["LEGACY16"]["CHS"]["TIGHTLEPVETO"]   = ( abs(eta)>3.0 && NEMF<0.90 && NumNeutralParticles>10 ); 

    jetIdMap["2.4"]["LEGACY16"]["PUPPI"]["LOOSE"]            = ( abs(eta)<=2.4 && NHF<0.99 && NEMF<0.99 && NumConst>1 && CHF>0 && CHM>0 && CEMF<0.99 );
    jetIdMap["2.7"]["LEGACY16"]["PUPPI"]["LOOSE"]            = ( abs(eta)>2.4 && abs(eta)<=2.7 && NHF<0.99 && NEMF<0.99 && NumConst>1 ); 
    jetIdMap["3.0"]["LEGACY16"]["PUPPI"]["LOOSE"]            = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 ); 
    jetIdMap["forward"]["LEGACY16"]["PUPPI"]["LOOSE"]        = ( abs(eta)>3.0 && NEMF<0.90 && NumNeutralParticles>10 );  

    jetIdMap["2.4"]["LEGACY16"]["PUPPI"]["TIGHT"]            = ( abs(eta)<=2.4 && NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 && CEMF<0.99 ); 
    jetIdMap["2.7"]["LEGACY16"]["PUPPI"]["TIGHT"]            = ( abs(eta)>2.4 && abs(eta)<=2.7 && NHF<0.90 && NEMF<0.90 && NumConst>1 ); 
    jetIdMap["3.0"]["LEGACY16"]["PUPPI"]["TIGHT"]            = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 ); 
    jetIdMap["forward"]["LEGACY16"]["PUPPI"]["TIGHT"]        = ( abs(eta)>3.0 && NEMF<0.90 && NumNeutralParticles>10 );  

    jetIdMap["2.4"]["LEGACY16"]["PUPPI"]["TIGHTLEPVETO"]     = ( abs(eta)<=2.4 && NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 && CHF>0 && CHM>0 && CEMF<0.90 ); 
    jetIdMap["2.7"]["LEGACY16"]["PUPPI"]["TIGHTLEPVETO"]     = ( abs(eta)>2.4 && abs(eta)<=2.7 && NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 ); 
    jetIdMap["3.0"]["LEGACY16"]["PUPPI"]["TIGHTLEPVETO"]     = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 ); 
    jetIdMap["forward"]["LEGACY16"]["PUPPI"]["TIGHTLEPVETO"] = ( abs(eta)>3.0 && NEMF<0.90 && NumNeutralParticles>10 ); 


    //RERECO17 CHS/PUPPI TIGHT/TIGHTLEPVETO
    jetIdMap["2.4"]["RERECO17"]["CHS"]["TIGHT"]              = ( abs(eta)<=2.4 && NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 ); 
    jetIdMap["2.7"]["RERECO17"]["CHS"]["TIGHT"]              = ( abs(eta)>2.4 && abs(eta)<=2.7 && NHF<0.90 && NEMF<0.90 && NumConst>1 ); 
    jetIdMap["3.0"]["RERECO17"]["CHS"]["TIGHT"]              = ( abs(eta)>2.7 && abs(eta)<=3.0 && NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 );
    jetIdMap["forward"]["RERECO17"]["CHS"]["TIGHT"]          = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10 ); 

    jetIdMap["2.4"]["RERECO17"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)<=2.4 && NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 && CHF>0 && CHM>0 && CEMF<0.8 ); 
    jetIdMap["2.7"]["RERECO17"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)>2.4 && abs(eta)<=2.7 && NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 ); 
    jetIdMap["3.0"]["RERECO17"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)>2.7 && abs(eta)<=3.0 && NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 );
    jetIdMap["forward"]["RERECO17"]["CHS"]["TIGHTLEPVETO"]   = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10 ); 

    jetIdMap["2.4"]["RERECO17"]["PUPPI"]["TIGHT"]            = ( abs(eta)<=2.4 && NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 ); 
    jetIdMap["2.7"]["RERECO17"]["PUPPI"]["TIGHT"]            = ( abs(eta)>2.4 && abs(eta)<=2.7 && NHF<0.90 && NEMF<0.90 && NumConst>1 ); 
    jetIdMap["3.0"]["RERECO17"]["PUPPI"]["TIGHT"]            = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<0.99 );
    jetIdMap["forward"]["RERECO17"]["PUPPI"]["TIGHT"]        = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticles>2 && NumNeutralParticles<15 ); 

    jetIdMap["2.4"]["RERECO17"]["PUPPI"]["TIGHTLEPVETO"]     = ( abs(eta)<=2.4 && NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 && CHF>0 && CHM>0 && CEMF<0.8 ); 
    jetIdMap["2.7"]["RERECO17"]["PUPPI"]["TIGHTLEPVETO"]     = ( abs(eta)>2.4 && abs(eta)<=2.7 && NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 ); 
    jetIdMap["3.0"]["RERECO17"]["PUPPI"]["TIGHTLEPVETO"]     = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<0.99 );
    jetIdMap["forward"]["RERECO17"]["PUPPI"]["TIGHTLEPVETO"] = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticles>2 && NumNeutralParticles<15 );

    
    //RERECO18 CHS/PUPPI TIGHT/TIGHTLEPVETO
    jetIdMap["2.6"]["RERECO18"]["CHS"]["TIGHT"]              = ( abs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 ); 
    jetIdMap["2.7"]["RERECO18"]["CHS"]["TIGHT"]              = ( abs(eta)>2.6 && abs(eta)<=2.7 && CHM>0 && NEMF<0.99 && NHF < 0.9 ); 
    jetIdMap["3.0"]["RERECO18"]["CHS"]["TIGHT"]              = ( abs(eta)>2.7 && abs(eta)<=3.0 && NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 ); 
    jetIdMap["forward"]["RERECO18"]["CHS"]["TIGHT"]          = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 ); 

    jetIdMap["2.6"]["RERECO18"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["2.7"]["RERECO18"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)>2.6 && abs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["3.0"]["RERECO18"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)>2.7 && abs(eta)<=3.0 && NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 ); 
    jetIdMap["forward"]["RERECO18"]["CHS"]["TIGHTLEPVETO"]   = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 ); 

    jetIdMap["2.6"]["RERECO18"]["PUPPI"]["TIGHT"]              = ( abs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 ); 
    jetIdMap["2.7"]["RERECO18"]["PUPPI"]["TIGHT"]              = ( abs(eta)>2.6 && abs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 ); 
    jetIdMap["3.0"]["RERECO18"]["PUPPI"]["TIGHT"]              = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<0.99 ); 
    jetIdMap["forward"]["RERECO18"]["PUPPI"]["TIGHT"]          = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticles>2 && NumNeutralParticles<15 ); 

    jetIdMap["2.6"]["RERECO18"]["PUPPI"]["TIGHTLEPVETO"]       = ( abs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["2.7"]["RERECO18"]["PUPPI"]["TIGHTLEPVETO"]       = ( abs(eta)>2.6 && abs(eta)<=2.7 && CEMF<0.8 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["3.0"]["RERECO18"]["PUPPI"]["TIGHTLEPVETO"]       = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<0.99 ); 
    jetIdMap["forward"]["RERECO18"]["PUPPI"]["TIGHTLEPVETO"]   = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticles>2 && NumNeutralParticles<15 );  


    //UL16 CHS/PUPPI TIGHT/TIGHTLEPVETO
    jetIdMap["2.4"]["UL16"]["CHS"]["TIGHT"]                  = ( abs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );
    jetIdMap["2.7"]["UL16"]["CHS"]["TIGHT"]                  = ( abs(eta)>2.4 && abs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 ); 
    jetIdMap["3.0"]["UL16"]["CHS"]["TIGHT"]                  = ( abs(eta)>2.7 && abs(eta)<=3.0 && NEMF>0.0 && NEMF<0.99 && NHF<0.9 && NumNeutralParticles>1 ); 
    jetIdMap["forward"]["UL16"]["CHS"]["TIGHT"]              = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 ); 

    jetIdMap["2.4"]["UL16"]["CHS"]["TIGHTLEPVETO"]           = ( abs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 );
    jetIdMap["2.7"]["UL16"]["CHS"]["TIGHTLEPVETO"]           = ( abs(eta)>2.4 && abs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 );
    jetIdMap["3.0"]["UL16"]["CHS"]["TIGHTLEPVETO"]           = ( abs(eta)>2.7 && abs(eta)<=3.0 && NEMF>0.0 && NEMF<0.99 && NHF<0.9 && NumNeutralParticles>1 ); 
    jetIdMap["forward"]["UL16"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 );

    jetIdMap["2.4"]["UL16"]["PUPPI"]["TIGHT"]                = ( abs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 ); 
    jetIdMap["2.7"]["UL16"]["PUPPI"]["TIGHT"]                = ( abs(eta)>2.4 && abs(eta)<=2.7 && NEMF<0.99 && NHF < 0.98 ); 
    jetIdMap["3.0"]["UL16"]["PUPPI"]["TIGHT"]                = ( abs(eta)>2.7 && abs(eta)<=3.0 && NumNeutralParticles>=1 ); 
    jetIdMap["forward"]["UL16"]["PUPPI"]["TIGHT"]            = ( abs(eta)>3.0 && NEMF<0.90 && NumNeutralParticles>2 );  

    jetIdMap["2.4"]["UL16"]["PUPPI"]["TIGHTLEPVETO"]         = ( abs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["2.7"]["UL16"]["PUPPI"]["TIGHTLEPVETO"]         = ( abs(eta)>2.4 && abs(eta)<=2.7 && NEMF<0.99 && NHF < 0.98 ); 
    jetIdMap["3.0"]["UL16"]["PUPPI"]["TIGHTLEPVETO"]         = ( abs(eta)>2.7 && abs(eta)<=3.0 && NumNeutralParticles>=1 ); 
    jetIdMap["forward"]["UL16"]["PUPPI"]["TIGHTLEPVETO"]     = ( abs(eta)>3.0 && NEMF<0.90 && NumNeutralParticles>2 );  

    
    //UL17 CHS/PUPPI TIGHT/TIGHTLEPVETO
    jetIdMap["2.6"]["UL17"]["CHS"]["TIGHT"]                  = ( abs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 ); 
    jetIdMap["2.7"]["UL17"]["CHS"]["TIGHT"]                  = ( abs(eta)>2.6 && abs(eta)<=2.7 && CHM>0 && NEMF<0.99 && NHF < 0.9 ); 
    jetIdMap["3.0"]["UL17"]["CHS"]["TIGHT"]                  = ( abs(eta)>2.7 && abs(eta)<=3.0 && NEMF>0.01 && NEMF<0.99 && NumNeutralParticles>1 );  
    jetIdMap["forward"]["UL17"]["CHS"]["TIGHT"]              = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 );

    jetIdMap["2.6"]["UL17"]["CHS"]["TIGHTLEPVETO"]           = ( abs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["2.7"]["UL17"]["CHS"]["TIGHTLEPVETO"]           = ( abs(eta)>2.6 && abs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["3.0"]["UL17"]["CHS"]["TIGHTLEPVETO"]           = ( abs(eta)>2.7 && abs(eta)<=3.0 && NEMF>0.01 && NEMF<0.99 && NumNeutralParticles>1 );   
    jetIdMap["forward"]["UL17"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 );

    jetIdMap["2.6"]["UL17"]["PUPPI"]["TIGHT"]                = ( abs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );  
    jetIdMap["2.7"]["UL17"]["PUPPI"]["TIGHT"]                = ( abs(eta)>2.6 && abs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 ); 
    jetIdMap["3.0"]["UL17"]["PUPPI"]["TIGHT"]                = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<1.0 ); 
    jetIdMap["forward"]["UL17"]["PUPPI"]["TIGHT"]            = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticles>2 ); 

    jetIdMap["2.6"]["UL17"]["PUPPI"]["TIGHTLEPVETO"]         = ( abs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["2.7"]["UL17"]["PUPPI"]["TIGHTLEPVETO"]         = ( abs(eta)>2.6 && abs(eta)<=2.7 && CEMF<0.8 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["3.0"]["UL17"]["PUPPI"]["TIGHTLEPVETO"]         = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<1.0 ); 
    jetIdMap["forward"]["UL17"]["PUPPI"]["TIGHTLEPVETO"]     = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticles>2 );  

    
    //UL18 CHS/PUPPI TIGHT/TIGHTLEPVETO
    jetIdMap["2.6"]["UL18"]["CHS"]["TIGHT"]                  = ( abs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 ); 
    jetIdMap["2.7"]["UL18"]["CHS"]["TIGHT"]                  = ( abs(eta)>2.6 && abs(eta)<=2.7 && CHM>0 && NEMF<0.99 && NHF < 0.9 ); 
    jetIdMap["3.0"]["UL18"]["CHS"]["TIGHT"]                  = ( abs(eta)>2.7 && abs(eta)<=3.0 && NEMF>0.01 && NEMF<0.99 && NumNeutralParticles>1 );  
    jetIdMap["forward"]["UL18"]["CHS"]["TIGHT"]              = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 );

    jetIdMap["2.6"]["UL18"]["CHS"]["TIGHTLEPVETO"]           = ( abs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["2.7"]["UL18"]["CHS"]["TIGHTLEPVETO"]           = ( abs(eta)>2.6 && abs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["3.0"]["UL18"]["CHS"]["TIGHTLEPVETO"]           = ( abs(eta)>2.7 && abs(eta)<=3.0 && NEMF>0.01 && NEMF<0.99 && NumNeutralParticles>1 );   
    jetIdMap["forward"]["UL18"]["CHS"]["TIGHTLEPVETO"]       = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 );

    jetIdMap["2.6"]["UL18"]["PUPPI"]["TIGHT"]                = ( abs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 );  
    jetIdMap["2.7"]["UL18"]["PUPPI"]["TIGHT"]                = ( abs(eta)>2.6 && abs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 ); 
    jetIdMap["3.0"]["UL18"]["PUPPI"]["TIGHT"]                = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<1.0 ); 
    jetIdMap["forward"]["UL18"]["PUPPI"]["TIGHT"]            = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticles>2 ); 

    jetIdMap["2.6"]["UL18"]["PUPPI"]["TIGHTLEPVETO"]         = ( abs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["2.7"]["UL18"]["PUPPI"]["TIGHTLEPVETO"]         = ( abs(eta)>2.6 && abs(eta)<=2.7 && CEMF<0.8 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 ); 
    jetIdMap["3.0"]["UL18"]["PUPPI"]["TIGHTLEPVETO"]         = ( abs(eta)>2.7 && abs(eta)<=3.0 && NHF<1.0 ); 
    jetIdMap["forward"]["UL18"]["PUPPI"]["TIGHTLEPVETO"]     = ( abs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticles>2 );  


    if(version=="LEGACY16" || version=="RERECO17" || version=="UL16"){
       if(fabs(eta)<=2.4) return jetIdMap["2.4"][version][algo][quality];
       if(fabs(eta)<=2.7) return jetIdMap["2.7"][version][algo][quality];
       if(fabs(eta)<=3.0) return jetIdMap["3.0"][version][algo][quality]; 
       if(fabs(eta)>3.0)  return jetIdMap["forward"][version][algo][quality]; 
    }else if(version=="RERECO18" || version=="UL18" || version=="UL17"){ 
       if(fabs(eta)<=2.6) return jetIdMap["2.6"][version][algo][quality];
       if(fabs(eta)<=2.7) return jetIdMap["2.7"][version][algo][quality];
       if(fabs(eta)<=3.0) return jetIdMap["3.0"][version][algo][quality]; 
       if(fabs(eta)>3.0)  return jetIdMap["forward"][version][algo][quality];  
    }
    return false;
    
}

float RecoSimDumper::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoSimDumper);

