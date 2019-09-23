#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/LooperFactory.h"
#include "FWCore/Framework/interface/ESProducerLooper.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
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
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "RecoSimStudies/Dumpers/plugins/RecoSimDumper.h"

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
RecoSimDumper::RecoSimDumper(const edm::ParameterSet& iConfig)
{

   genToken_                = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   caloPartToken_           = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
   PCaloHitEBToken_         = consumes< std::vector<PCaloHit> >(iConfig.getParameter<edm::InputTag>("PCaloHitEBCollection"));
   PCaloHitEEToken_         = consumes< std::vector<PCaloHit> >(iConfig.getParameter<edm::InputTag>("PCaloHitEECollection"));
   ebRechitToken_           = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_           = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
   pfRecHitToken_           = consumes<std::vector<reco::PFRecHit> >(iConfig.getParameter<edm::InputTag>("pfRechitCollection")); 
   pfClusterToken_          = consumes<std::vector<reco::PFCluster> >(iConfig.getParameter<edm::InputTag>("pfClusterCollection")); 
   ebSuperClusterToken_     = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebSuperClusterCollection"));
   eeSuperClusterToken_     = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeSuperClusterCollection"));
   rhoToken_                = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoTag"));

   doCompression_           = iConfig.getParameter<bool>("doCompression");
   nBits_                   = iConfig.getParameter<int>("nBits");
   saveCalohits_            = iConfig.getParameter<bool>("saveCalohits");
   saveSimhits_             = iConfig.getParameter<bool>("saveSimhits");
   saveRechits_             = iConfig.getParameter<bool>("saveRechits");
   savePFRechits_           = iConfig.getParameter<bool>("savePFRechits"); 
   savePFCluster_           = iConfig.getParameter<bool>("savePFCluster");
   saveSuperCluster_        = iConfig.getParameter<bool>("saveSuperCluster");
   useEnergyRegression_     = iConfig.getParameter<bool>("useEnergyRegression"); 
   genID_                   = iConfig.getParameter<std::vector<int>>("genID");

   if(nBits_>23 && doCompression_){
      cout << "WARNING: float compression bits > 23 ---> Using 23 (i.e. no compression) instead!" << endl;
      nBits_=23;
   }

   //output file, historgrams and trees

   //gInterpreter->GenerateDictionary("vector<vector<bool>>","vector;vector");
   //gInterpreter->GenerateDictionary("vector<map<int,int>>","map;vector"); 

   tree = iFile->make<TTree>("caloTree","caloTree"); 
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/I");
   tree->Branch("runId", &runId, "runId/I");
   tree->Branch("rho", &rho, "rho/D");
   tree->Branch("genParticle_id","std::vector<int>",&genParticle_id);
   tree->Branch("genParticle_energy","std::vector<float>",&genParticle_energy);
   tree->Branch("genParticle_pt","std::vector<float>",&genParticle_pt);
   tree->Branch("genParticle_eta","std::vector<float>",&genParticle_eta);
   tree->Branch("genParticle_phi","std::vector<float>",&genParticle_phi);
   tree->Branch("caloParticle_energy","std::vector<float>",&caloParticle_energy);
   tree->Branch("caloParticle_simEnergy","std::vector<float>",&caloParticle_simEnergy); 
   tree->Branch("caloParticle_pt","std::vector<float>",&caloParticle_pt);
   tree->Branch("caloParticle_eta","std::vector<float>",&caloParticle_eta);
   tree->Branch("caloParticle_phi","std::vector<float>",&caloParticle_phi);
   tree->Branch("caloParticle_ieta","std::vector<int>",&caloParticle_ieta);
   tree->Branch("caloParticle_iphi","std::vector<int>",&caloParticle_iphi);
   tree->Branch("caloParticle_iz","std::vector<int>",&caloParticle_iz);
   if(saveCalohits_){
      tree->Branch("caloHit_energy","std::vector<std::vector<float> >",&caloHit_energy);
      tree->Branch("caloHit_time","std::vector<std::vector<float> >",&caloHit_time);
      tree->Branch("caloHit_eta","std::vector<std::vector<float> >",&caloHit_eta);
      tree->Branch("caloHit_phi","std::vector<std::vector<float> >",&caloHit_phi);
      tree->Branch("caloHit_ieta","std::vector<std::vector<int> >",&caloHit_ieta);
      tree->Branch("caloHit_iphi","std::vector<std::vector<int> >",&caloHit_iphi);
      tree->Branch("caloHit_iz","std::vector<std::vector<int> >",&caloHit_iz);
   }
   if(saveSimhits_){
      tree->Branch("simHit_energy","std::vector<std::vector<float> >",&simHit_energy);
      tree->Branch("simHit_eta","std::vector<std::vector<float> >",&simHit_eta);
      tree->Branch("simHit_phi","std::vector<std::vector<float> >",&simHit_phi);
      tree->Branch("simHit_ieta","std::vector<std::vector<int> >",&simHit_ieta);
      tree->Branch("simHit_iphi","std::vector<std::vector<int> >",&simHit_iphi);
      tree->Branch("simHit_iz","std::vector<std::vector<int> >",&simHit_iz);
   }
   if(saveRechits_){
      tree->Branch("recHit_energy","std::vector<std::vector<float> >",&recHit_energy);
      if(!saveSimhits_){  
         tree->Branch("recHit_eta","std::vector<std::vector<float> >",&recHit_eta); 
         tree->Branch("recHit_phi","std::vector<std::vector<float> >",&recHit_phi);
         tree->Branch("recHit_ieta","std::vector<std::vector<int> >",&recHit_ieta); 
         tree->Branch("recHit_iphi","std::vector<std::vector<int> >",&recHit_iphi);
         tree->Branch("recHit_iz","std::vector<std::vector<int> >",&recHit_iz);  
      }   
   }
   if(savePFRechits_){ 
      tree->Branch("pfRecHit_unMatched_eta","std::vector<float>",&pfRecHit_unMatched_eta); 
      tree->Branch("pfRecHit_unMatched_phi","std::vector<float>",&pfRecHit_unMatched_phi); 
      tree->Branch("pfRecHit_unMatched_ieta","std::vector<int>",&pfRecHit_unMatched_ieta); 
      tree->Branch("pfRecHit_unMatched_iphi","std::vector<int>",&pfRecHit_unMatched_iphi);
      tree->Branch("pfRecHit_unMatched_iz","std::vector<int>",&pfRecHit_unMatched_iz);  
      if(saveRechits_) tree->Branch("pfRecHit_isMatched","std::vector<std::vector<bool> >",&pfRecHit_isMatched);
      if(!saveRechits_) tree->Branch("pfRecHit_energy","std::vector<std::vector<float> >",&pfRecHit_energy);  
      if(!saveSimhits_ && !saveRechits_){
         tree->Branch("pfRecHit_eta","std::vector<std::vector<float> >",&pfRecHit_eta); 
         tree->Branch("pfRecHit_phi","std::vector<std::vector<float> >",&pfRecHit_phi); 
         tree->Branch("pfRecHit_ieta","std::vector<std::vector<int> >",&pfRecHit_ieta); 
         tree->Branch("pfRecHit_iphi","std::vector<std::vector<int> >",&pfRecHit_iphi);
         tree->Branch("pfRecHit_iz","std::vector<std::vector<int> >",&pfRecHit_iz);  
      }  
   }
   if(savePFCluster_){
      tree->Branch("pfClusterHit_energy","std::vector<std::vector<std::map<int, float> >",&pfClusterHit_energy);
      if(!saveSimhits_ && !saveRechits_ && !savePFRechits_){ 
         tree->Branch("pfClusterHit_eta","std::vector<std::vector<float> >",&pfClusterHit_eta);
         tree->Branch("pfClusterHit_phi","std::vector<std::vector<float> >",&pfClusterHit_phi);
         tree->Branch("pfClusterHit_ieta","std::vector<std::vector<int> >",&pfClusterHit_ieta);
         tree->Branch("pfClusterHit_iphi","std::vector<std::vector<int> >",&pfClusterHit_iphi);
         tree->Branch("pfClusterHit_iz","std::vector<std::vector<int> >",&pfClusterHit_iz);
      }
      tree->Branch("pfClusterHit_noCaloPart_energy","std::vector<std::vector<float> >",&pfClusterHit_noCaloPart_energy);
      tree->Branch("pfClusterHit_noCaloPart_eta","std::vector<std::vector<float> >",&pfClusterHit_noCaloPart_eta);
      tree->Branch("pfClusterHit_noCaloPart_phi","std::vector<std::vector<float> >",&pfClusterHit_noCaloPart_phi);
      tree->Branch("pfClusterHit_noCaloPart_ieta","std::vector<std::vector<int> >",&pfClusterHit_noCaloPart_ieta);
      tree->Branch("pfClusterHit_noCaloPart_iphi","std::vector<std::vector<int> >",&pfClusterHit_noCaloPart_iphi);
      tree->Branch("pfClusterHit_noCaloPart_iz","std::vector<std::vector<int> >",&pfClusterHit_noCaloPart_iz);
      tree->Branch("pfCluster_energy","std::vector<float>",&pfCluster_energy);
      tree->Branch("pfCluster_eta","std::vector<float>",&pfCluster_eta);
      tree->Branch("pfCluster_phi","std::vector<float>",&pfCluster_phi);   
      tree->Branch("pfCluster_ieta","std::vector<int>",&pfCluster_ieta);
      tree->Branch("pfCluster_iphi","std::vector<int>",&pfCluster_iphi);   
      tree->Branch("pfCluster_iz","std::vector<int>",&pfCluster_iz);   
   } 
   if(saveSuperCluster_){
      tree->Branch("superClusterHit_energy","std::vector<std::vector<std::map<int, float> >",&superClusterHit_energy);
      if(!saveSimhits_ && !saveRechits_ && !savePFRechits_ && !savePFCluster_){ 
         tree->Branch("superClusterHit_eta","std::vector<std::vector<float> >",&superClusterHit_eta);
         tree->Branch("superClusterHit_phi","std::vector<std::vector<float> >",&superClusterHit_phi);
         tree->Branch("superClusterHit_ieta","std::vector<std::vector<int> >",&superClusterHit_ieta);
         tree->Branch("superClusterHit_iphi","std::vector<std::vector<int> >",&superClusterHit_iphi);
         tree->Branch("superClusterHit_iz","std::vector<std::vector<int> >",&superClusterHit_iz);
      }
      tree->Branch("superClusterHit_noCaloPart_energy","std::vector<std::vector<float> >",&superClusterHit_noCaloPart_energy);
      tree->Branch("superClusterHit_noCaloPart_eta","std::vector<std::vector<float> >",&superClusterHit_noCaloPart_eta);
      tree->Branch("superClusterHit_noCaloPart_phi","std::vector<std::vector<float> >",&superClusterHit_noCaloPart_phi);
      tree->Branch("superClusterHit_noCaloPart_ieta","std::vector<std::vector<int> >",&superClusterHit_noCaloPart_ieta);
      tree->Branch("superClusterHit_noCaloPart_iphi","std::vector<std::vector<int> >",&superClusterHit_noCaloPart_iphi);
      tree->Branch("superClusterHit_noCaloPart_iz","std::vector<std::vector<int> >",&superClusterHit_noCaloPart_iz);
      tree->Branch("superCluster_energy","std::vector<float> ",&superCluster_energy);
      tree->Branch("superCluster_eta","std::vector<float>",&superCluster_eta);
      tree->Branch("superCluster_phi","std::vector<float>",&superCluster_phi);  
      tree->Branch("superCluster_ieta","std::vector<int>",&superCluster_ieta);
      tree->Branch("superCluster_iphi","std::vector<int>",&superCluster_iphi);  
      tree->Branch("superCluster_iz","std::vector<int>",&superCluster_iz);   
   }
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
   edm::ESHandle<CaloGeometry> caloGeometry;
   iSetup.get<CaloGeometryRecord>().get(caloGeometry);
   const CaloGeometry *geometry = caloGeometry.product();
   _ebGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
   _eeGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
   _esGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
   if (_esGeom) {
    for (uint32_t ic = 0; ic < _esGeom->getValidDetIds().size() && (!_esPlus || !_esMinus); ++ic) {
      const double z = _esGeom->getGeometry(_esGeom->getValidDetIds()[ic])->getPosition().z();
      _esPlus = _esPlus || (0 < z);
      _esMinus = _esMinus || (0 > z);
    }
   }
   
   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   ev.getByToken(genToken_,genParticles);
   if (!genParticles.isValid()) {
       std::cerr << "Analyze --> genParticles not found" << std::endl; 
       return;
   }

   edm::Handle<std::vector<CaloParticle> > caloParticles;
   ev.getByToken(caloPartToken_,caloParticles);
   if (!caloParticles.isValid()) {
       std::cerr << "Analyze --> caloParticles not found" << std::endl; 
       return;
   }

   edm::Handle<std::vector<PCaloHit> > PCaloHitsEB;
   ev.getByToken(PCaloHitEBToken_, PCaloHitsEB);
   if(saveCalohits_){
      if (!PCaloHitsEB.isValid()) {
          std::cerr << "Analyze --> PCaloHitsEB not found" << std::endl; 
          return;
      }
   }

   edm::Handle<std::vector<PCaloHit> > PCaloHitsEE;
   ev.getByToken(PCaloHitEEToken_, PCaloHitsEE);
   if(saveCalohits_){
      if (!PCaloHitsEE.isValid()) {
          std::cerr << "Analyze --> PCaloHitsEE not found" << std::endl; 
          return;
      }
   }

   edm::Handle<EcalRecHitCollection> recHitsEB;
   ev.getByToken(ebRechitToken_, recHitsEB);
   if(saveRechits_) {
      if (!recHitsEB.isValid()) {
          std::cerr << "Analyze --> recHitsEB not found" << std::endl; 
          return;
      }
   }

   edm::Handle<EcalRecHitCollection> recHitsEE;
   ev.getByToken(eeRechitToken_, recHitsEE);
   if(saveRechits_) {
      if (!recHitsEE.isValid()) {
          std::cerr << "Analyze --> recHitsEE not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::PFRecHit> > pfRecHits;
   ev.getByToken(pfRecHitToken_, pfRecHits);
   if(savePFRechits_) {
      if (!pfRecHits.isValid()) {
          std::cerr << "Analyze --> pfRecHits not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::PFCluster> > pfClusters;
   ev.getByToken(pfClusterToken_, pfClusters);
   if(savePFCluster_) {
      if (!pfClusters.isValid()) {
          std::cerr << "Analyze --> pfClusters not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEB;
   ev.getByToken(ebSuperClusterToken_, superClusterEB);
   if(saveSuperCluster_) {
      if (!superClusterEB.isValid()) {
          std::cerr << "Analyze --> superClusterEB not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEE;
   ev.getByToken(eeSuperClusterToken_, superClusterEE);
   if(saveSuperCluster_) {
      if (!superClusterEE.isValid()) {
          std::cerr << "Analyze --> superClusterEE not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<double> rhoHandle;
   ev.getByToken(rhoToken_, rhoHandle);
   if (! rhoHandle.isValid()) {
      std::cerr << "Analyze --> rho handle not found" << std::endl; 
      return;
   }

   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();
   rho = *(rhoHandle.product());

   caloParticleXtals_.clear();
   caloParticleXtals_ = caloParticleXtals(caloParticles,&genID_);
   int nCaloParticles = caloParticleXtals_.size();
  
   genParticle_id.clear();
   genParticle_energy.clear();
   genParticle_pt.clear();
   genParticle_eta.clear();
   genParticle_phi.clear();
   caloParticle_energy.clear();
   caloParticle_simEnergy.clear();
   caloParticle_pt.clear();
   caloParticle_eta.clear();
   caloParticle_phi.clear();
   caloParticle_ieta.clear();
   caloParticle_iphi.clear();
   caloParticle_iz.clear();
 
   caloHit_energy.clear();
   caloHit_time.clear();
   caloHit_eta.clear();
   caloHit_phi.clear();
   caloHit_ieta.clear();
   caloHit_iphi.clear();
   caloHit_iz.clear();  
   caloHit_energy.resize(nCaloParticles);
   caloHit_time.resize(nCaloParticles);
   caloHit_eta.resize(nCaloParticles);
   caloHit_phi.resize(nCaloParticles);
   caloHit_ieta.resize(nCaloParticles);
   caloHit_iphi.resize(nCaloParticles);
   caloHit_iz.resize(nCaloParticles);
 
   simHit_energy.clear();
   simHit_eta.clear();
   simHit_phi.clear();
   simHit_ieta.clear();
   simHit_iphi.clear();
   simHit_iz.clear();
   simHit_energy.resize(nCaloParticles);
   simHit_eta.resize(nCaloParticles);
   simHit_phi.resize(nCaloParticles);
   simHit_ieta.resize(nCaloParticles);
   simHit_iphi.resize(nCaloParticles);
   simHit_iz.resize(nCaloParticles);

   recHit_energy.clear();
   recHit_eta.clear();
   recHit_phi.clear();
   recHit_ieta.clear();
   recHit_iphi.clear();
   recHit_iz.clear();  
   recHit_energy.resize(nCaloParticles);
   recHit_eta.resize(nCaloParticles);
   recHit_phi.resize(nCaloParticles);
   recHit_ieta.resize(nCaloParticles);
   recHit_iphi.resize(nCaloParticles);
   recHit_iz.resize(nCaloParticles);  

   pfRecHit_isMatched.clear();
   pfRecHit_energy.clear();
   pfRecHit_eta.clear();
   pfRecHit_phi.clear();
   pfRecHit_ieta.clear();
   pfRecHit_iphi.clear();
   pfRecHit_iz.clear();
   pfRecHit_isMatched.resize(nCaloParticles);  
   pfRecHit_energy.resize(nCaloParticles);  
   pfRecHit_eta.resize(nCaloParticles);  
   pfRecHit_phi.resize(nCaloParticles);  
   pfRecHit_ieta.resize(nCaloParticles);
   pfRecHit_iphi.resize(nCaloParticles);
   pfRecHit_iz.resize(nCaloParticles);

   pfRecHit_unMatched_energy.clear();
   pfRecHit_unMatched_eta.clear();
   pfRecHit_unMatched_phi.clear();
   pfRecHit_unMatched_ieta.clear();
   pfRecHit_unMatched_iphi.clear();
   pfRecHit_unMatched_iz.clear();
   
   pfClusterHit_energy.clear();
   pfClusterHit_eta.clear();
   pfClusterHit_phi.clear();   
   pfClusterHit_ieta.clear();
   pfClusterHit_iphi.clear(); 
   pfClusterHit_iz.clear();           
   pfClusterHit_energy.resize(nCaloParticles);  
   pfClusterHit_eta.resize(nCaloParticles);  
   pfClusterHit_phi.resize(nCaloParticles);     
   pfClusterHit_ieta.resize(nCaloParticles);  
   pfClusterHit_iphi.resize(nCaloParticles);   
   pfClusterHit_iz.resize(nCaloParticles);   

   int nPFClusters = (pfClusters.product())->size();
   pfClusterHit_noCaloPart_energy.clear();
   pfClusterHit_noCaloPart_eta.clear();
   pfClusterHit_noCaloPart_phi.clear();   
   pfClusterHit_noCaloPart_ieta.clear();
   pfClusterHit_noCaloPart_iphi.clear(); 
   pfClusterHit_noCaloPart_iz.clear();           
   pfClusterHit_noCaloPart_energy.resize(nPFClusters);  
   pfClusterHit_noCaloPart_eta.resize(nPFClusters);  
   pfClusterHit_noCaloPart_phi.resize(nPFClusters);     
   pfClusterHit_noCaloPart_ieta.resize(nPFClusters);  
   pfClusterHit_noCaloPart_iphi.resize(nPFClusters);   
   pfClusterHit_noCaloPart_iz.resize(nPFClusters);   
        
   pfCluster_energy.clear();
   pfCluster_eta.clear();
   pfCluster_phi.clear();
   pfCluster_ieta.clear();
   pfCluster_iphi.clear();
   pfCluster_iz.clear();
   
   superClusterHit_energy.clear();
   superClusterHit_eta.clear();
   superClusterHit_phi.clear();
   superClusterHit_ieta.clear();
   superClusterHit_iphi.clear();
   superClusterHit_iz.clear();
   superClusterHit_energy.resize(nCaloParticles);
   superClusterHit_eta.resize(nCaloParticles);
   superClusterHit_phi.resize(nCaloParticles);
   superClusterHit_ieta.resize(nCaloParticles);
   superClusterHit_iphi.resize(nCaloParticles);
   superClusterHit_iz.resize(nCaloParticles);

   int nSuperClusters = (superClusterEB.product())->size() + (superClusterEE.product())->size();
   superClusterHit_noCaloPart_energy.clear();
   superClusterHit_noCaloPart_eta.clear();
   superClusterHit_noCaloPart_phi.clear();
   superClusterHit_noCaloPart_ieta.clear();
   superClusterHit_noCaloPart_iphi.clear();
   superClusterHit_noCaloPart_iz.clear();
   superClusterHit_noCaloPart_energy.resize(nSuperClusters);
   superClusterHit_noCaloPart_eta.resize(nSuperClusters);
   superClusterHit_noCaloPart_phi.resize(nSuperClusters);
   superClusterHit_noCaloPart_ieta.resize(nSuperClusters);
   superClusterHit_noCaloPart_iphi.resize(nSuperClusters);
   superClusterHit_noCaloPart_iz.resize(nSuperClusters);

   superCluster_energy.clear(); 
   superCluster_eta.clear(); 
   superCluster_phi.clear();  
   superCluster_ieta.clear(); 
   superCluster_iphi.clear();    
   superCluster_iz.clear();   
    
   GlobalPoint cell;

   std::cout << "CaloParticles size  : " << nCaloParticles << std::endl;
   for(const auto& iCalo : *(caloParticles.product()))
   {
       bool isGoodParticle = false; 
       for(unsigned int id=0; id<genID_.size(); id++) 
           //M.G. 
           if(iCalo.pdgId()==genID_.at(id) || genID_.at(id)==0) isGoodParticle=true;
           //iisGoodParticle=true;

       if(!isGoodParticle) continue;     

       const auto& genParticles_caloPart = iCalo.genParticles();
       genParticle_id.push_back(iCalo.pdgId());
       if(genParticles_caloPart.empty()){
          cout << "WARNING: no associated genParticle found, making standard dR matching" << endl;
          float dR=999.;
          int igen_tmp=-1; 
          int igen=0; 
          for(const auto& iGen : *(genParticles.product()))
          {
              float dR_tmp = deltaR(iCalo.eta(),iCalo.phi(),iGen.eta(),iGen.phi());  
              if(dR_tmp<dR && iGen.status()==1){
                 dR=dR_tmp;
                 igen_tmp=igen;
              }  
              igen++;
          } 
          const auto& genParticles_tmp = *(genParticles.product());
          auto genParticle = genParticles_tmp[igen_tmp]; 
          genParticle_energy.push_back(reduceFloat(genParticle.energy(),nBits_));
          genParticle_pt.push_back(reduceFloat(genParticle.pt(),nBits_));
          genParticle_eta.push_back(reduceFloat(genParticle.eta(),nBits_));
          genParticle_phi.push_back(reduceFloat(genParticle.phi(),nBits_));
       }else{
          genParticle_energy.push_back(reduceFloat((*genParticles_caloPart.begin())->energy(),nBits_));
          genParticle_pt.push_back(reduceFloat((*genParticles_caloPart.begin())->pt(),nBits_));
          genParticle_eta.push_back(reduceFloat((*genParticles_caloPart.begin())->eta(),nBits_));
          genParticle_phi.push_back(reduceFloat((*genParticles_caloPart.begin())->phi(),nBits_));
       }
 
       caloParticle_energy.push_back(reduceFloat(iCalo.energy(),nBits_));
       caloParticle_pt.push_back(reduceFloat(iCalo.pt(),nBits_));
      
       std::vector<std::pair<DetId, float> >* hitsAndEnergies_CP = getHitsAndEnergiesCaloPart(iCalo);
       GlobalPoint caloParticle_position = calculateAndSetPositionActual(hitsAndEnergies_CP, 7.4, 3.1, 1.2, 4.2, 0.89, 0., geometry, false);
       caloParticle_eta.push_back(reduceFloat(caloParticle_position.eta(),nBits_));
       caloParticle_phi.push_back(reduceFloat(caloParticle_position.phi(),nBits_));
       if(std::abs(caloParticle_position.eta()) < 1.479){  
          EBDetId eb_id(_ebGeom->getClosestCell(caloParticle_position));  
          caloParticle_ieta.push_back(eb_id.ieta());
          caloParticle_iphi.push_back(eb_id.iphi());
          caloParticle_iz.push_back(0); 
       }else{            
          int iz=-99;
          EEDetId ee_id(_eeGeom->getClosestCell(caloParticle_position));   
          caloParticle_ieta.push_back(ee_id.ix());
          caloParticle_iphi.push_back(ee_id.iy());
          if(ee_id.zside()<0) iz=-1;
          if(ee_id.zside()>0) iz=1;  
          caloParticle_iz.push_back(iz); 
       }   
   }

   for(unsigned int iCaloCount=0; iCaloCount<caloParticleXtals_.size(); iCaloCount++) 
   {
       //Get hits from caloParticle, and associated recHits, pfRechits, PFClusterhit and superClusterhit 

       // This indexes are not used but could be useful
       int simHit_index=-1; 
       int recHit_index=-1; 
       int pfRecHit_index=-1; 
       int superClusterHit_index=-1;

       float calo_simEnergy=0.;
       for(auto const& hit: caloParticleXtals_[iCaloCount])
       {
           DetId id(hit.first);
           if(id.subdetId()!=EcalBarrel && id.subdetId()!=EcalEndcap) continue;
               
           calo_simEnergy += hit.second; 

          // PfCluster and superClusters associated with this detid
           map<int, float> map_pfCluster_energy;
           map<int, float> map_superCluster_energy;

           cell = geometry->getPosition(id);
           float eta = cell.eta();  
           float phi = cell.phi();  
           int ieta = -99; 
           int iphi = -99;
           int iz = -99;  
           if(id.subdetId()==EcalBarrel){
              EBDetId eb_id(id);
              ieta = eb_id.ieta(); 
              iphi = eb_id.iphi();  
              iz = 0;   
           }
           if(id.subdetId()==EcalEndcap){
              EEDetId ee_id(id);
              ieta = ee_id.ix(); 
              iphi = ee_id.iy();
              if(ee_id.zside()<0) iz=-1;
              if(ee_id.zside()>0) iz=1;   
           }

           if(saveSimhits_){
              simHit_index++;
              simHit_energy[iCaloCount].push_back(reduceFloat(hit.second,nBits_));
              simHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_));
              simHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_));
              simHit_ieta[iCaloCount].push_back(ieta);
              simHit_iphi[iCaloCount].push_back(iphi);
              simHit_iz[iCaloCount].push_back(iz); 
           }

           float recHit_energy_ = -1.;
           if(id.subdetId()==EcalBarrel){

              //Save associated caloHit energy 
              if(saveCalohits_){
                 for(auto &ipCaloHit : *(PCaloHitsEB.product())){
                     if(ipCaloHit.id() == id.rawId()){
                        caloHit_energy[iCaloCount].push_back(reduceFloat(ipCaloHit.energy(),nBits_));
                        caloHit_time[iCaloCount].push_back(reduceFloat(ipCaloHit.time(),nBits_));
                        caloHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_));
                        caloHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_)); 
                        caloHit_ieta[iCaloCount].push_back(ieta);
                        caloHit_iphi[iCaloCount].push_back(iphi);          
                     }
                 } 
              }

              //Save associated recHit energy
              for(auto &iRecHit : *(recHitsEB.product())){
                  if(iRecHit.id().rawId() == id.rawId()){
                     recHit_index++;  
                     recHit_energy_ = iRecHit.energy();;
                     if(!saveSimhits_ && saveRechits_){
                        recHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_));
                        recHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_));   
                        recHit_ieta[iCaloCount].push_back(ieta);
                        recHit_iphi[iCaloCount].push_back(iphi);
                        recHit_iz[iCaloCount].push_back(0); 
                     } 
                     break;
                  }   
              }    
              if(saveRechits_){
                 recHit_energy[iCaloCount].push_back(reduceFloat(recHit_energy_,nBits_));
              }
   
              //Save SuperClusterHit energy
              float superClusterHit_energy_ = -1.;
              int superCluster_index_tmp=0;
              if(saveSuperCluster_){
                 for(const auto& iSuperCluster : *(superClusterEB.product())){
                     for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
                         const std::vector<std::pair<DetId,float> > &seedrechits = ( *iBC )->hitsAndFractions();
                         for(unsigned int i = 0; i < seedrechits.size(); i++){      
                             if(seedrechits[i].first.rawId() == id.rawId()){   
                                //for matched SuperClusterHit   
                                if(useEnergyRegression_) superClusterHit_energy_ = recHit_energy_*seedrechits[i].second;
                                else superClusterHit_energy_ = (iSuperCluster.rawEnergy()/iSuperCluster.energy())*recHit_energy_*seedrechits[i].second;
                                map_superCluster_energy.insert(pair<int,float>(superCluster_index_tmp,reduceFloat(superClusterHit_energy_,nBits_) ));
                                break;
                             }                   
                         }
                     } 
                     superCluster_index_tmp++;
                 }
                 superClusterHit_index++; 
                 superClusterHit_energy[iCaloCount].push_back(map_superCluster_energy);
                 if(!saveSimhits_ && !saveRechits_ && !savePFRechits_ && !savePFCluster_){  
                    superClusterHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_)); 
                    superClusterHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_)); 
                    superClusterHit_ieta[iCaloCount].push_back(ieta);
                    superClusterHit_iphi[iCaloCount].push_back(iphi);
                    superClusterHit_iz[iCaloCount].push_back(iz); 
                 }      
              }

           }else if(id.subdetId()==EcalEndcap){

              //Save associated caloHit energy 
              if(saveCalohits_){
                 for(auto &ipCaloHit : *(PCaloHitsEE.product())){
                     if(ipCaloHit.id() == id.rawId()){
                        caloHit_energy[iCaloCount].push_back(reduceFloat(ipCaloHit.energy(),nBits_));
                        caloHit_time[iCaloCount].push_back(reduceFloat(ipCaloHit.time(),nBits_));
                        caloHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_));
                        caloHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_)); 
                        caloHit_ieta[iCaloCount].push_back(ieta);
                        caloHit_iphi[iCaloCount].push_back(iphi);   
                        caloHit_iz[iCaloCount].push_back(iz);   
                     }
                 } 
              }
                  
              //Save associated recHit energy
              for(auto &iRecHit : *(recHitsEE.product())){
                  if(iRecHit.id().rawId() == id.rawId()){
                     recHit_index++;   
                     recHit_energy_ = iRecHit.energy();;
                     if(!saveSimhits_ && saveRechits_){
                        recHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_));
                        recHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_));   
                        recHit_ieta[iCaloCount].push_back(ieta);
                        recHit_iphi[iCaloCount].push_back(iphi);
                        recHit_iz[iCaloCount].push_back(0); 
                     } 
                     break;
                  }
              }     
              if(saveRechits_){
                 recHit_energy[iCaloCount].push_back(reduceFloat(recHit_energy_,nBits_));
              }  
                   
              //Save SuperClusterHit energy
              float superClusterHit_energy_ = -1.;
              int superCluster_index_tmp=0; 
              if(saveSuperCluster_){
                 for(const auto& iSuperCluster : *(superClusterEE.product())){
                     for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
                         const std::vector<std::pair<DetId,float> > &seedrechits = ( *iBC )->hitsAndFractions();
                         for(unsigned int i = 0; i < seedrechits.size(); i++){      
                             if(seedrechits[i].first.rawId() == id.rawId()){      
                                if(useEnergyRegression_) superClusterHit_energy_ = recHit_energy_*seedrechits[i].second;
                                else superClusterHit_energy_ = (iSuperCluster.rawEnergy()/iSuperCluster.energy())*recHit_energy_*seedrechits[i].second;
                                map_superCluster_energy.insert(pair<int,float>(superCluster_index_tmp,reduceFloat(superClusterHit_energy_,nBits_) ));
                                break;
                             }                 
                         }
                     }    
                     superCluster_index_tmp++;
                 }
                 superClusterHit_index++; 
                 superClusterHit_energy[iCaloCount].push_back(map_superCluster_energy);
                 if(!saveSimhits_ && !saveRechits_ && !savePFRechits_ && !savePFCluster_){  
                    superClusterHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_)); 
                    superClusterHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_)); 
                    superClusterHit_ieta[iCaloCount].push_back(ieta);
                    superClusterHit_iphi[iCaloCount].push_back(iphi);
                    superClusterHit_iz[iCaloCount].push_back(iz); 
                 }      
              }
           }
                
           //Save associated pfRechit energy
           bool pfRecHit_isMatched_ = false;
           float pfRecHit_energy_ = -1.; 
           if(savePFRechits_){
              for(const auto& iPFRechit : *(pfRecHits.product())){
                  if(iPFRechit.detId() == id.rawId()){
                     pfRecHit_index++;  
                     pfRecHit_isMatched_ = true;
                     pfRecHit_energy_ = iPFRechit.energy();
                     break;
                  }   
              }  
              if(saveRechits_ && savePFRechits_) pfRecHit_isMatched[iCaloCount].push_back(pfRecHit_isMatched_);
              if(!saveRechits_ && savePFRechits_) pfRecHit_energy[iCaloCount].push_back(reduceFloat(pfRecHit_energy_,nBits_));
              if(!saveSimhits_ && !saveRechits_){  
                 pfRecHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_));
                 pfRecHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_)); 
                 pfRecHit_ieta[iCaloCount].push_back(ieta);
                 pfRecHit_iphi[iCaloCount].push_back(iphi);
                 pfRecHit_iz[iCaloCount].push_back(iz);  
              }   
           } 
                
           //Save PFClusterHit energy
           float pfClusterHit_energy_ = -1.;
           int pfCluster_index_tmp = 0;
           if(savePFCluster_){                         
              for(const auto& iPFCluster : *(pfClusters.product())){
                  reco::CaloCluster caloBC(iPFCluster);
                  const std::vector<std::pair<DetId,float> > &hitsAndFractions = caloBC.hitsAndFractions();
                  for(unsigned int i = 0; i < hitsAndFractions.size(); i++){
                      if(hitsAndFractions[i].first.rawId() == id.rawId()){       
                         if(!useEnergyRegression_) pfClusterHit_energy_ = recHit_energy_*hitsAndFractions[i].second;
                         else pfClusterHit_energy_ = (iPFCluster.correctedEnergy()/iPFCluster.energy())*recHit_energy_*hitsAndFractions[i].second;
                         // Save clusterHit with the cluster id
                         map_pfCluster_energy.insert(pair<int, float>(pfCluster_index_tmp,reduceFloat(pfClusterHit_energy_,nBits_)));
                         break; 
                      }    
                  }  
                  pfCluster_index_tmp++;
              }   
              pfClusterHit_energy[iCaloCount].push_back(map_pfCluster_energy);
              if(!saveSimhits_ && !saveRechits_ && !savePFRechits_){  
                 pfClusterHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_));
                 pfClusterHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_)); 
                 pfClusterHit_ieta[iCaloCount].push_back(ieta);
                 pfClusterHit_iphi[iCaloCount].push_back(iphi);
                 pfClusterHit_iz[iCaloCount].push_back(iz);  
              } 
           } 
                  
       } // --> End of loop on simhits
       
       caloParticle_simEnergy.push_back(reduceFloat(calo_simEnergy,nBits_));
   }  // --> End of loop on caloparticles
   
   //Save PFClusters
   if(savePFCluster_){
      int iPFCl=0;
      std::cout << "PFClusters size     : " << (pfClusters.product())->size() << std::endl;
      for(const auto& iPFCluster : *(pfClusters.product())){      
          pfCluster_energy.push_back(reduceFloat(iPFCluster.energy(),nBits_));
          pfCluster_eta.push_back(reduceFloat(iPFCluster.eta(),nBits_));
          pfCluster_phi.push_back(reduceFloat(iPFCluster.phi(),nBits_));
          reco::CaloCluster caloBC(iPFCluster);
          math::XYZPoint caloPos = caloBC.position();
          if(std::abs(iPFCluster.eta()) < 1.479){  
             EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
             pfCluster_ieta.push_back(eb_id.ieta());
             pfCluster_iphi.push_back(eb_id.iphi());
             pfCluster_iz.push_back(0); 
          }else{  
             int iz=-99;
             EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
             if(ee_id.zside()<0) iz=-1;
             if(ee_id.zside()>0) iz=1;     
             pfCluster_ieta.push_back(ee_id.ix());
             pfCluster_iphi.push_back(ee_id.iy());
             pfCluster_iz.push_back(iz); 
          }   
          
          //for unmatched PFClusterHit      
          const std::vector<std::pair<DetId,float> > &hitsAndFractions = caloBC.hitsAndFractions();
          for(unsigned int i = 0; i < hitsAndFractions.size(); i++){

              bool isMatched = false;
              for(unsigned int iCalo=0; iCalo<caloParticleXtals_.size(); iCalo++){  
                  std::map<uint32_t,float>::iterator it;
                  it = caloParticleXtals_[iCalo].find(hitsAndFractions[i].first.rawId());                                   
                  if(it != caloParticleXtals_[iCalo].end()) isMatched = true;
              }
              if(isMatched == true) continue;  

              float clusterHit_noCaloPart_energy_ = -1.;
              float energy_tmp_ = -1.;
              if(hitsAndFractions[i].first.subdetId()==EcalBarrel) energy_tmp_ = (*(recHitsEB.product())->find(hitsAndFractions[i].first)).energy();
              else if(hitsAndFractions[i].first.subdetId()==EcalEndcap) energy_tmp_ = (*(recHitsEE.product())->find(hitsAndFractions[i].first)).energy();  
 
              if(useEnergyRegression_) clusterHit_noCaloPart_energy_ = energy_tmp_*hitsAndFractions[i].second;
              else clusterHit_noCaloPart_energy_ = (iPFCluster.correctedEnergy()/iPFCluster.energy())*energy_tmp_*hitsAndFractions[i].second;

              cell = geometry->getPosition(hitsAndFractions[i].first);
              pfClusterHit_noCaloPart_energy[iPFCl].push_back(reduceFloat(clusterHit_noCaloPart_energy_,nBits_));
              pfClusterHit_noCaloPart_eta[iPFCl].push_back(reduceFloat(cell.eta(),nBits_));
              pfClusterHit_noCaloPart_phi[iPFCl].push_back(reduceFloat(cell.phi(),nBits_));
              if(hitsAndFractions[i].first.subdetId()==EcalBarrel){ 
                 EBDetId eb_id(hitsAndFractions[i].first);  
                 pfClusterHit_noCaloPart_ieta[iPFCl].push_back(eb_id.ieta());
                 pfClusterHit_noCaloPart_iphi[iPFCl].push_back(eb_id.iphi());
                 pfClusterHit_noCaloPart_iz[iPFCl].push_back(0); 
              }else if(hitsAndFractions[i].first.subdetId()==EcalEndcap){  
                 int iz=-99;
                 EEDetId ee_id(hitsAndFractions[i].first);  
                 pfClusterHit_noCaloPart_ieta[iPFCl].push_back(ee_id.ix());
                 pfClusterHit_noCaloPart_iphi[iPFCl].push_back(ee_id.iy());
                 if(ee_id.zside()<0) iz=-1;
                 if(ee_id.zside()>0) iz=1;   
                 pfClusterHit_noCaloPart_iz[iPFCl].push_back(iz); 
              } 
          }
          iPFCl++;        
      } 
   }

   //Save SuperClusters 
   if(saveSuperCluster_){
      int iSC=0;
      std::cout << "SuperClustersEB size: " << (superClusterEB.product())->size() << std::endl;
      for(const auto& iSuperCluster : *(superClusterEB.product())){    
          superCluster_energy.push_back(reduceFloat(iSuperCluster.energy(),nBits_));
          superCluster_eta.push_back(reduceFloat(iSuperCluster.eta(),nBits_));
          superCluster_phi.push_back(reduceFloat(iSuperCluster.phi(),nBits_));
          math::XYZPoint caloPos = iSuperCluster.seed()->position();
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          superCluster_ieta.push_back(eb_id.ieta());
          superCluster_iphi.push_back(eb_id.iphi());
          superCluster_iz.push_back(0);   
          //for unmatched SuperClusterHit 
          for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
              const std::vector<std::pair<DetId,float> > &seedrechits = ( *iBC )->hitsAndFractions();
              for(unsigned int i = 0; i < seedrechits.size(); i++){    
    
                   bool isMatched = false;
                   for(unsigned int iCalo=0; iCalo<caloParticleXtals_.size(); iCalo++){  
                      std::map<uint32_t,float>::iterator it;
                      it = caloParticleXtals_[iCalo].find(seedrechits[i].first.rawId());                                   
                      if(it != caloParticleXtals_[iCalo].end()) isMatched = true;
                   }
                   if(isMatched == true) continue;  
                   
                   float superClusterHit_noCaloPart_energy_ = -1.;
                   float energy_tmp_ = (*(recHitsEB.product())->find(seedrechits[i].first)).energy();
                   if(useEnergyRegression_) superClusterHit_noCaloPart_energy_ = energy_tmp_*seedrechits[i].second;
                   else superClusterHit_noCaloPart_energy_ = (iSuperCluster.rawEnergy()/iSuperCluster.energy())*energy_tmp_*seedrechits[i].second;

                   cell = geometry->getPosition(seedrechits[i].first);
                   EBDetId eb_id(seedrechits[i].first);  
                   superClusterHit_noCaloPart_energy[iSC].push_back(reduceFloat(superClusterHit_noCaloPart_energy_,nBits_));
                   superClusterHit_noCaloPart_eta[iSC].push_back(reduceFloat(cell.eta(),nBits_));
                   superClusterHit_noCaloPart_phi[iSC].push_back(reduceFloat(cell.phi(),nBits_));
                   superClusterHit_noCaloPart_ieta[iSC].push_back(eb_id.ieta());
                   superClusterHit_noCaloPart_iphi[iSC].push_back(eb_id.iphi());
                   superClusterHit_noCaloPart_iz[iSC].push_back(0); 
                                     
              } 
          }      
          iSC++;  
      } 
      iSC=0;
      std::cout << "SuperClustersEE size: " << (superClusterEE.product())->size() << std::endl;
      for(const auto& iSuperCluster : *(superClusterEE.product())){    
          superCluster_energy.push_back(reduceFloat(iSuperCluster.energy(),nBits_));
          superCluster_eta.push_back(reduceFloat(iSuperCluster.eta(),nBits_));
          superCluster_phi.push_back(reduceFloat(iSuperCluster.phi(),nBits_));
          math::XYZPoint caloPos = iSuperCluster.seed()->position(); 
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          superCluster_ieta.push_back(ee_id.ix());
          superCluster_iphi.push_back(ee_id.iy());
          superCluster_iz.push_back(ee_id.zside());       
          //for unmatched SuperClusterHit 
          for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
              const std::vector<std::pair<DetId,float> > &seedrechits = ( *iBC )->hitsAndFractions();
              for(unsigned int i = 0; i < seedrechits.size(); i++){    

                   bool isMatched = false;
                   for(unsigned int iCalo=0; iCalo<caloParticleXtals_.size(); iCalo++){  
                      std::map<uint32_t,float>::iterator it;
                      it = caloParticleXtals_[iCalo].find(seedrechits[i].first.rawId());                                   
                      if(it != caloParticleXtals_[iCalo].end()) isMatched = true;
                   }
                   if(isMatched == true) continue;     
                   
                   float superClusterHit_noCaloPart_energy_ = -1.;
                   float energy_tmp_ = (*(recHitsEE.product())->find(seedrechits[i].first)).energy();
                   if(useEnergyRegression_) superClusterHit_noCaloPart_energy_ = energy_tmp_*seedrechits[i].second;
                   else superClusterHit_noCaloPart_energy_ = (iSuperCluster.rawEnergy()/iSuperCluster.energy())*energy_tmp_*seedrechits[i].second;
                      
                   int iz=-99;
                   cell = geometry->getPosition(seedrechits[i].first);
                   EEDetId ee_id(seedrechits[i].first);  
                   superClusterHit_noCaloPart_energy[iSC].push_back(reduceFloat(superClusterHit_noCaloPart_energy_,nBits_));
                   superClusterHit_noCaloPart_eta[iSC].push_back(reduceFloat(cell.eta(),nBits_));
                   superClusterHit_noCaloPart_phi[iSC].push_back(reduceFloat(cell.phi(),nBits_));
                   superClusterHit_noCaloPart_ieta[iSC].push_back(ee_id.ix());
                   superClusterHit_noCaloPart_iphi[iSC].push_back(ee_id.iy());
                   if(ee_id.zside()<0) iz=-1;
                   if(ee_id.zside()>0) iz=1; 
                   superClusterHit_noCaloPart_iz[iSC].push_back(iz); 
                                     
              } 
          }      
          iSC++;     
      } 
   }

   //Save unMatched pfRechits 
   if(savePFRechits_){
      for(const auto& iPFRechit : *(pfRecHits.product())){

          DetId pf_id(iPFRechit.detId());
          bool pfRecHit_isMatched_ = false;

          for(unsigned int iCaloCount=0; iCaloCount<caloParticleXtals_.size(); iCaloCount++) 
          {
              for(auto const& hit: caloParticleXtals_[iCaloCount])
              {
                  DetId id(hit.first);
                  if(iPFRechit.detId() == id.rawId()) pfRecHit_isMatched_ = true;
                  break;
              }
          }

          for(const auto& iPFCluster : *(pfClusters.product())){ 
              reco::CaloCluster caloBC(iPFCluster);     
              const std::vector<std::pair<DetId,float> > &hitsAndFractions = caloBC.hitsAndFractions();
              for(unsigned int i = 0; i < hitsAndFractions.size(); i++){
                  if(iPFRechit.detId() == hitsAndFractions[i].first.rawId()) pfRecHit_isMatched_ = true;
                  break;   
              }
          }

          if(pf_id.subdetId()==EcalBarrel){
             for(const auto& iSuperCluster : *(superClusterEB.product())){ 
                 for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
                     const std::vector<std::pair<DetId,float> > &seedrechits = ( *iBC )->hitsAndFractions();
                     for(unsigned int i = 0; i < seedrechits.size(); i++){  
                         if(iPFRechit.detId() == seedrechits[i].first.rawId()) pfRecHit_isMatched_ = true;
                         break;  
                     }
                 }
             }       
          }else if(pf_id.subdetId()==EcalEndcap){
             for(const auto& iSuperCluster : *(superClusterEE.product())){ 
                 for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
                     const std::vector<std::pair<DetId,float> > &seedrechits = ( *iBC )->hitsAndFractions();
                     for(unsigned int i = 0; i < seedrechits.size(); i++){  
                         if(iPFRechit.detId() == seedrechits[i].first.rawId()) pfRecHit_isMatched_ = true;
                         break;  
                     }
                 }
             }     
          }

          if(pfRecHit_isMatched_) continue;

          cell = geometry->getPosition(pf_id); 
          pfRecHit_unMatched_energy.push_back(iPFRechit.energy());    
          pfRecHit_unMatched_eta.push_back(cell.eta());  
          pfRecHit_unMatched_phi.push_back(cell.phi()); 
          if(pf_id.subdetId()==EcalBarrel){ 
             EBDetId eb_id(pf_id);  
             pfRecHit_unMatched_ieta.push_back(eb_id.ieta());  
             pfRecHit_unMatched_iphi.push_back(eb_id.iphi());  
             pfRecHit_unMatched_iz.push_back(0);     
          }else if(pf_id.subdetId()==EcalEndcap){
             int iz=-99;
             EEDetId ee_id(pf_id);  
             if(ee_id.zside()<0) iz=-1;
             if(ee_id.zside()>0) iz=1; 
             pfRecHit_unMatched_ieta.push_back(ee_id.ix());  
             pfRecHit_unMatched_iphi.push_back(ee_id.iy());  
             pfRecHit_unMatched_iz.push_back(iz);    
          } 
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
std::vector<std::map<uint32_t,float> > RecoSimDumper::caloParticleXtals(edm::Handle<std::vector<CaloParticle> > caloParticles, std::vector<int>* genID_)
{
    std::vector<std::map<uint32_t,float> > xtals;
    for(const auto& iCalo : *(caloParticles.product()))
    {
       bool isGoodParticle = false; 
       for(unsigned int id=0; id<genID_->size(); id++)
           //std::cout << "i=" << id << " genID=" <<genID_->at(id) << " iCalo.pdgId()=" << iCalo.pdgId() << std::endl;
           // M.G.
           if(iCalo.pdgId()==genID_->at(id) || genID_->at(id)==0) isGoodParticle=true;
      
       if(!isGoodParticle) continue;

       std::map<uint32_t,float> xtals_tmp;
       const auto& simClusters = iCalo.simClusters();
       for( unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
            auto simCluster = simClusters[iSC];  
            auto hitsAndEnergies = simCluster->hits_and_energies(); 
            for(unsigned int iHit= 0; iHit<hitsAndEnergies.size(); iHit++)
                xtals_tmp[hitsAndEnergies[iHit].first]+=hitsAndEnergies[iHit].second;
       }  

       xtals.push_back(xtals_tmp);                
    } 

    return xtals;
}

std::vector<std::pair<DetId, float> >* RecoSimDumper::getHitsAndEnergiesCaloPart(CaloParticle iCaloParticle)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_CP = new std::vector<std::pair<DetId, float> >;
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_tmp = new std::vector<std::pair<DetId, float> >;
    std::map<DetId, float> HitsAndEnergies_map;
    
    const auto& simClusters = iCaloParticle.simClusters();
    for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
        auto simCluster = simClusters[iSC];  
        auto hits_and_energies = simCluster->hits_and_energies();
        for(unsigned int i = 0; i < hits_and_energies.size(); i++){  
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
         HitsAndEnergies_CP->push_back(make_pair(hit.first,hit.second));

    return HitsAndEnergies_CP;
}

GlobalPoint RecoSimDumper::calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, const CaloGeometry *geometry, bool useES)
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
  }else{
      throw cms::Exception("InvalidLayer") << "ECAL Position Calc only accepts ECAL_BARREL or ECAL_ENDCAP";
  }

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

float RecoSimDumper::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoSimDumper);

