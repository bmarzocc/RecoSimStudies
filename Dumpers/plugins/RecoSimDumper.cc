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
      tree->Branch("pfClusterHit_energy","std::vector<std::vector<float> >",&pfClusterHit_energy);
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
      if(saveSimhits_)tree->Branch("map_simHit_pfCluster","std::vector<std::map<int,int> >",&map_simHit_pfCluster); 
      if(!saveSimhits_ && saveRechits_)tree->Branch("map_recHit_pfCluster","std::vector<std::map<int,int> >",&map_recHit_pfCluster); 
      if(!saveSimhits_ && !saveRechits_ && savePFRechits_)tree->Branch("map_pfRecHit_pfCluster","std::vector<std::map<int,int> >",&map_pfRecHit_pfCluster);    
      if(!saveSimhits_ && !saveRechits_ && !savePFRechits_)tree->Branch("map_pfClusterHit_pfCluster","std::vector<std::map<int,int> >",&map_pfClusterHit_pfCluster);        
   } 
   if(saveSuperCluster_){
      tree->Branch("superClusterHit_energy","std::vector<std::vector<float> >",&superClusterHit_energy);
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
      if(saveSimhits_)tree->Branch("map_simHit_superCluster","std::vector<std::map<int,int> >",&map_simHit_superCluster); 
      if(!saveSimhits_ && saveRechits_)tree->Branch("map_recHit_superCluster","std::vector<std::map<int,int> >",&map_recHit_superCluster); 
      if(!saveSimhits_ && !saveRechits_ && savePFRechits_)tree->Branch("map_pfRecHit_superCluster","std::vector<std::map<int,int> >",&map_pfRecHit_pfCluster);    
      if(!saveSimhits_ && !saveRechits_ && !savePFRechits_)tree->Branch("map_superClusterHit_superCluster","std::vector<std::map<int,int> >",&map_superClusterHit_superCluster);       
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

   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();
   //std::cout << "Event=" << eventId << " run=" << runId << std::endl;
   int nCaloParticles = nSkimmedCaloParticles(caloParticles,&genID_);
  
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
   
   map_simHit_pfCluster.clear(); 
   map_recHit_pfCluster.clear();
   map_pfRecHit_pfCluster.clear();
   map_pfClusterHit_pfCluster.clear();   
   map_simHit_superCluster.clear(); 
   map_recHit_superCluster.clear();
   map_pfRecHit_superCluster.clear();
   map_superClusterHit_superCluster.clear();    
   map_simHit_pfCluster.resize(nCaloParticles);
   map_recHit_pfCluster.resize(nCaloParticles);
   map_pfRecHit_pfCluster.resize(nCaloParticles);
   map_pfClusterHit_pfCluster.resize(nCaloParticles);   
   map_simHit_superCluster.resize(nCaloParticles); 
   map_recHit_superCluster.resize(nCaloParticles);
   map_pfRecHit_superCluster.resize(nCaloParticles);
   map_superClusterHit_superCluster.resize(nCaloParticles);    
  
   unmatchedHit_energy_Cl_.clear();
   unmatchedHit_energy_SCEB_.clear();
   unmatchedHit_energy_SCEE_.clear();

   int iCaloCount = 0;
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
       caloParticle_eta.push_back(reduceFloat(iCalo.eta(),nBits_));
       caloParticle_phi.push_back(reduceFloat(iCalo.phi(),nBits_));

       float calo_simEnergy=-1.;

       //Get hits from simClusters, and associated recHits, pfRechits, PFClusterhit and superClusterhit 
       int simHit_index=-1; 
       int recHit_index=-1; 
       int pfRecHit_index=-1; 
       int pfClusterHit_index=-1;
       int superClusterHit_index=-1;
       const auto& simClusters = iCalo.simClusters();
       for( unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
            auto simCluster = simClusters[iSC];  
            auto hits_and_fractions = simCluster->hits_and_fractions();
            auto hits_and_energies = simCluster->hits_and_energies();
            for(unsigned int iHit= 0; iHit<hits_and_fractions.size(); iHit++){
                DetId id(hits_and_fractions[iHit].first);
                if(id.subdetId()!=EcalBarrel && id.subdetId()!=EcalEndcap) continue;
                
                int pfCluster_index=-1;
                int superCluster_index=-1;

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

                calo_simEnergy += hits_and_energies[iHit].second; 
                if(saveSimhits_){
                   simHit_index++;
                   simHit_energy[iCaloCount].push_back(reduceFloat(hits_and_energies[iHit].second,nBits_));
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
                       if(iRecHit.id() == id.rawId()){
                          recHit_index++;  
                          recHit_energy_ = iRecHit.energy();;
                          if(saveRechits_){
                             recHit_energy[iCaloCount].push_back(reduceFloat(recHit_energy_,nBits_));
                             if(!saveSimhits_){
                                recHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_));
                                recHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_));   
                                recHit_ieta[iCaloCount].push_back(ieta);
                                recHit_iphi[iCaloCount].push_back(iphi);
                                recHit_iz[iCaloCount].push_back(0); 
                             } 
                          }
                       }       
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
                                     superCluster_index = superCluster_index_tmp;
                                  }else{
                                     //for unmatched SuperClusterHit 
                                     float superClusterHit_noCaloPart_energy_ = -1.;
                                     float energy_tmp_ = (*(recHitsEB.product())->find(seedrechits[i].first)).energy();
                                     if(useEnergyRegression_) superClusterHit_noCaloPart_energy_ = energy_tmp_*seedrechits[i].second;
                                     else superClusterHit_noCaloPart_energy_ = (iSuperCluster.rawEnergy()/iSuperCluster.energy())*energy_tmp_*seedrechits[i].second;
                                     unmatchedHit_energy_SCEB_[superCluster_index_tmp][seedrechits[i].first] = superClusterHit_noCaloPart_energy_;
                                  }                   
                              }
                          } 
                          superCluster_index_tmp++;
                      }
                      superClusterHit_index++; 
                      superClusterHit_energy[iCaloCount].push_back(reduceFloat(superClusterHit_energy_,nBits_));
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
                       if(iRecHit.id() == id.rawId()){
                          recHit_index++;   
                          recHit_energy_ = iRecHit.energy();;
                          if(saveRechits_){
                             recHit_energy[iCaloCount].push_back(reduceFloat(recHit_energy_,nBits_));
                             if(!saveSimhits_){
                                recHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_));
                                recHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_));   
                                recHit_ieta[iCaloCount].push_back(ieta);
                                recHit_iphi[iCaloCount].push_back(iphi);
                                recHit_iz[iCaloCount].push_back(0); 
                             } 
                          }
                       }       
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
                                     superCluster_index = superCluster_index_tmp;
                                     break;
                                  }else{
                                     //for unmatched SuperClusterHit 
                                     float superClusterHit_noCaloPart_energy_ = -1.;
                                     float energy_tmp_ = (*(recHitsEE.product())->find(seedrechits[i].first)).energy();
                                     if(useEnergyRegression_) superClusterHit_noCaloPart_energy_ = energy_tmp_*seedrechits[i].second;
                                     else superClusterHit_noCaloPart_energy_ = (iSuperCluster.rawEnergy()/iSuperCluster.energy())*energy_tmp_*seedrechits[i].second;
                                     unmatchedHit_energy_SCEE_[superCluster_index_tmp][seedrechits[i].first] = superClusterHit_noCaloPart_energy_;
                                  }                        
                              }
                          }    
                          superCluster_index_tmp++;
                      }
                      superClusterHit_index++; 
                      superClusterHit_energy[iCaloCount].push_back(reduceFloat(superClusterHit_energy_,nBits_));
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
                int pfCluster_index_tmp=0;
                if(savePFCluster_){                         
                   for(const auto& iPFCluster : *(pfClusters.product())){
                       reco::PFCluster pfCl(iPFCluster);
                       const std::vector<std::pair<DetId,float> > &hitsAndFractions = pfCl.hitsAndFractions();
                       for(unsigned int i = 0; i < hitsAndFractions.size(); i++){
                           if(hitsAndFractions[i].first.rawId() == id.rawId()){       
                              if(!useEnergyRegression_) pfClusterHit_energy_ = recHit_energy_*hitsAndFractions[i].second;
                              else pfClusterHit_energy_ = (iPFCluster.correctedEnergy()/iPFCluster.energy())*recHit_energy_*hitsAndFractions[i].second;
                              pfCluster_index = pfCluster_index_tmp;
                              break;
                           }else{
                              //for unmatched PFClusterHit 
                              float clusterHit_noCaloPart_energy_ = -1.;
                              float energy_tmp_ = -1;
                              if(hitsAndFractions[i].first.subdetId()==EcalBarrel) energy_tmp_ = (*(recHitsEB.product())->find(hitsAndFractions[i].first)).energy();
                              else if(hitsAndFractions[i].first.subdetId()==EcalEndcap) energy_tmp_ = (*(recHitsEE.product())->find(hitsAndFractions[i].first)).energy();   
                              if(useEnergyRegression_) clusterHit_noCaloPart_energy_ = energy_tmp_*hitsAndFractions[i].second;
                              else clusterHit_noCaloPart_energy_ = (iPFCluster.correctedEnergy()/iPFCluster.energy())*energy_tmp_*hitsAndFractions[i].second;
                              unmatchedHit_energy_Cl_[pfCluster_index_tmp][hitsAndFractions[i].first] = clusterHit_noCaloPart_energy_;
                           }            
                       }
                       pfCluster_index_tmp++;  
                   }  
                   pfClusterHit_index++;  
                   pfClusterHit_energy[iCaloCount].push_back(reduceFloat(pfClusterHit_energy_,nBits_));
                   if(!saveSimhits_ && !saveRechits_ && !savePFRechits_){  
                      pfClusterHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_));
                      pfClusterHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_)); 
                      pfClusterHit_ieta[iCaloCount].push_back(ieta);
                      pfClusterHit_iphi[iCaloCount].push_back(iphi);
                      pfClusterHit_iz[iCaloCount].push_back(iz);  
                   } 
                } 

                //map hits to PFCluster and to SuperCluster 
                if(saveSimhits_){
                   if(savePFCluster_) map_simHit_pfCluster[iCaloCount].insert(pair<int,int>(simHit_index,pfCluster_index));
                   if(saveSuperCluster_) map_simHit_superCluster[iCaloCount].insert(pair<int,int>(simHit_index,superCluster_index));
                }   
                if(!saveSimhits_ && saveRechits_){
                   if(savePFCluster_) map_recHit_pfCluster[iCaloCount].insert(pair<int,int>(recHit_index,pfCluster_index));
                   if(saveSuperCluster_) map_recHit_superCluster[iCaloCount].insert(pair<int,int>(recHit_index,superCluster_index));
                }    
                if(!saveSimhits_ && !saveRechits_ && savePFRechits_){
                   if(savePFCluster_) map_pfRecHit_pfCluster[iCaloCount].insert(pair<int,int>(pfRecHit_index,pfCluster_index));
                   if(saveSuperCluster_) map_pfRecHit_superCluster[iCaloCount].insert(pair<int,int>(pfRecHit_index,superCluster_index));
                }   
                if(!saveSimhits_ && !saveRechits_ && !savePFRechits_){
                   if(savePFCluster_) map_pfClusterHit_pfCluster[iCaloCount].insert(pair<int,int>(pfClusterHit_index,pfCluster_index));
                   if(saveSuperCluster_) map_superClusterHit_superCluster[iCaloCount].insert(pair<int,int>(superClusterHit_index,superCluster_index));
                }     
            }  
       }
       
       caloParticle_simEnergy.push_back(reduceFloat(calo_simEnergy,nBits_));

       iCaloCount++; 
   }
   
   //Save PFClusters
   if(savePFCluster_){
      int iPFCl=0;
      std::cout << "PFClusters size     : " << (pfClusters.product())->size() << std::endl;
      for(const auto& iPFCluster : *(pfClusters.product())){      
          pfCluster_energy.push_back(reduceFloat(iPFCluster.energy(),nBits_));
          pfCluster_eta.push_back(reduceFloat(iPFCluster.eta(),nBits_));
          pfCluster_phi.push_back(reduceFloat(iPFCluster.phi(),nBits_));
          for(auto const& iDetId : unMatchedHits(&unmatchedHit_energy_Cl_,(int)(pfClusters.product())->size(),nCaloParticles)[iPFCl]){
              cell = geometry->getPosition(iDetId.first);
              pfClusterHit_noCaloPart_energy[iPFCl].push_back(reduceFloat(iDetId.second,nBits_));
              pfClusterHit_noCaloPart_eta[iPFCl].push_back(reduceFloat(cell.eta(),nBits_));
              pfClusterHit_noCaloPart_phi[iPFCl].push_back(reduceFloat(cell.phi(),nBits_));
              if(iDetId.first.subdetId()==EcalBarrel){ 
                 EBDetId eb_id(iDetId.first);  
                 pfClusterHit_noCaloPart_ieta[iPFCl].push_back(eb_id.ieta());
                 pfClusterHit_noCaloPart_iphi[iPFCl].push_back(eb_id.iphi());
                 pfClusterHit_noCaloPart_iz[iPFCl].push_back(0); 
              }else if(iDetId.first.subdetId()==EcalEndcap){  
                 EEDetId ee_id(iDetId.first);  
                 pfClusterHit_noCaloPart_ieta[iPFCl].push_back(ee_id.ix());
                 pfClusterHit_noCaloPart_iphi[iPFCl].push_back(ee_id.iy());
                 pfClusterHit_noCaloPart_iz[iPFCl].push_back(ee_id.zside()); 
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
          for(auto const& iDetId : unMatchedHits(&unmatchedHit_energy_SCEB_,(int)(superClusterEB.product())->size(),nCaloParticles)[iSC]){
              cell = geometry->getPosition(iDetId.first);
              EBDetId eb_id(iDetId.first);  
              superClusterHit_noCaloPart_energy[iSC].push_back(reduceFloat(iDetId.second,nBits_));
              superClusterHit_noCaloPart_eta[iSC].push_back(reduceFloat(cell.eta(),nBits_));
              superClusterHit_noCaloPart_phi[iSC].push_back(reduceFloat(cell.phi(),nBits_));
              superClusterHit_noCaloPart_ieta[iSC].push_back(eb_id.ieta());
              superClusterHit_noCaloPart_iphi[iSC].push_back(eb_id.iphi());
              superClusterHit_noCaloPart_iz[iSC].push_back(0); 
          }       
          iSC++;  
       } 
       iSC=0;
       std::cout << "SuperClustersEE size: " << (superClusterEE.product())->size() << std::endl;
       for(const auto& iSuperCluster : *(superClusterEE.product())){    
           superCluster_energy.push_back(reduceFloat(iSuperCluster.energy(),nBits_));
           superCluster_eta.push_back(reduceFloat(iSuperCluster.eta(),nBits_));
           superCluster_phi.push_back(reduceFloat(iSuperCluster.phi(),nBits_));
           for(auto const& iDetId : unMatchedHits(&unmatchedHit_energy_SCEE_,(int)(superClusterEE.product())->size(),nCaloParticles)[iSC]){
              cell = geometry->getPosition(iDetId.first);
              EEDetId ee_id(iDetId.first);  
              superClusterHit_noCaloPart_energy[iSC].push_back(reduceFloat(iDetId.second,nBits_));
              superClusterHit_noCaloPart_eta[iSC].push_back(reduceFloat(cell.eta(),nBits_));
              superClusterHit_noCaloPart_phi[iSC].push_back(reduceFloat(cell.phi(),nBits_));
              superClusterHit_noCaloPart_ieta[iSC].push_back(ee_id.ix());
              superClusterHit_noCaloPart_iphi[iSC].push_back(ee_id.iy());
              superClusterHit_noCaloPart_iz[iSC].push_back(ee_id.zside()); 
          }       
          iSC++;     
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
int RecoSimDumper::nSkimmedCaloParticles(edm::Handle<std::vector<CaloParticle> > caloParticles, std::vector<int>* genID_)
{
    int nCaloParticles=0;
    for(const auto& iCalo : *(caloParticles.product()))
    {
       bool isGoodParticle = false; 
       for(unsigned int id=0; id<genID_->size(); id++)
           //std::cout << "i=" << id << " genID=" <<genID_->at(id) << " iCalo.pdgId()=" << iCalo.pdgId() << std::endl;
           // M.G.
           if(iCalo.pdgId()==genID_->at(id) || genID_->at(id)==0) isGoodParticle=true;
      
       if(isGoodParticle) nCaloParticles++;  
    } 

    return nCaloParticles;
}

std::map<int,std::map<DetId,float> > RecoSimDumper::unMatchedHits(std::map<int,std::map<DetId,float> >* unmatchedHit_energy_, int nClusters, int nCaloParts)
{
    std::map<int,std::map<DetId,float> > unMatchedHits_tmp;
    
    cry_IDs_.clear();
    cry_IDs_.resize(nClusters); 
    int iCluster=0;
    for(auto const& iCl : *unmatchedHit_energy_){   
        for(auto const& iDetId : iCl.second)   
            cry_IDs_[iCluster].push_back(iDetId.first.rawId()); 
            
        std::vector<uint32_t> cry_IDs_tmp = cry_IDs_[iCluster];
        std::sort( cry_IDs_[iCluster].begin(), cry_IDs_[iCluster].end() );
        cry_IDs_[iCluster].erase( std::unique( cry_IDs_[iCluster].begin(), cry_IDs_[iCluster].end() ), cry_IDs_[iCluster].end() );

        for(unsigned int iCry =0; iCry<cry_IDs_[iCluster].size(); iCry++){
            if(std::count(cry_IDs_tmp.begin(), cry_IDs_tmp.end(), cry_IDs_[iCluster].at(iCry))<nCaloParts) continue;
            unMatchedHits_tmp[iCluster][DetId(cry_IDs_[iCluster].at(iCry))] = (*unmatchedHit_energy_)[iCluster][DetId(cry_IDs_[iCluster].at(iCry))];  
        }

        iCluster++; 
    }
    
    return unMatchedHits_tmp;
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

